nextflow.enable.dsl = 2

/*
 * 1. Import Refactored Modules
 */
include { PORECHOP_TRIM as PREPROCESS_ADAPTERS; SEQKIT_CLEAN } from './module/preprocess.nf'
include { CUTADAPT_MARK; BIOPYTHON_EXTRACT                 } from './module/demux.nf'
include { KRAKEN2_CLASSIFY; KRAKENTOOLS_EXTRACT             } from './module/kraken2_filter.nf'
include { AS_PRE_STATS; AMPLICON_SORTER; AS_DEDUPLICATE      } from './module/amplicon_sorter.nf'
include { OTU_COUNT_TABLE                                   } from './module/otu_count.nf'
include { OTU_MERGE                                         } from './module/otu_merge.nf'
include { BLAST_ANNOTATE                                    } from './module/blast_annotate.nf'
include { TAXONKIT_LINEAGE; BIOM_GENERATE                   } from './module/create_biom.nf'
include { SUMMARY_REPORT                                    } from './module/summary_report.nf'
include { FINALIZE_RESULTS                                  } from './module/finalize_results.nf'
include { GENERATE_PROVENANCE                               } from './module/generate_provenance.nf'

workflow {

    // --- 2. Input Channel Creation ---
    
    // Parse sample sheet and create a metadata map for each sample row
    ch_sample_sheet = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id:  row.sample ?: row.sample_id,
                fastq_dir:  row.fastq_dir ?: row.data_dir ?: row.fastq_passDir ?: row.path,
                min_len:    row.min_len ? row.min_len.toInteger() : (row.min ? row.min.toInteger() : 0),
                max_len:    row.max_len ? row.max_len.toInteger() : (row.max ? row.max.toInteger() : 0),
                f_idx:      row.fwd_index ?: row.iF ?: "",
                f_prm:      row.fwd_primer ?: row.FwPrimer ?: "",
                r_idx:      row.rev_index ?: row.iR ?: "",
                r_prm:      row.rev_primer ?: row.RvPrimer ?: ""
            ]
            return meta
        }

    // Group by fastq_dir so PREPROCESS runs once per directory.
    // Deduplicate by full path to avoid staging collisions when sample-sheet rows
    // reference the same directory and FASTQ basenames are repeated.
    ch_for_preprocess = ch_sample_sheet
        .map { meta -> [ meta.fastq_dir, file("${meta.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}") ] }
        .groupTuple()
        .map { dir, files ->
            def dir_id = dir.split('/')[-1]
            def unique_fastq_files = files.flatten().unique { it.toString() }
            [ dir_id, unique_fastq_files ]
        }

    // --- 3. Parallel Processing Phase (Per-Directory / Per-Sample) ---

    // [A] Pre-processing (Adapter trimming via Porechop)
    PREPROCESS_ADAPTERS(ch_for_preprocess)

    // [B] Length Filtering (SeqKit Clean) - Remove reads marked by Porechop
    SEQKIT_CLEAN(
        PREPROCESS_ADAPTERS.out.round1_fastq,
        PREPROCESS_ADAPTERS.out.log_round2
    )

    // [C] Redistribute trimmed reads back to individual samples
    ch_demux_input = ch_sample_sheet
        .map { meta -> [ meta.fastq_dir.split('/')[-1], meta ] }
        .combine(SEQKIT_CLEAN.out.chopped_fastq, by: 0)
        .map { dir_id, meta, fastq ->
            [ meta.sample_id, fastq, meta.min_len, meta.max_len, meta.f_idx, meta.f_prm, meta.r_idx, meta.r_prm ]
        }

    // [D] High-precision Demultiplexing - Cutadapt marking
    CUTADAPT_MARK(ch_demux_input)

    // [E] Extract sequences between primers with Biopython
    BIOPYTHON_EXTRACT(CUTADAPT_MARK.out.marked_data)

    // [F] Host Removal / Target Extraction (Kraken2 Filter Split)
    KRAKEN2_CLASSIFY(
        BIOPYTHON_EXTRACT.out.reads,
        params.kraken2_db
    )
    
    ch_for_extraction = KRAKEN2_CLASSIFY.out.output
        .join(KRAKEN2_CLASSIFY.out.report)
        .join(BIOPYTHON_EXTRACT.out.reads)

    KRAKENTOOLS_EXTRACT(
        ch_for_extraction,
        params.target_taxid
    )

    // [E] Clustering & Variant Identification (Amplicon_Sorter Split)
    // E-1: Count reads before sorting (Denominator for count table)
    AS_PRE_STATS(KRAKENTOOLS_EXTRACT.out.filtered_fastq)

    // E-2: Core Clustering (Requires metadata re-join for min/max lengths)
    ch_as_input = KRAKENTOOLS_EXTRACT.out.filtered_fastq
        .join(ch_sample_sheet.map { m -> [m.sample_id, m.min_len, m.max_len] })
    
    AMPLICON_SORTER(
        ch_as_input,
        params.max_reads,
        file(params.as_script_path)
    )

    // E-3: Post-clustering Deduplication (Using SeqKit)
    AS_DEDUPLICATE(AMPLICON_SORTER.out.raw_consensus)

    // [F] Per-sample Read Counting
    ch_otu_count_input = AS_DEDUPLICATE.out.consensus
        .join(AS_PRE_STATS.out.count)

    OTU_COUNT_TABLE(ch_otu_count_input)

    // --- 4. Global Consolidation Phase ---

    ch_all_consensus = AS_DEDUPLICATE.out.consensus.map { it[1] }.collect()
    ch_all_counts    = OTU_COUNT_TABLE.out.final_table.collect()

    // [G] Merge OTUs across all samples (100% identity)
    OTU_MERGE(ch_all_consensus, ch_all_counts)

    // [H] Batch BLAST Annotation
    BLAST_ANNOTATE(
        OTU_MERGE.out.otu_fasta,
        params.blast_db_dir,
        params.blast_db_name,
        params.blast_type
    )

    // --- 5. Result Finalization & Provenance ---

    // [I] Taxonomic Lineage Extraction from BLAST results (TaxonKit)
    TAXONKIT_LINEAGE(
        BLAST_ANNOTATE.out.blast_results,
        params.tax_db_dir
    )

    // [J] BIOM Format Generation (OTU counts + taxonomic lineage)
    BIOM_GENERATE(
        OTU_MERGE.out.count_matrix,
        TAXONKIT_LINEAGE.out.lineage,
        BLAST_ANNOTATE.out.blast_results
    )

    // [K] Summary Report Generation (Phylum-level aggregation)
    ch_kraken_reports = KRAKEN2_CLASSIFY.out.report.map { it[1] }.collect()
    
    SUMMARY_REPORT(
        ch_kraken_reports,
        OTU_MERGE.out.count_matrix,
        BLAST_ANNOTATE.out.blast_results
    )

    // [L] Organize final results
    FINALIZE_RESULTS(
        OTU_MERGE.out.count_matrix,
        OTU_MERGE.out.otu_fasta,
        ch_all_consensus,
        BLAST_ANNOTATE.out.blast_results
    )

    // Aggregate version info from all utilized containers
    ch_versions = PREPROCESS_ADAPTERS.out.versions
        .mix(SEQKIT_CLEAN.out.versions)
        .mix(CUTADAPT_MARK.out.versions)
        .mix(BIOPYTHON_EXTRACT.out.versions)
        .mix(KRAKEN2_CLASSIFY.out.versions)
        .mix(KRAKENTOOLS_EXTRACT.out.versions)
        .mix(AS_PRE_STATS.out.versions)
        .mix(AMPLICON_SORTER.out.versions)
        .mix(AS_DEDUPLICATE.out.versions)
        .mix(OTU_COUNT_TABLE.out.versions)
        .mix(OTU_MERGE.out.versions)
        .mix(BLAST_ANNOTATE.out.versions)
        .mix(TAXONKIT_LINEAGE.out.versions)
        .mix(BIOM_GENERATE.out.versions)
        .mix(SUMMARY_REPORT.out.versions)
        .collect()

    GENERATE_PROVENANCE(ch_versions, workflow)
}