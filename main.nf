nextflow.enable.dsl = 2

/*
 * 1. Import Refactored Modules
 */
include { PORECHOP_TRIM as PREPROCESS_ADAPTERS; SEQKIT_CLEAN } from './modules/preprocess.nf'
include { CUTADAPT_MARK; BIOPYTHON_EXTRACT                 } from './modules/demux.nf'
include { KRAKEN2_CLASSIFY; KRAKENTOOLS_EXTRACT             } from './modules/kraken2_filter.nf'
include { AS_PRE_STATS; AMPLICON_SORTER; AS_DEDUPLICATE      } from './modules/amplicon_sorter.nf'
include { OTU_COUNT_TABLE                                   } from './modules/otu_count.nf'
include { OTU_MERGE                                         } from './modules/otu_merge.nf'
include { BLAST_ANNOTATE                                    } from './modules/blast_annotate.nf'
include { FINALIZE_RESULTS                                  } from './modules/finalize_results.nf'
include { GENERATE_PROVENANCE                               } from './modules/generate_provenance.nf'

workflow {

    // --- 2. Input Channel Creation ---
    
    // Parse sampleSheet.csv and create a metadata map for each sample
    ch_sample_sheet = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id:  row.sample,     // Updated to match your CSV column 'sample'
                fastq_dir:  row.fastq_dir,
                min_len:    row.min_len.toInteger(),
                max_len:    row.max_len.toInteger(),
                f_idx:      row.fwd_index ?: "",
                f_prm:      row.fwd_primer,
                r_idx:      row.rev_index ?: "",
                r_prm:      row.rev_primer
            ]
            return meta
        }

    // Grouping by fastq_dir to run PREPROCESS only once per directory
    ch_for_preprocess = ch_sample_sheet
        .map { meta -> [ meta.fastq_dir, file("${meta.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}") ] }
        .groupTuple()
        .map { dir, files -> [ dir.split('/')[-1], files.flatten() ] }

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
    OTU_COUNT_TABLE(
        AS_DEDUPLICATE.out.consensus,
        AS_PRE_STATS.out.count
    )

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

    // Organize final results
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
        .collect()

    GENERATE_PROVENANCE(ch_versions, workflow)
}