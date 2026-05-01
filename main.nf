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

def isHttpReachable(String url, int timeoutMs = 8000) {
    try {
        def conn = new URL(url).openConnection() as HttpURLConnection
        conn.setRequestMethod('HEAD')
        conn.setConnectTimeout(timeoutMs)
        conn.setReadTimeout(timeoutMs)
        conn.setInstanceFollowRedirects(true)
        conn.connect()
        int code = conn.responseCode
        conn.disconnect()
        return code >= 200 && code < 400
    } catch (all) {
        return false
    }
}

def fetchDepotContainerIndex(String indexUrl, int timeoutMs = 8000) {
    try {
        def conn = new URL(indexUrl).openConnection() as HttpURLConnection
        conn.setRequestMethod('GET')
        conn.setConnectTimeout(timeoutMs)
        conn.setReadTimeout(timeoutMs)
        conn.setInstanceFollowRedirects(true)
        def html = conn.inputStream.getText('UTF-8')
        conn.disconnect()
        def matcher = (html =~ /href="([^"]+)"/)
        return matcher.collect { it[1] }
    } catch (all) {
        return []
    }
}

def suggestContainerTags(String url, List<String> indexEntries, int limit = 5) {
    if (!url?.startsWith('https://depot.galaxyproject.org/singularity/') || !indexEntries) {
        return []
    }

    def fileName = url.tokenize('/')[-1]
    def imageName = fileName.contains(':') ? fileName.split(':', 2)[0] : fileName
    if (!imageName) {
        return []
    }

    return indexEntries
        .findAll { it.startsWith("${imageName}:") }
        .sort()
        .reverse()
        .take(limit)
}

def commandExists(String command, int timeoutMs = 3000) {
    try {
        def proc = new ProcessBuilder(command, '--version').redirectErrorStream(true).start()
        def finished = proc.waitFor(timeoutMs as long, java.util.concurrent.TimeUnit.MILLISECONDS)
        if (!finished) {
            proc.destroyForcibly()
            return false
        }
        return proc.exitValue() == 0
    } catch (all) {
        return false
    }
}

def isDockerRefReachable(String dockerRef, int timeoutMs = 8000) {
    try {
        def proc = new ProcessBuilder('skopeo', 'inspect', '--no-tags', dockerRef).redirectErrorStream(true).start()
        def finished = proc.waitFor(timeoutMs as long, java.util.concurrent.TimeUnit.MILLISECONDS)
        if (!finished) {
            proc.destroyForcibly()
            return false
        }
        return proc.exitValue() == 0
    } catch (all) {
        return false
    }
}

def runContainerPreflightChecks(params) {
    def imageMap = (params.container_images ?: [:]) as Map
    def timeoutMs = (params.container_preflight_timeout_ms ?: 8000) as int
    def strictMode = (params.container_preflight_strict ?: false) as boolean
    def suggestTags = (params.container_preflight_suggest_tags ?: true) as boolean
    def checkDockerRefs = (params.container_preflight_check_docker ?: true) as boolean
    def suggestionLimit = (params.container_suggestion_limit ?: 5) as int
    def indexUrl = (params.container_registry_index_url ?: 'https://depot.galaxyproject.org/singularity/').toString()

    if (!imageMap) {
        log.warn "[container-preflight] No entries found in params.container_images. Skipping checks."
        return
    }

    def uniqueUrls = imageMap.values().findAll { it != null && it.toString().trim() }.collect { it.toString() }.unique()
    def dockerRefs = uniqueUrls.findAll { it.startsWith('docker://') }
    def depotIndex = suggestTags ? fetchDepotContainerIndex(indexUrl, timeoutMs) : []
    def failures = []

    def skopeoAvailable = true
    if (checkDockerRefs && dockerRefs) {
        skopeoAvailable = commandExists('skopeo')
        if (!skopeoAvailable) {
            def msg = '[container-preflight] docker:// reference check requires skopeo, but it was not found on PATH.'
            failures << msg
            log.warn msg
        }
    }

    uniqueUrls.each { url ->
        if (url.startsWith('docker://')) {
            if (!checkDockerRefs) {
                log.info "[container-preflight] Skip docker:// check by config: ${url}"
                return
            }
            if (!skopeoAvailable) {
                return
            }
            if (!isDockerRefReachable(url, timeoutMs)) {
                def msg = "[container-preflight] Unreachable docker reference: ${url}"
                failures << msg
                log.warn msg
            } else {
                log.info "[container-preflight] OK (docker registry): ${url}"
            }
            return
        }

        if (!(url.startsWith('http://') || url.startsWith('https://'))) {
            log.info "[container-preflight] Skip connectivity check for non-HTTP container ref: ${url}"
            return
        }

        if (!isHttpReachable(url, timeoutMs)) {
            def msg = "[container-preflight] Unreachable container URL: ${url}"
            if (suggestTags) {
                def candidates = suggestContainerTags(url, depotIndex, suggestionLimit)
                if (candidates) {
                    msg += " | candidate tags: ${candidates.join(', ')}"
                }
            }
            failures << msg
            log.warn msg
        } else {
            log.info "[container-preflight] OK: ${url}"
        }
    }

    if (strictMode && failures) {
        throw new IllegalStateException("Container preflight failed for ${failures.size()} URL(s). Set --container_preflight_strict false to continue with warnings.")
    }
}

workflow {

    def preflightOnly = (params.preflight_only ?: false) as boolean
    def runPreflight = (params.container_preflight_check ?: false) as boolean || preflightOnly

    if (runPreflight) {
        runContainerPreflightChecks(params)
    }

    if (preflightOnly) {
        log.info "[container-preflight] Preflight-only mode enabled. Exiting workflow after checks."
        return
    }

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

    def workflow_info = [
        nextflow_version: workflow.nextflow.version.toString(),
        command_line: workflow.commandLine?.toString() ?: "N/A"
    ]

    GENERATE_PROVENANCE(ch_versions, workflow_info)
}