/*
 * KRAKEN2_FILTER Module
 * * Description:
 * This module is split into two processes to handle different container environments.
 * 1. KRAKEN2_CLASSIFY: Assigns taxonomy to consensus sequences.
 * 2. KRAKENTOOLS_EXTRACT: Extracts specific reads (e.g., target parasites) based on taxonomy.
 */

/*
 * Process 1: Taxonomic Classification
 */
process KRAKEN2_CLASSIFY {
    tag "${sample_id}"
    label 'process_high'
    // Using verified stable URI from Galaxy Depot
    container 'https://depot.galaxyproject.org/singularity/kraken2:2.1.3--pl5321h077b44d_4'

    input:
    tuple val(sample_id), path(fastq)
    path kraken2_db

    output:
    tuple val(sample_id), path("${sample_id}.kraken2.out"),    emit: output
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report
    path "versions_kraken2.yml",                               emit: versions

    script:
    """
    # Execute Kraken2 classification
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample_id}.kraken2.out \\
        --report ${sample_id}.kraken2.report \\
        ${fastq}

    # Capture Kraken2 version info
    echo "Kraken2: \$(kraken2 --version | head -n 1 | awk '{print \$3}')" > versions_kraken2.yml
    """
}

/*
 * Process 2: Read Extraction by Taxonomy ID
 */
process KRAKENTOOLS_EXTRACT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/04_kraken_filter/${sample_id}", mode: 'copy'
    // Using verified stable URI from Docker Hub (Nanozoo)
    container 'docker://nanozoo/krakentools:1.2--13d5ba5'

    input:
    tuple val(sample_id), path(kraken_out), path(kraken_report), path(fastq)
    val taxid_to_extract // Target TaxID (e.g., from params.target_taxid)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered_fastq
    path "versions_krakentools.yml",                             emit: versions

    script:
    """
    # Extract specific reads matching the TaxID and its children
    extract_kraken_reads.py \\
        -k ${kraken_out} \\
        -r ${kraken_report} \\
        -s ${fastq} \\
        -t ${taxid_to_extract} \\
        --include-children \\
        -o ${sample_id}_filtered.fastq

    # Compress the output fastq
    gzip ${sample_id}_filtered.fastq

    # Record KrakenTools version
    echo "KrakenTools (extract_kraken_reads.py): 1.2" > versions_krakentools.yml
    """
}