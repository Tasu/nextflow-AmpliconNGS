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
    // Use Galaxy Depot version-only tag
    container 'depot.galaxyproject.org/singularity/kraken2:2.1.3'

    input:
    tuple val(sample_id), path(fastq)
    path kraken2_db

    output:
    tuple val(sample_id), path("${sample_id}.kraken2.out"),    emit: output
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report
    path "versions_kraken2.yml",                               emit: versions

    script:
    """
    # Skip Kraken2 when the input FASTQ has zero reads.
    if [[ "${fastq}" == *.gz ]]; then
        input_reads=\$(gzip -cd "${fastq}" | awk 'END{print int(NR/4)}')
    else
        input_reads=\$(awk 'END{print int(NR/4)}' "${fastq}")
    fi

    if [[ "\${input_reads}" -eq 0 ]]; then
        : > ${sample_id}.kraken2.out
        printf "0.00\t0\t0\tU\t0\tunclassified\n" > ${sample_id}.kraken2.report
    else
        # Execute Kraken2 classification
        kraken2 \\
            --db ${kraken2_db} \\
            --threads ${task.cpus} \\
            --output ${sample_id}.kraken2.out \\
            --report ${sample_id}.kraken2.report \\
            ${fastq}
    fi

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
    // Use Docker version-only tag
    container 'docker://nanozoo/krakentools:1.2'

    input:
    tuple val(sample_id), path(kraken_out), path(kraken_report), path(fastq)
    val taxid_to_extract // Target TaxID (e.g., from params.target_taxid)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered_fastq
    path "versions_krakentools.yml",                             emit: versions

    script:
    """
    # Count classified reads directly from kraken2 output.
    classified_reads=\$(awk 'BEGIN{c=0} /^C\t/{c++} END{print c+0}' ${kraken_out})

    if [[ "\${classified_reads}" -eq 0 ]]; then
        : > ${sample_id}_filtered.fastq
    else
        # Extract specific reads matching the TaxID and its children
        extract_kraken_reads.py \\
            -k ${kraken_out} \\
            -r ${kraken_report} \\
            -s ${fastq} \\
            -t ${taxid_to_extract} \\
            --include-children \\
            -o ${sample_id}_filtered.fastq
    fi

    # Compress the output FASTQ (empty file is valid and represents zero reads)
    gzip ${sample_id}_filtered.fastq

    # Record KrakenTools version
    echo "KrakenTools (extract_kraken_reads.py): 1.2" > versions_krakentools.yml
    """
}