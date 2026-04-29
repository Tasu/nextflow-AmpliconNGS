/*
 * PREPROCESS Module: Separated into Porechop and SeqKit
 */

process PORECHOP_TRIM {
    tag "${fastq_dir_id}"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/porechop:0.2.4'

    input:
    tuple val(fastq_dir_id), path(fastq_files)

    output:
    tuple val(fastq_dir_id), path("round1.fastq.gz"), emit: round1_fastq
    path "${fastq_dir_id}_porechop_round2.log",       emit: log_round2
    path "versions_porechop.yml",                     emit: versions

    script:
    def n_threads = task.cpus > 4 ? 4 : task.cpus
    """
    zcat -f ${fastq_files} | gzip > merged.fastq.gz

    # Round 1: Trim and discard middle
    porechop -i merged.fastq.gz -o round1.fastq.gz --threads ${n_threads} --extra_end_trim -1 --discard_middle

    # Round 2: Detection only for strict filtering later
    porechop -i round1.fastq.gz -o /dev/null --threads ${n_threads} -v 2 > ${fastq_dir_id}_porechop_round2.log

    echo "Porechop: \$(porechop --version 2>&1)" > versions_porechop.yml
    """
}

process SEQKIT_CLEAN {
    tag "${fastq_dir_id}"
    label 'process_low'
    publishDir "${params.outdir}/00_preprocess/${fastq_dir_id}", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0'

    input:
    tuple val(fastq_dir_id), path(round1_fastq)
    path log_round2

    output:
    tuple val(fastq_dir_id), path("${fastq_dir_id}_chopped.fastq.gz"), emit: chopped_fastq
    path "versions_seqkit.yml", emit: versions

    script:
    """
    grep "trimmed" ${log_round2} | awk '{print \$1}' | sed 's/@//' > discarded_ids.txt || touch discarded_ids.txt

    if [ -s discarded_ids.txt ]; then
        seqkit grep -v -f discarded_ids.txt ${round1_fastq} -o ${fastq_dir_id}_chopped.fastq.gz
    else
        cp ${round1_fastq} ${fastq_dir_id}_chopped.fastq.gz
    fi

    echo "SeqKit: \$(seqkit version | awk '{print \$2}')" > versions_seqkit.yml
    """
}