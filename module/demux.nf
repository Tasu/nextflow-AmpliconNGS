/*
 * DEMULTIPLEX Module: Separated into Cutadapt and Python(Biopython)
 */

process CUTADAPT_MARK {
    tag "${sample_id}"
    container 'https://depot.galaxyproject.org/singularity/cutadapt:5.0'

    input:
    tuple val(sample_id), path(fastq), val(min_len), val(max_len), val(f_idx), val(f_prm), val(r_idx), val(r_prm)

    output:
    tuple val(sample_id), path("marked.fastq"), val(min_len), val(max_len), val(f_idx), val(r_idx), emit: marked_data
    path "versions_cutadapt.yml", emit: versions

    script:
    """
    # Orientation search and lowercase marking
    cutadapt -j ${task.cpus} -g "${f_prm}" --discard-untrimmed ${fastq} | \
    cutadapt -j ${task.cpus} -a "\$(echo ${r_prm} | rev | tr ATCG TAGC)" --discard-untrimmed - --action=lowercase -o fwd_marked.fastq
    
    # (Reverse search and merge logic omitted for brevity, but same principle applies)
    cat fwd_marked.fastq > marked.fastq
    echo "Cutadapt: \$(cutadapt --version)" > versions_cutadapt.yml
    """
}

process BIOPYTHON_EXTRACT {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_demux/${sample_id}", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'

    input:
    tuple val(sample_id), path(marked_fastq), val(min_len), val(max_len), val(f_idx), val(r_idx)

    output:
    tuple val(sample_id), path("${sample_id}_final.fastq.gz"), emit: reads
    path "versions_biopython.yml", emit: versions

    script:
    def f_idx_len = f_idx ? f_idx.length() : 0
    """
    python3 -c "
import sys, re
from Bio import SeqIO
# Logic to extract uppercase insert between lowercase primers
# Apply min_len / max_len filtering here
    " > ${sample_id}_final.fastq
    gzip ${sample_id}_final.fastq
    echo "BioPython: 1.79" > versions_biopython.yml
    """
}