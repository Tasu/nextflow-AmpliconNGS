/*
 * DEMULTIPLEX Module: Separated into Cutadapt and Python(Biopython)
 */

process CUTADAPT_MARK {
    tag "${sample_id}"
    container "${params.container_images.cutadapt}"

    input:
    tuple val(sample_id), path(fastq), val(min_len), val(max_len), val(f_idx), val(f_prm), val(r_idx), val(r_prm)

    output:
    tuple val(sample_id), path("marked.fastq"), val(min_len), val(max_len), val(f_idx), val(r_idx), emit: marked_data
    path "versions_cutadapt.yml", emit: versions

    script:
    // Reverse complement the primers for the "reverse" orientation search
    """
    # 1. Forward orientation: F-primer at Start, R-primer (RC) at End
    # We use --action=lowercase to mask primers
    cutadapt -j ${task.cpus} -g "${f_prm}" --discard-untrimmed ${fastq} | \
    cutadapt -j ${task.cpus} -a "\$(echo ${r_prm} | rev | tr ATCG TAGC)" --discard-untrimmed --action=lowercase -o fwd_marked.fastq -

    # 2. Reverse orientation: R-primer at Start, F-primer (RC) at End
    cutadapt -j ${task.cpus} -g "${r_prm}" --discard-untrimmed ${fastq} | \
    cutadapt -j ${task.cpus} -a "\$(echo ${f_prm} | rev | tr ATCG TAGC)" --discard-untrimmed --action=lowercase -o rev_marked.fastq -

    # Combine results
    cat fwd_marked.fastq rev_marked.fastq > marked.fastq

    cat <<-END_VERSIONS > versions_cutadapt.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}


process BIOPYTHON_EXTRACT {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_demux/${sample_id}", mode: 'copy'
    container "${params.container_images.biopython}"

    input:
    tuple val(sample_id), path(marked_fastq), val(min_len), val(max_len), val(f_idx), val(r_idx)

    output:
    tuple val(sample_id), path("${sample_id}_final.fastq.gz"), emit: reads
    path "versions_biopython.yml", emit: versions

    script:
    // Normalize null/empty values on the Groovy side
    def f_idx_val = f_idx ?: ""
    def r_idx_val = r_idx ?: ""
    """
    python3 -c "
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

# 1. Regex to identify uppercase insert sequence flanked by lowercase primers
# group(1): F-primer, group(2): Insert, group(3): R-primer
pattern = re.compile(r'([a-z]+)([A-Z]+)([a-z]+)')

f_target = '${f_idx_val}'
raw_r_idx = '${r_idx_val}'
# r_idx is provided in 5'->3', so convert it for comparison against read sequence (RC)
r_target = str(Seq(raw_r_idx).reverse_complement()) if raw_r_idx else ''

f_len = len(f_target)
r_len = len(r_target)

with open('${sample_id}_final.fastq', 'w') as out_handle:
    for record in SeqIO.parse('${marked_fastq}', 'fastq'):
        seq_str = str(record.seq)
        match = pattern.search(seq_str)
        
        if match:
            # Get boundary coordinates for each region (0-based)
            f_prm_start = match.start(1)
            insert_start, insert_end = match.span(2)
            r_prm_end = match.end(3)
            
            is_valid = True
            
            # --- Index Validation Logic ---
            # Forward index: check exact match immediately before the primer
            if f_len > 0:
                if f_prm_start < f_len:
                    is_valid = False
                else:
                    detected_f = seq_str[f_prm_start - f_len : f_prm_start]
                    if detected_f != f_target:
                        is_valid = False
            
            # Reverse index: check exact match (RC) immediately after the primer
            if is_valid and r_len > 0:
                if (len(seq_str) - r_prm_end) < r_len:
                    is_valid = False
                else:
                    detected_r = seq_str[r_prm_end : r_prm_end + r_len]
                    if detected_r != r_target:
                        is_valid = False
            
            # --- Final Extraction ---
            if is_valid:
                # Slice only the insert region (quality scores are preserved automatically)
                insert_record = record[insert_start:insert_end]
                # Apply length filtering before writing output
                if ${min_len} <= len(insert_record.seq) <= ${max_len}:
                    SeqIO.write(insert_record, out_handle, 'fastq')

    "
    gzip -f ${sample_id}_final.fastq

    cat <<-END_VERSIONS > versions_biopython.yml
    '${task.process}':
        python: \$(python3 --version | cut -d ' ' -f 2)
        biopython: 1.79
    END_VERSIONS
    """
}