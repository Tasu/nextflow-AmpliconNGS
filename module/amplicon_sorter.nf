/*
 * AMPLICON_SORTER Workflow Components
 * 1. AS_PRE_STATS: Count reads before clustering (using SeqKit)
 * 2. AMPLICON_SORTER: Core clustering (using Local SIF + Python Script)
 * 3. AS_DEDUPLICATE: Merge and remove duplicate OTUs (using SeqKit)
 */

/*
 * Step 1: Pre-processing stats to count total reads used
 */
process AS_PRE_STATS {
    tag "${sample_id}"
    label 'process_low'
    container "${params.container_images.seqkit}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("read_count_actual.txt"), emit: count
    path "versions_as_pre.yml",                          emit: versions

    script:
    """
    # Count actual reads (denominator for downstream analysis) with a safe zero fallback.
    if [[ ! -s "${reads}" ]]; then
        echo 0 > read_count_actual.txt
    else
        read_count=\$(seqkit stats -T "${reads}" 2>/dev/null | awk 'NR==2 {print \$4}')
        if [[ -z "\${read_count}" ]]; then
            echo 0 > read_count_actual.txt
        else
            echo "\${read_count}" > read_count_actual.txt
        fi
    fi
    
    echo "SeqKit (AS_PRE): \$(seqkit version | awk '{print \$2}')" > versions_as_pre.yml
    """
}

/*
 * Step 2: Core clustering process
 */
process AMPLICON_SORTER {
    tag "${sample_id}"
    label 'process_high'
    // Local SIF container with edlib, biopython, matplotlib, numpy
    container "${params.as_container_path}"

    input:
    tuple val(sample_id), path(reads), val(min_len), val(max_len)
    val max_reads
    path as_script  // Path to local amplicon_sorter.py (2025-10-09)

    output:
    tuple val(sample_id), path("out_dir/*_consensussequences.fasta"), emit: raw_consensus, optional: true
    path "out_dir/results.txt",                                      emit: summary,       optional: true
    path "out_dir/*",                                                emit: all_outputs
    path "versions_as_main.yml",                                     emit: versions

    script:
    """
    # Skip clustering when there are zero reads and emit empty placeholders.
    if [[ "${reads}" == *.gz ]]; then
        input_reads=\$(gzip -cd "${reads}" | awk 'END{print int(NR/4)}')
    else
        input_reads=\$(awk 'END{print int(NR/4)}' "${reads}")
    fi

    mkdir -p out_dir
    if [[ "\${input_reads}" -eq 0 ]]; then
        touch out_dir/${sample_id}_empty_consensussequences.fasta
        printf "No reads available for clustering (sample=%s)\n" "${sample_id}" > out_dir/results.txt
    else
        # Execute the local clustering script
        python3 ${as_script} \\
            -i ${reads} \\
            -n ${max_reads} \\
            -min ${min_len} \\
            -max ${max_len} \\
            -o out_dir
    fi

    echo "Amplicon_Sorter_Script: 2025-10-09 (Local)" > versions_as_main.yml
    """
}

/*
 * Step 3: Post-processing to merge and deduplicate OTUs
 */
process AS_DEDUPLICATE {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/03_amplicon_sorter/${sample_id}", mode: 'copy'
    container "${params.container_images.seqkit}"

    input:
    tuple val(sample_id), path(raw_fastas)

    output:
    tuple val(sample_id), path("${sample_id}_clustered_consensus.fasta"), emit: consensus
    path "versions_as_post.yml",                                         emit: versions

    script:
    """
    # Merge all cluster consensus files and remove exact duplicates.
    # If all inputs are empty, create an empty consensus file for downstream zero-count handling.
    has_content=0
    for f in ${raw_fastas}; do
        if [[ -s "\$f" ]]; then
            has_content=1
            break
        fi
    done

    if [[ "\${has_content}" -eq 1 ]]; then
        cat ${raw_fastas} > merged_tmp.fasta
        seqkit rmdup -s merged_tmp.fasta -o ${sample_id}_clustered_consensus.fasta
    else
        touch ${sample_id}_clustered_consensus.fasta
    fi

    echo "SeqKit (AS_POST): \$(seqkit version | awk '{print \$2}')" > versions_as_post.yml
    """
}