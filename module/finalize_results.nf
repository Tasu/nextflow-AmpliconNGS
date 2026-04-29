/*
 * FINALIZE_RESULTS:
 * Consolidates all key outputs into a structured directory for the user.
 */
process FINALIZE_RESULTS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path count_matrix        // From OTU_MERGE
    path unique_otu_fasta    // From OTU_MERGE
    path all_consensus_fasta // From AS_DEDUPLICATE (collected)
    path blast_results       // From BLAST_ANNOTATE

    output:
    path "06_final_results/*"
    path "06_final_results/sequences/*"
    path "06_final_results/blast/*"

    script:
    """
    # Create directory structure
    mkdir -p 06_final_results/sequences 06_final_results/blast

    # 1. Final OTU Count Matrix
    cp ${count_matrix} 06_final_results/otu_count_matrix.tsv

    # 2. Sequences
    cp ${unique_otu_fasta} 06_final_results/sequences/integrated_unique_otus.fasta
    # Merge all sample consensus sequences
    cat ${all_consensus_fasta} > 06_final_results/sequences/all_samples_consensus.fasta

    # 3. Annotation Results
    cp ${blast_results} 06_final_results/blast/blast_annotation.tsv
    """
}