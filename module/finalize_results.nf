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
    path "01_Final_Analysis/*"
    path "02_Sequences/*"
    path "03_Annotation/*"

    script:
    """
    # Create directory structure
    mkdir -p 01_Final_Analysis 02_Sequences 03_Annotation

    # 1. Final OTU Count Matrix
    cp ${count_matrix} 01_Final_Analysis/otu_count_matrix.tsv

    # 2. Sequences
    cp ${unique_otu_fasta} 02_Sequences/integrated_unique_otus.fasta
    # Merge all sample consensus sequences
    cat ${all_consensus_fasta} > 02_Sequences/all_samples_consensus.fasta

    # 3. Annotation Results
    cp ${blast_results} 03_Annotation/blast_annotation.tsv
    """
}