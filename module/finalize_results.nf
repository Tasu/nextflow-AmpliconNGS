/*
 * FINALIZE_RESULTS:
 * Consolidates all key outputs into a structured directory for the user.
 */
process FINALIZE_RESULTS {
    label 'process_low'
    publishDir "\${params.outdir}", mode: 'copy'

    input:
    path merged_biom         // From CREATE_BIOM
    path merged_report       // From CREATE_BIOM
    path unique_otu_fasta    // From OTU_MERGE
    path all_consensus_fasta // From AMPLICON_SORTER (collected)
    path stats_tables        // From SUMMARY_REPORT (collected)
    path plots               // From SUMMARY_REPORT (collected)

    output:
    path "01_Final_Analysis/*"
    path "02_Sequences/*"
    path "03_Statistics/*"

    script:
    """
    # Create directory structure
    mkdir -p 01_Final_Analysis 02_Sequences 03_Statistics

    # 1. Final Tables & BIOM
    cp \${merged_biom} 01_Final_Analysis/
    cp \${merged_report} 01_Final_Analysis/

    # 2. Sequences
    cp \${unique_otu_fasta} 02_Sequences/
    # 統合コンセンサスを作成（全サンプルのfastaを一つにまとめる）
    cat \${all_consensus_fasta} > 02_Sequences/all_samples_consensus.fasta

    # 3. Statistics & Plots
    cp \${stats_tables} 03_Statistics/
    cp \${plots} 03_Statistics/
    """
}