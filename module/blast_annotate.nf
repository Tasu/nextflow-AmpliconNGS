/*
 * BLAST_ANNOTATE: 
 * Perform batch BLAST search for the integrated unique OTUs.
 * Minimizes redundant computations by annotating each unique sequence only once.
 * Software: BLAST+ v2.16.0
 */
process BLAST_ANNOTATE {
    label 'process_high'
    publishDir "${params.outdir}/05_annotation", mode: 'copy'

    // Updated to the verified stable URI
    container 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'

    input:
    path unique_otus_fasta  // Integrated FASTA file from OTU_MERGE
    path db_dir             // Directory containing BLAST database (mounted via nextflow.config or --blast_db_dir)
    val db_name             // Name of the database (e.g., 'nt')
    val blast_type          // BLAST algorithm (e.g., 'blastn')

    output:
    path "all_otus_blast_results.tsv", emit: blast_results
    path "versions_blast.yml",         emit: versions

    script:
    """
    # If the input FASTA is empty (all samples had zero OTUs), skip BLAST and
    # emit an empty result file to allow downstream processes to continue.
    if [[ ! -s "${unique_otus_fasta}" ]]; then
        echo "WARNING: Input FASTA is empty — skipping BLAST annotation." >&2
        touch all_otus_blast_results.tsv
    else
        # Execute batch BLAST search
        # Params: min 70% query coverage, max 1 target, evalue 1e-10
        ${blast_type} \\
            -query ${unique_otus_fasta} \\
            -db ${db_dir}/${db_name} \\
            -max_target_seqs 1 \\
            -qcov_hsp_perc 70 \\
            -evalue 1e-10 \\
            -num_threads ${task.cpus} \\
            -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \\
            -out all_otus_blast_results.tsv
    fi

    # Capture BLAST version for provenance
    echo "BLAST+: \$(${blast_type} -version | head -n 1 | awk '{print \$2}')" > versions_blast.yml
    """
}