/*
 * OTU_COUNT Module
 * * Description:
 * Processes the OTU count matrix and incorporates total read counts per sample
 * to provide a normalized or summarized view of the community composition.
 * * Requirements:
 * - Container with Biopython and Python 3.11+
 */

process OTU_COUNT_TABLE {
    label 'process_medium'
    publishDir "${params.outdir}/04_otu_merge", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container "${params.container_images.biopython}"

    input:
    tuple val(sample_id), path(consensus_fasta), path(read_count_file)

    output:
    path "${sample_id}_otu_table_final.tsv",     emit: final_table
    path "${sample_id}_otu_table_summary.txt",   emit: summary
    path "${sample_id}_versions_otu_count.yml",  emit: versions

    script:
    """
    # Count OTUs in the per-sample consensus FASTA and record total reads used
    python3 << 'EOF'
from Bio import SeqIO

sample_id = "${sample_id}"

# Count OTUs (sequences) in consensus FASTA
otu_count = sum(1 for _ in SeqIO.parse("${consensus_fasta}", "fasta"))

# Read total reads used
with open("${read_count_file}", 'r') as fh:
    total_reads = int(fh.read().strip())

# Write per-sample OTU table
with open(f"{sample_id}_otu_table_final.tsv", 'w') as f_out:
    f_out.write("Sample\\tOTU_Count\\tTotal_Reads\\n")
    f_out.write(f"{sample_id}\\t{otu_count}\\t{total_reads}\\n")

# Write summary
with open(f"{sample_id}_otu_table_summary.txt", 'w') as summary:
    summary.write("Sample\\tOTU_Count\\tTotal_Reads_Used\\n")
    summary.write(f"{sample_id}\\t{otu_count}\\t{total_reads}\\n")

EOF

    # Capture versions
    echo "Python (Biopython container): \$(python3 --version)" > ${sample_id}_versions_otu_count.yml
    echo "Biopython: 1.79" >> ${sample_id}_versions_otu_count.yml
    """
}