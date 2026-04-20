/*
 * OTU_MERGE Module
 * * Description:
 * Integrates clustered consensus sequences (OTUs) from all samples,
 * performs 100% identity clustering to identify unique global OTUs,
 * and generates a count matrix (OTU table).
 * * Requirements:
 * - Container with Biopython and Python 3.11+
 */

process OTU_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}/04_otu_merge", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79--py311h1425ee9_1'

    input:
    tuple val(sample_ids), path(consensus_fastas) // List of all consensus files
    tuple val(sample_ids_count), path(count_files)   // List of read_count_actual.txt files

    output:
    path "integrated_unique_otus.fasta", emit: otu_fasta
    path "otu_count_matrix.tsv",         emit: count_matrix
    path "versions_otu_merge.yml",       emit: versions

    script:
    """
    # 1. Integrate all sequences into a single global FASTA
    # Use python script to parse fasta and track sample origins
    python3 << 'EOF'
import os
from Bio import SeqIO
from collections import defaultdict

otu_counts = defaultdict(lambda: defaultdict(int))
unique_seqs = {}

# Iterate through each sample's consensus file
fastas = "${consensus_fastas}".split()
for f in fastas:
    sample_id = os.path.basename(f).replace("_clustered_consensus.fasta", "")
    for record in SeqIO.parse(f, "fasta"):
        seq = str(record.seq).upper()
        # Use sequence itself as key for 100% identity merging
        if seq not in unique_seqs:
            otu_id = f"OTU_{len(unique_seqs) + 1}"
            unique_seqs[seq] = otu_id
        
        target_otu = unique_seqs[seq]
        # In this logic, each cluster from Amplicon_Sorter is treated as 1 unit or weighted by its size if available
        # Here we increment by 1 for the existence of the cluster in the sample
        otu_counts[target_otu][sample_id] += 1

# Write integrated unique OTUs to FASTA
with open("integrated_unique_otus.fasta", "w") as fa:
    for seq, otu_id in unique_seqs.items():
        fa.write(f">{otu_id}\\n{seq}\\n")

# Write OTU count matrix (TSV)
samples = sorted("${sample_ids}".split())
with open("otu_count_matrix.tsv", "w") as tsv:
    tsv.write("OTU_ID\\t" + "\\t".join(samples) + "\\n")
    for otu_id in sorted(unique_seqs.values(), key=lambda x: int(x.split('_')[1])):
        row = [otu_id]
        for s in samples:
            row.append(str(otu_counts[otu_id][s]))
        tsv.write("\\t".join(row) + "\\n")
EOF

    # 2. Capture versions for provenance
    echo "Python (Biopython container): \$(python3 --version)" > versions_otu_merge.yml
    echo "Biopython: 1.79" >> versions_otu_merge.yml
    """
}