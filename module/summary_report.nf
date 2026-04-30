/*
 * SUMMARY_REPORT Module
 * * Description:
 * Generates summary tables at the Phylum level by aggregating 
 * Kraken2 taxonomy reports and merging OTU count matrices with BLAST annotations.
 * * Requirements:
 * - Container with Python 3.11+ (Biopython container)
 */

process SUMMARY_REPORT {
    label 'process_medium'
    publishDir "${params.outdir}/07_summary_report", mode: 'copy'

    // Using the consistent Python 3.11 based container
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'

    input:
    path kraken_reports    // Collected from KRAKEN2_CLASSIFY.out.report
    path otu_matrix        // From OTU_MERGE.out.count_matrix
    path blast_results     // From BLAST_ANNOTATE.out.blast_results

    output:
    path "summary_phylum_kraken2.tsv", emit: kraken_summary
    path "summary_phylum_otu.tsv",     emit: otu_summary
    path "versions_summary.yml",       emit: versions

    script:
    """
    python3 << 'EOF'
import os
from collections import defaultdict

# --- 1. Aggregate Kraken2 Reports at Phylum Level ---
kraken_phylum_data = defaultdict(lambda: defaultdict(int))
samples = []

report_files = "${kraken_reports}".split()
for f in report_files:
    # Infer sample_id from filename (e.g., sample_id.kraken2.report.txt)
    sample_id = os.path.basename(f).split('.')[0]
    samples.append(sample_id)

    with open(f, 'r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue
            # Kraken2 format: percentage, count, self_count, rank, taxid, name
            rank = cols[3]
            count = int(cols[1])
            name = cols[5].strip()

            if rank == 'P':  # Phylum level
                kraken_phylum_data[name][sample_id] += count

samples = sorted(list(set(samples)))
with open("summary_phylum_kraken2.tsv", "w") as out:
    out.write("Phylum\t" + "\t".join(samples) + "\n")
    for phylum in sorted(kraken_phylum_data.keys()):
        row = [phylum] + [str(kraken_phylum_data[phylum][s]) for s in samples]
        out.write("\t".join(row) + "\n")


# --- 2. Aggregate OTU Table at Phylum Level using BLAST Annotations ---
# Load BLAST results: OTU_ID -> taxonomy string
otu_to_phylum = {}
with open("${blast_results}", 'r') as bf:
    for line in bf:
        cols = line.strip().split('\t')
        if len(cols) < 2:
            continue
        otu_id = cols[0]
        full_tax = cols[1]
        otu_to_phylum[otu_id] = full_tax

# Load OTU Matrix and collapse
# OTUs with no BLAST hit are assigned "NA" (not excluded)
otu_phylum_summary = defaultdict(lambda: defaultdict(int))
otu_samples = []

with open("${otu_matrix}", 'r') as mf:
    header = mf.readline().strip().split('\t')
    otu_samples = header[1:]
    for line in mf:
        cols = line.strip().split('\t')
        otu_id = cols[0]
        counts = [int(x) for x in cols[1:]]

        phylum = otu_to_phylum.get(otu_id, "NA")
        for i, s in enumerate(otu_samples):
            otu_phylum_summary[phylum][s] += counts[i]

with open("summary_phylum_otu.tsv", "w") as out:
    out.write("Phylum\t" + "\t".join(otu_samples) + "\n")
    for phylum in sorted(otu_phylum_summary.keys()):
        row = [phylum] + [str(otu_phylum_summary[phylum][s]) for s in otu_samples]
        out.write("\t".join(row) + "\n")

EOF

    # Capture versions
    echo "Python (Summary): \$(python3 --version)" > versions_summary.yml
    echo "Biopython: 1.79" >> versions_summary.yml
    """
}