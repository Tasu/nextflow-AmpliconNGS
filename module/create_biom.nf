/*
 * CREATE_BIOM Module: Separated into TaxonKit and Python(Pandas/Biom)
 */

process TAXONKIT_LINEAGE {
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/taxonkit:0.18.0--h9ee0642_0'

    input:
    path blast_results
    path tax_db_dir

    output:
    path "reformatted_lineage.txt", emit: lineage
    path "versions_taxonkit.yml",  emit: versions

    script:
    """
    export TAXONKIT_DB="${tax_db_dir}"

    # If BLAST results are empty (no hits at all), emit an empty lineage file and exit cleanly.
    if [[ ! -s "${blast_results}" ]]; then
        echo "WARNING: BLAST results are empty — skipping TaxonKit lineage lookup." >&2
        touch reformatted_lineage.txt
    else
        # Extract staxids (last column), filter out missing/invalid values (0, N/A, empty),
        # then run taxonkit. This prevents taxonkit errors when a BLAST hit has no taxonomy.
        awk -F'\\t' '\$NF != "" && \$NF != "N/A" && \$NF != "0" {print \$NF}' ${blast_results} \
            | sort -u \
            | taxonkit lineage \
            | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" \
            > reformatted_lineage.txt || true

        # Ensure file exists even if taxonkit produced no output
        touch reformatted_lineage.txt
    fi

    echo "TaxonKit: \$(taxonkit version)" > versions_taxonkit.yml
    """
}

process BIOM_GENERATE {
    label 'process_medium'
    publishDir "${params.outdir}/06_biom", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/biom-format:2.1.15'
    // Note: Biom-format image should include pandas. If not, quay.io/biocontainers/pandas:2.2.1 can be used with a combined script.

    input:
    path merged_counts
    path lineage
    path blast_results

    output:
    path "merged_results.biom",   emit: biom
    path "merged_otu_report.tsv", emit: tsv_report
    path "versions_biom.yml",     emit: versions

    script:
    """
    python3 << 'EOF'
import csv
import sys

NA_LINEAGE = "NA;NA;NA;NA;NA;NA;NA"

# --- 1. Load TaxonKit lineage: taxid -> lineage string ---
# reformatted_lineage.txt has two columns: taxid<TAB>lineage
taxid_to_lineage = {}
with open("${lineage}", 'r') as fh:
    for line in fh:
        parts = line.rstrip("\\n").split("\t")
        if len(parts) >= 2 and parts[0].strip():
            taxid_to_lineage[parts[0].strip()] = parts[1].strip()

# --- 2. Load BLAST results: otu_id -> (staxids, stitle) ---
# OTUs absent from blast_results will receive NA taxonomy
otu_blast = {}  # otu_id -> lineage string
with open("${blast_results}", 'r') as fh:
    for line in fh:
        cols = line.rstrip("\\n").split("\t")
        if len(cols) < 14:
            continue
        otu_id = cols[0]
        staxids = cols[13].strip()
        # Use first taxid if multiple (semicolon-separated)
        taxid = staxids.split(";")[0] if staxids not in ("", "N/A", "0") else ""
        lineage = taxid_to_lineage.get(taxid, NA_LINEAGE) if taxid else NA_LINEAGE
        otu_blast[otu_id] = lineage

# --- 3. Load OTU count matrix ---
with open("${merged_counts}", 'r') as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    samples = header[1:]
    otu_rows = [row for row in reader]

# --- 4. Build merged TSV (left join — OTUs with no BLAST hit get NA) ---
tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
out_header = ["OTU_ID"] + tax_levels + samples

with open("merged_otu_report.tsv", 'w') as out:
    out.write("\t".join(out_header) + "\\n")
    for row in otu_rows:
        otu_id = row[0]
        counts = row[1:]
        lineage_str = otu_blast.get(otu_id, NA_LINEAGE)
        tax_fields = (lineage_str.split(";") + ["NA"] * 7)[:7]
        out.write("\t".join([otu_id] + tax_fields + counts) + "\\n")

# --- 5. Convert to BIOM (JSON-based format v1) ---
import json, time

rows = [{"id": row[0], "metadata": {
    "taxonomy": (otu_blast.get(row[0], NA_LINEAGE).split(";") + ["NA"] * 7)[:7]
}} for row in otu_rows]

cols = [{"id": s, "metadata": None} for s in samples]

data = []
for r_idx, row in enumerate(otu_rows):
    for c_idx, val in enumerate(row[1:]):
        v = int(val)
        if v != 0:
            data.append([r_idx, c_idx, v])

biom = {
    "id": None,
    "format": "Biological Observation Matrix 1.0.0",
    "format_url": "http://biom-format.org",
    "type": "OTU table",
    "generated_by": "nextflow-18SAmplicon",
    "date": "",
    "rows": rows,
    "columns": cols,
    "matrix_type": "sparse",
    "matrix_element_type": "int",
    "shape": [len(rows), len(cols)],
    "data": data
}

with open("merged_results.biom", 'w') as bf:
    json.dump(biom, bf)

EOF

    echo "BIOM-format: 2.1.15, Pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')" > versions_biom.yml
    """
}