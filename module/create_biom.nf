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
    awk -F'\\t' '{print \$NF}' ${blast_results} | \
        taxonkit lineage | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" > reformatted_lineage.txt
    echo "TaxonKit: \$(taxonkit version)" > versions_taxonkit.yml
    """
}

process BIOM_GENERATE {
    label 'process_medium'
    publishDir "${params.outdir}/01_Final_Outputs", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/biom-format:2.1.15'
    // Note: Biom-format image should include pandas. If not, quay.io/biocontainers/pandas:2.2.1 can be used with a combined script.

    input:
    path merged_counts
    path lineage
    path blast_results

    output:
    path "merged_results.biom",   emit: biom
    path "merged_otu_report.tsv", emit: tsv_report

    script:
    """
    # Python script for final matrix joining and BIOM export
    """
}