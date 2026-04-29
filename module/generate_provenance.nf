/*
 * GENERATE_PROVENANCE Module
 * * Description:
 * Aggregates version information from all processes and creates a 
 * Markdown report and a plain text file for audit trails.
 * * Requirements:
 * - Container with Python 3.11+
 */

process GENERATE_PROVENANCE {
    label 'process_low'
    publishDir "${params.outdir}/99_provenance", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'

    input:
    // Stage each collected version file into a separate numbered subdirectory
    // to avoid filename collisions (e.g., multiple versions_seqkit.yml).
    path v_files, stageAs: 'versions??/*'
    val workflow_info // workflow object from main.nf

    output:
    path "software_versions.md",  emit: md_report
    path "software_versions.txt", emit: txt_report
    path "provenance.yml",        emit: yaml_summary

    script:
    """
    python3 << 'EOF'
    import datetime

    # 1. Collect versions from all input files
    version_data = {}
    files = "${v_files}".split()
    for f in files:
        with open(f, 'r') as fh:
            for line in fh:
                if ':' in line:
                    key, val = line.strip().split(':', 1)
                    version_data[key.strip()] = val.strip()

    # 2. Generate Markdown Report
    with open("software_versions.md", "w") as md:
        md.write("# Software Provenance Report\\n\\n")
        md.write(f"**Analysis Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        md.write(f"**Nextflow Version:** ${workflow_info.nextflow.version}\\n")
        md.write(f"**Command Line:** `${workflow_info.commandLine}`\\n\\n")
        md.write("| Software/Tool | Version |\\n")
        md.write("| :--- | :--- |\\n")
        for tool, ver in sorted(version_data.items()):
            md.write(f"| {tool} | {ver} |\\n")

    # 3. Generate plain text and yaml summary
    with open("software_versions.txt", "w") as txt:
        for tool, ver in sorted(version_data.items()):
            txt.write(f"{tool}: {ver}\\n")

    with open("provenance.yml", "w") as yml:
        yml.write("software_versions:\\n")
        for tool, ver in sorted(version_data.items()):
            yml.write(f"  {tool}: {ver}\\n")
    EOF
    """
}