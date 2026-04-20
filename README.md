# Nextflow 18S Amplicon Analysis Pipeline

This Nextflow pipeline performs comprehensive analysis of 18S amplicon sequencing data, including preprocessing, demultiplexing, taxonomic classification, OTU clustering, and BLAST annotation.

## Features

- Adapter trimming and quality filtering
- Demultiplexing based on sample barcodes
- Taxonomic classification using Kraken2
- OTU clustering with custom amplicon_sorter
- BLAST annotation for taxonomic assignment
- Final result summarization

## Prerequisites

- tested in following env.
- Nextflow (version 22.04.0)
- Singularity (version 3.8.6)
- Linux environment with sufficient computational resources (64 cores+, 32GB+ memory)

## Setup

### 1. Clone the Repository

```bash
git clone https://github.com/Tasu/nextflow-18SAmplicon.git
cd nextflow-18SAmplicon
```

### 2. Build Singularity Container

Build the required Singularity container for the amplicon sorter process:

```bash
chmod +x sif/build_amplicon_sorter.sh
./sif/build_amplicon_sorter.sh
```

This creates `sif/amplicon_sorter_v2.sif` with Python 3.9.7 and required libraries (edlib 1.3.9, biopython 1.79, matplotlib 3.6.2, numpy 1.21.2).

### 3. Configure Parameters

Edit `params.yaml` to set paths and parameters for your environment:

#### Required Database Paths
- `kraken2_db`: Path to Kraken2 database (e.g., `/path/to/kraken2_db`)
- `blast_db_dir`: Directory containing BLAST databases (e.g., `/path/to/blast_db`)
- `tax_db_dir`: Directory containing taxonomy databases (e.g., `/path/to/taxonomy`)

#### Container and Script Paths
- `as_container_path`: Path to the built Singularity container (e.g., `/absolute/path/to/sif/amplicon_sorter_v2.sif`)
- `as_script_path`: Path to the amplicon_sorter.py script (e.g., `/path/to/scripts/amplicon_sorter.py`)

#### Other Parameters
- `sample_sheet`: Path to sample sheet CSV (default: "sampleSheet.csv")
- `outdir`: Output directory (default: "results")
- `filter_taxids`: TaxIDs to filter out (default: "9606 10088")
- `max_reads`: Maximum reads for clustering (default: 10000)
- `blast_db_name`: BLAST database name (default: "nt")
- `max_cpus`, `max_memory`, `max_time`: Resource limits

### 4. Prepare Input Files

- **Sample Sheet**: CSV file with columns: sample, fastq_dir, min_len, max_len, fwd_index, fwd_primer, rev_index, rev_primer
- **FastQ Files**: Place in directories specified in sample_sheet
- **Databases**: Ensure Kraken2, BLAST, and taxonomy databases are accessible

#### Helper Script for Sample Sheet

Use the provided helper scripts to generate sample sheets:

**For scanning existing FASTQ files:**
```bash
python helperScript/create_sample_sheet.py <fastq_directory> [output_csv]
```

**For generating target combinations:**
```bash
python helperScript/generate_target_sample_sheet.py
```
This will prompt for iF and iR index arrays and generate all combinations for the specified barcode range.

Example:
```bash
python helperScript/create_sample_sheet.py demo/fastq samplesheet/samplesheet.csv
```

## Usage

Run the pipeline with Singularity profile:

```bash
nextflow run main.nf -profile singularity -params-file params.yaml
```

### Singularity Configuration

If simlinked directory is used for data, or current directory, please make sure to mount the destination directry. i.e. The pipeline automatically binds the network-mounted drive `/pigeon:/pigeon` for all Singularity containers. This ensures access to shared resources on your server.

If you need additional bind mounts, modify `nextflow.config`:

```groovy
singularity {
    runOptions = '--bind /pigeon:/pigeon --bind /additional:/path'
}
```

## Troubleshooting

- Ensure all database paths are accessible and correct
- Check Singularity container paths in `params.yaml`
- Verify network mounts are properly configured
- Monitor resource usage with `max_cpus`, `max_memory`, `max_time`

## License

To be determined - full permission for reuse and redistribution with proper citation.