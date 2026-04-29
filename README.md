# Nextflow 18S Amplicon Analysis Pipeline

This Nextflow pipeline performs comprehensive analysis of 18S amplicon sequencing data, including preprocessing, demultiplexing, taxonomic classification, OTU clustering, and BLAST annotation.

## Features

- Adapter trimming and quality filtering
- Demultiplexing based on sample barcodes
- Taxonomic classification using Kraken2
- OTU clustering with custom amplicon_sorter
- BLAST annotation for taxonomic assignment
- Final result summarization
- Robust FASTQ input handling when identical FASTQ file names exist in different barcode directories

## Pipeline Flow and Output Structure

The pipeline runs in the following order:

1. `00_preprocess`: Adapter trimming and basic read cleanup per input `fastq_dir`
2. `01_demux`: Sample-level demultiplexing and primer-region extraction
3. `02_kraken_filter`: Kraken2 classification + KrakenTools target read extraction
4. `03_amplicon_sorter`: Clustering and consensus sequence generation per sample
5. `04_otu_merge`: Per-sample OTU counting and global OTU merge
6. `05_blast_annotation`: BLAST annotation of merged unique OTUs
7. `06_biom`: BIOM/TSV generation with taxonomy lineage integration
8. `07_summary_report`: Phylum-level summary tables
9. `06_final_results`: Consolidated deliverables for downstream use
10. `99_provenance`: Workflow run metadata and versions

Published output directory structure (`--outdir`):

```text
results/
├── 00_preprocess/
├── 01_demux/
├── 02_kraken_filter/
├── 03_amplicon_sorter/
├── 04_otu_merge/
├── 05_blast_annotation/
├── 06_biom/
│   ├── merged_results.biom
│   └── merged_otu_report.tsv
├── 06_final_results/
│   ├── otu_count_matrix.tsv
│   ├── sequences/
│   │   ├── integrated_unique_otus.fasta
│   │   └── all_samples_consensus.fasta
│   └── blast/
│       └── blast_annotation.tsv
├── 07_summary_report/
│   ├── summary_phylum_kraken2.tsv
│   └── summary_phylum_otu.tsv
└── 99_provenance/
```

## Prerequisites

- Tested with the following software:
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
- `as_script_path`: Path to the local `amplicon_sorter.py` script

`amplicon_sorter.py` source repository:
- https://github.com/avierstr/amplicon_sorter

#### Other Parameters
- `sample_sheet`: Path to sample sheet CSV (default in repo: `samplesheet/samplesheet_test_generated.csv`)
- `outdir`: Output directory (default: `results`)
- `target_taxid`: TaxID used in KrakenTools extraction (example: `7711`)
- `filter_action`: Read filtering mode for target taxonomy (example: `remove`)
- `max_reads`: Maximum reads for clustering (default: `10000`)
- `blast_db_name`: BLAST database name (default: `nt`)
- `blast_type`: BLAST algorithm (default: `blastn`)
- `max_cpus`, `max_memory`, `max_time`: Resource limits

### 4. Prepare Input Files

- **Sample Sheet**: CSV file with columns: sample, fastq_dir, min_len, max_len, fwd_index, fwd_primer, rev_index, rev_primer
- **FastQ Files**: Place in directories specified in sample_sheet. The pipeline preprocesses each `fastq_dir` once and deduplicates staged inputs by full file path, so identical FASTQ basenames across different directories are supported.
- **Databases**: Ensure Kraken2, BLAST, and taxonomy databases are accessible

#### Helper Script for Sample Sheet

Use the provided helper script to generate sample sheets:

**To generate target combinations from a template:**
```bash
python helperScript/generate_target_sample_sheet.py
```
This reads index/primer combinations from `template/18SV4-9_index.tsv` and prompts for F/R index selections to generate all combinations for the specified barcode range.

## Usage

Run the pipeline with Singularity profile:

```bash
nextflow run main.nf -profile singularity -params-file params.yaml
```

### Recommended Run Layout (per analysis)

For reproducibility and easier cleanup, keep each run in its own analysis directory:
- `analysisDir/work`: Nextflow work directory
- `analysisDir/results`: published outputs (`params.outdir` target)

Example:

```bash
analysisDir="./analysis/run01"
analysisDir="$(realpath -m "$analysisDir")"
mkdir -p "$analysisDir/work" "$analysisDir/results"

nextflow run /path/to/nextflow-18SAmplicon/main.nf \
    -profile singularity \
    -params-file /path/to/nextflow-18SAmplicon/params.yaml \
    -work-dir "$analysisDir/work" \
    --outdir "$analysisDir/results"
```

Notes:
- `--outdir` overrides `outdir` in `params.yaml`.
- `-work-dir` controls where Nextflow stores task working files and cache.
- This pattern works whether `analysisDir` is specified as a relative path (from launch directory) or an absolute path.

### Singularity Configuration

If a symlinked or current directory is used for data, ensure it is mounted. The pipeline automatically binds the network-mounted drive `/pigeon:/pigeon` for all Singularity containers. This ensures access to shared resources on your server.

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