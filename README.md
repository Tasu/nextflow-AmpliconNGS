# Nextflow NGS Amplicon Analysis Pipeline

This Nextflow pipeline performs comprehensive analysis of amplicon sequencing data, including preprocessing, demultiplexing, taxonomic classification, OTU clustering, and BLAST annotation. Still under construction.

## Features
- helperScript to prepare sampleSheet for nextflow pipeline execute.
- Adapter trimming and quality filtering
- Demultiplexing based on nanopore barcode and index in primer
- Taxonomic classification using Kraken2 for filtering host sequence.
- OTU clustering with custom amplicon_sorter
- BLAST annotation for taxonomic assignment
- Final result summarization

## Pipeline Flow and Output Structure

The pipeline runs in the following order:

1. `00_preprocess`: Adapter trimming and basic read cleanup per input `fastq_dir`
2. `01_demux`: Sample-level demultiplexing and primer-region extraction
3. `02_kraken_filter`: Kraken2 classification + KrakenTools target read extraction/filtration
4. `03_amplicon_sorter`: Clustering and consensus sequence generation per sample
5. `04_otu_merge`: Per-sample OTU counting and global OTU merge
6. `05_blast_annotation`: BLAST annotation of merged unique OTUs
7. `07_biom`: BIOM/TSV generation with taxonomy lineage integration
8. `08_summary_report`: Phylum-level summary tables
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
├── 06_final_results/
│   ├── otu_count_matrix.tsv
│   ├── sequences/
│   │   ├── integrated_unique_otus.fasta
│   │   └── all_samples_consensus.fasta
│   └── blast/
│       └── blast_annotation.tsv
├── 07_biom/
│   ├── merged_results.biom
│   └── merged_otu_report.tsv
├── 08_summary_report/
│   ├── summary_phylum_kraken2.tsv
│   └── summary_phylum_otu.tsv
└── 99_provenance/
```

## Prerequisites

- Tested with the following software:
- Nextflow (version 22.10.6)
- Singularity (version 3.8.6)
- Skopeo (version 1.22.2, required for `docker://` preflight checks)
- Linux environment with sufficient computational resources (64 cores+, 32GB+ memory)

Example conda environment setup (`nf-env`, based on `version.ipynb`):

```bash
conda create -n nf-env -c bioconda -c conda-forge \
    nextflow=22.10.6 \
    singularity=3.8.6
# if you want to do full preflight check,
conda install -n nf-env -c conda-forge skopeo=1.22.2
```

Environment check commands:

```bash
conda run -n nf-env nextflow -version
conda run -n nf-env singularity --version
conda run -n nf-env skopeo --version
```

## Setup

### 1. Clone the Repository

```bash
git clone https://github.com/Tasu/nextflow-18SAmplicon.git
cd nextflow-18SAmplicon
```

### 2. Build required additional Singularity Container

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

#### Container Image URLs

- All module container URIs are centrally managed in `nextflow.config` via `params.container_images`.
- To update image tags, edit only this map (module files do not need updates).
- For simple operation, set only selected keys in `params.yaml` under `container_images` (unspecified keys keep defaults).
- Optional pre-run checks:
- `container_preflight_check`: enable URL connectivity checks before workflow execution (default: `true`)
- `container_preflight_strict`: fail early on unreachable URLs (default: `false`, warn only)
- `container_preflight_suggest_tags`: suggest candidate tags from Galaxy Depot index when a URL is unreachable (default: `true`)
- `container_preflight_check_docker`: verify `docker://` references with `skopeo inspect` (default: `true`)
- `container_registry_index_url`: index used for tag suggestions (default: `https://depot.galaxyproject.org/singularity/`)
- `preflight_only`: run only preflight checks and exit before analysis starts (default: `false`)

Example (`params.yaml`):

```yaml
container_images:
    blast: "https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4"
    biopython: "https://depot.galaxyproject.org/singularity/biopython:1.79"
```

`amplicon_sorter.py` source repository:

- <https://github.com/avierstr/amplicon_sorter>

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
# custom fastq dir and output path (optional)
python helperScript/generate_target_sample_sheet.py --fastq-dir demo/fastq --output samplesheet/sampleSheet.target.csv
# custom template path (optional)
python helperScript/generate_target_sample_sheet.py --template /path/to/18SV4-9_index.tsv
```

This reads index/primer combinations from `helperScript/18SV4-9_index.tsv` by default and prompts for F/R index selections to generate all combinations for the specified barcode range.

## Usage

Run the pipeline with Singularity profile:

```bash
nextflow run main.nf -profile singularity -params-file params.yaml
```

If you use Nextflow from conda environment `nf-env`, run:

```bash
conda run -n nf-env nextflow run main.nf -profile singularity -params-file params.yaml
```

### Preflight-Only Mode

Use this mode to run only container URL preflight checks and exit before analysis tasks start.

```bash
nextflow run main.nf -profile singularity -params-file params.yaml --preflight_only true
```

With `nf-env`:

```bash
conda run -n nf-env nextflow run main.nf -profile singularity -params-file params.yaml --preflight_only true
```

CI/release validation example (fail on unreachable container URLs):

```bash
conda run -n nf-env nextflow run main.nf \
    -profile singularity \
    -params-file params.yaml \
    --preflight_only true \
    --container_preflight_strict true
```

Expected behavior:

- Preflight logs are printed for each configured container URL.
- Workflow exits before any analysis process is scheduled.
- Exit code is non-zero only when strict mode is enabled and checks fail.

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

If a symlinked or current directory is used for data, ensure it is mounted. The pipeline automatically binds the network-mounted drive `/dataDrive:/dataDrive` for all Singularity containers. This ensures access to shared resources on your server.

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

This repository is licensed under the MIT License. See the `LICENSE` file.

Third-party tools and databases used by this workflow are licensed separately by their original authors. Their terms still apply when you run this pipeline, including when tools are executed via Singularity/Apptainer containers.

See `THIRD_PARTY_LICENSES.md` for usage and redistribution notes.
