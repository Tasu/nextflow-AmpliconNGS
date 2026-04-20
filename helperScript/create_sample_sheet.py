#!/usr/bin/env python3
"""
Helper script to generate sample sheet CSV for Nextflow 18S Amplicon Pipeline.

Usage:
    python helperScript/create_sample_sheet.py <fastq_directory> [output_csv]

This script scans the specified directory for FASTQ files and generates a sample sheet
with default parameters. Edit the generated CSV manually for primers, indices, etc.

Example:
    python helperScript/create_sample_sheet.py demo/fastq samplesheet/samplesheet.csv
"""

import os
import sys
import glob
import csv

def create_sample_sheet(fastq_dir, output_csv):
    """
    Create a sample sheet CSV from FASTQ files in the directory.

    Assumes paired-end reads with naming convention: sample_R1.fastq.gz, sample_R2.fastq.gz
    or sample_1.fastq.gz, sample_2.fastq.gz
    """

    if not os.path.exists(fastq_dir):
        print(f"Error: Directory {fastq_dir} does not exist.")
        sys.exit(1)

    # Find all FASTQ files
    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) + \
                  glob.glob(os.path.join(fastq_dir, "*.fq.gz")) + \
                  glob.glob(os.path.join(fastq_dir, "*.fastq")) + \
                  glob.glob(os.path.join(fastq_dir, "*.fq"))

    if not fastq_files:
        print(f"No FASTQ files found in {fastq_dir}")
        sys.exit(1)

    # Group by sample name (remove _R1, _R2, _1, _2 suffixes)
    samples = {}
    for f in fastq_files:
        basename = os.path.basename(f)
        # Remove extensions
        name = basename.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
        # Remove read suffixes
        for suffix in ['_R1', '_R2', '_1', '_2']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        if name not in samples:
            samples[name] = []
        samples[name].append(f)

    # Create CSV
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['sample', 'fastq_dir', 'min_len', 'max_len', 'fwd_index', 'fwd_primer', 'rev_index', 'rev_primer']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for sample, files in samples.items():
            # Sort files to ensure R1/R2 order
            files.sort()
            writer.writerow({
                'sample': sample,
                'fastq_dir': fastq_dir,
                'min_len': 300,  # Default values - adjust as needed
                'max_len': 500,
                'fwd_index': '',  # Fill in manually
                'fwd_primer': 'CCTACGGGNGGCWGCAG',  # Example 18S primers - adjust
                'rev_index': '',
                'rev_primer': 'GACTACHVGGGTATCTAATCC'
            })

    print(f"Sample sheet created: {output_csv}")
    print(f"Found {len(samples)} samples with FASTQ files.")
    print("Please edit the CSV to add correct primer sequences and indices.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    fastq_dir = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "samplesheet/samplesheet.csv"

    create_sample_sheet(fastq_dir, output_csv)