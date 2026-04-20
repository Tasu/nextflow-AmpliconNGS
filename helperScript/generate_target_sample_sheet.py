#!/usr/bin/env python3
"""
Helper script to generate sample sheet CSV for Nextflow 18S Amplicon Pipeline.
This version reads index/primer combinations from template file and generates
sample sheets based on barcode ranges and selected index combinations.

Usage:
    python helperScript/generate_target_sample_sheet.py

This will read the template file, prompt for index combinations, and generate
a sample sheet with all combinations for the specified barcode range.
"""

import sys
import csv

def load_template(template_file):
    """Load index/primer combinations from template TSV file."""
    combinations = {}
    try:
        with open(template_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample_id = row['SampleID']
                # Parse F and R indices from SampleID (e.g., 01_05 -> F01, R05)
                f_num, r_num = sample_id.split('_')
                f_idx = f"F{f_num}"
                r_idx = f"R{r_num}"

                combinations[sample_id] = {
                    'f_idx': f_idx,
                    'r_idx': r_idx,
                    'fwd_primer': row['FwPrimer'],
                    'rev_primer': row['RvPrimer']
                }
    except FileNotFoundError:
        print(f"Error: Template file {template_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading template file: {e}")
        sys.exit(1)

    return combinations

def get_available_indices(combinations):
    """Get unique F and R indices from combinations."""
    f_indices = sorted(set([combo['f_idx'] for combo in combinations.values()]))
    r_indices = sorted(set([combo['r_idx'] for combo in combinations.values()]))
    return f_indices, r_indices

def select_combinations(combinations, selected_f, selected_r):
    """Get combinations for selected F and R indices."""
    selected_combos = {}
    for sample_id, combo in combinations.items():
        if combo['f_idx'] in selected_f and combo['r_idx'] in selected_r:
            selected_combos[sample_id] = combo
    return selected_combos

def generate_data(start_num, end_num, selected_combos, base_row):
    """Generate data rows for selected combinations."""
    rows = []

    for i in range(start_num, end_num + 1):
        formatted_num = f"{i:02d}"
        barcode_name = f"barcode{formatted_num}"

        # For each selected combination
        for sample_id, combo in selected_combos.items():
            # sample name as barcode_Fxx_Rxx
            full_sample_name = f"{barcode_name}_{combo['f_idx']}_{combo['r_idx']}"

            row = [
                base_row["date"],
                base_row["desc"],
                full_sample_name,       # sample
                base_row["flow_cell"],
                base_row["kit"],
                combo['f_idx'],         # iF
                combo['r_idx'],         # iR
                barcode_name,           # barcode
                base_row["caller"],
                base_row["path"],
                base_row["min"],
                base_row["max"],
                base_row["maxreads"],
                base_row["dummy"]
            ]
            rows.append(row)
    return rows

def main():
    # Configuration
    template_file = "template/18SV4-9_index.tsv"
    file_name = "samplesheet/sampleSheet.target.csv"
    start_num = 28  # Start number
    end_num = 31    # End number
    delimiter = ","  # CSV delimiter

    # Header labels
    header = [
        "#Date", "Project description", "sample", "flow cell", "lib kit",
        "iF", "iR", "barcode", "fastq_base_caller", "fastq_passDir",
        "min", "max", "maxreads", "DUMMY"
    ]

    # Fixed data - adjust as needed
    base_row = {
        "date": "20260415_TS_ShikabetaMammal",
        "desc": "shikabetaIndexed",
        "flow_cell": "FLO-FLG114",
        "kit": "SQK-NBD114.96",
        "caller": "8",
        "path": "/home/tsugi/work/rawread/Nanopore/fastq/20260415_TS_ShikabetaMammal/basecalling/pass",
        "min": "1000",
        "max": "1800",
        "maxreads": "10000",
        "dummy": "A"
    }

    # Load template
    combinations = load_template(template_file)
    f_indices, r_indices = get_available_indices(combinations)

    print(f"Available F indices: {', '.join(f_indices)}")
    print(f"Available R indices: {', '.join(r_indices)}")

    # Get user selection
    print("\nEnter F indices to use (comma-separated, e.g., F01,F02):")
    f_input = input().strip()
    selected_f = [x.strip() for x in f_input.split(',') if x.strip()]

    print("Enter R indices to use (comma-separated, e.g., R01,R05):")
    r_input = input().strip()
    selected_r = [x.strip() for x in r_input.split(',') if x.strip()]

    # Filter combinations
    selected_combos = select_combinations(combinations, selected_f, selected_r)

    if not selected_combos:
        print("No combinations found for selected indices.")
        sys.exit(1)

    print(f"\nSelected {len(selected_combos)} combinations:")
    for sample_id, combo in selected_combos.items():
        print(f"  {sample_id}: {combo['f_idx']}, {combo['r_idx']}")

    # Generate data
    all_data = generate_data(start_num, end_num, selected_combos, base_row)

    # Create output file
    try:
        with open(file_name, mode='w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(header)
            writer.writerows(all_data)

        print(f"\nSuccess: '{file_name}' created.")
        print(f"Details: {len(all_data)} rows generated.")
        barcode_count = end_num - start_num + 1
        combo_count = len(selected_combos)
        print(f"(Barcodes {start_num}-{end_num}: {barcode_count} × {combo_count} combinations = {barcode_count * combo_count} total)")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()