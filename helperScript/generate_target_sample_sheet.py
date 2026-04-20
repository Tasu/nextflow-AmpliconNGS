#!/usr/bin/env python3
"""
Helper script to generate sample sheet CSV for Nextflow 18S Amplicon Pipeline.
This version generates combinations based on barcode ranges and index arrays.

Usage:
    python helperScript/generate_target_sample_sheet.py

This will prompt for index arrays and generate a sample sheet with all combinations.
"""

import sys
import csv

def get_index_arrays():
    """Get iF and iR index arrays from user input."""
    print("Enter iF indices (comma-separated, e.g., F01,F02,F03):")
    if_input = input().strip()
    if_list = [x.strip() for x in if_input.split(',') if x.strip()]

    print("Enter iR indices (comma-separated, e.g., R01,R02,R03):")
    ir_input = input().strip()
    ir_list = [x.strip() for x in ir_input.split(',') if x.strip()]

    return if_list, ir_list

def generate_data(start_num, end_num, if_list, ir_list, base_row):
    """Generate data rows for all combinations."""
    rows = []

    for i in range(start_num, end_num + 1):
        formatted_num = f"{i:02d}"
        barcode_name = f"barcode{formatted_num}"

        # iF × iR combinations
        for f_val in if_list:
            for r_val in ir_list:
                # sample name as barcode_Fxx_Rxx
                full_sample_name = f"{barcode_name}_{f_val}_{r_val}"

                row = [
                    base_row["date"],
                    base_row["desc"],
                    full_sample_name,       # sample
                    base_row["flow_cell"],
                    base_row["kit"],
                    f_val,              # iF
                    r_val,              # iR
                    barcode_name,       # barcode
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

    # Get index arrays from user
    if_list, ir_list = get_index_arrays()

    # Generate data
    all_data = generate_data(start_num, end_num, if_list, ir_list, base_row)

    # Create output file
    try:
        with open(file_name, mode='w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(header)
            writer.writerows(all_data)

        print(f"Success: '{file_name}' created.")
        print(f"Details: {len(all_data)} rows generated.")
        barcode_count = end_num - start_num + 1
        combinations = len(if_list) * len(ir_list)
        print(f"(Barcodes {start_num}-{end_num}: {barcode_count} types × {combinations} combinations = {barcode_count * combinations} total)")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()