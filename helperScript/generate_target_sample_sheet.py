#!/usr/bin/env python3
"""
Helper script to generate sample sheet CSV for Nextflow 18S Amplicon Pipeline.
This version reads index/primer combinations from template file and generates
sample sheets based on barcode ranges and selected index combinations.

Requirement:
    - Python 3
    - pandas (for index sequence conversion)

Usage:
    python helperScript/generate_target_sample_sheet.py [fastq_directory] [output_csv]

Arguments:
    fastq_directory: Path to FASTQ directory (default: demo/fastq)
    output_csv: Path to output CSV file (default: samplesheet/sampleSheet.target.csv)

This will read the template file, auto-detect barcode range, prompt for index
combinations, and generate a sample sheet with all combinations.

Examples:
    python helperScript/generate_target_sample_sheet.py
    python helperScript/generate_target_sample_sheet.py demo/fastq
    python helperScript/generate_target_sample_sheet.py /path/to/fastq samplesheet/output.csv
"""

import sys
import csv
import os
import re
import glob
from io import StringIO
from pathlib import Path
import pandas as pd


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

def detect_barcode_range(fastq_dir):
    """Detect barcode range from directory structure."""
    barcode_pattern = re.compile(r'barcode(\d+)')
    barcode_numbers = []

    if not os.path.exists(fastq_dir):
        return None, None

    # Find all barcode** directories
    for item in os.listdir(fastq_dir):
        item_path = os.path.join(fastq_dir, item)
        if os.path.isdir(item_path):
            match = barcode_pattern.match(item)
            if match:
                barcode_numbers.append(int(match.group(1)))

    if barcode_numbers:
        return min(barcode_numbers), max(barcode_numbers)
    return None, None

def get_barcode_range_from_user(fastq_dir="demo/fastq"):
    """Get barcode range from user with auto-detection as default."""
    start_detected, end_detected = detect_barcode_range(fastq_dir)

    if start_detected is not None and end_detected is not None:
        print(f"\nDetected barcodes in {fastq_dir}: barcode{start_detected:02d} - barcode{end_detected:02d}")
        default_str = f"{start_detected}-{end_detected}"
        print(f"Enter barcode range [default: {default_str}] (e.g., 28-31):")
        user_input = input().strip()
        
        if not user_input:
            return start_detected, end_detected
        
        try:
            parts = user_input.split('-')
            return int(parts[0].strip()), int(parts[1].strip())
        except (ValueError, IndexError):
            print("Invalid input. Using detected range.")
            return start_detected, end_detected
    else:
        print(f"No barcodes found in {fastq_dir}")
        print("Enter barcode range (e.g., 28-31):")
        user_input = input().strip()
        try:
            parts = user_input.split('-')
            return int(parts[0].strip()), int(parts[1].strip())
        except (ValueError, IndexError):
            print("Invalid input.")
            sys.exit(1)
def convert_indices_to_sequence(indices):
    """Convert index names (e.g., F01) to actual sequences."""
    # This is a placeholder function. You would replace this with actual logic
    # to convert index names to sequences based on your specific requirements.
    

    script_dir = Path(__file__).parent
    template_file = script_dir / "18SV4-9_index.tsv"
    df = pd.read_csv(template_file, sep='\t')

    index_sequences = {}

    # 2. 各行を処理して辞書に登録していく
    for _, row in df.iterrows():
        # SampleID (例: 01_02) を "_" で分割
        sample_parts = str(row['SampleID']).split('_')
        if len(sample_parts) != 2:
            continue
            
        f_num, r_num = sample_parts
        
        # F01, R02 のようなキー名を作成
        f_key = f"F{f_num}"
        r_key = f"R{r_num}"
        
        # まだ辞書になければ、配列（Index）を登録
        if f_key not in index_sequences:
            index_sequences[f_key] = row['FwIndex']
        
        if r_key not in index_sequences:
            index_sequences[r_key] = row['RvIndex']


    # index_sequences = {
    #     "F01": "ACGT",
    #     "F02": "TGCA",
    #     "R01": "GCTA",
    #     "R02": "CAGT",
    #     # Add more mappings as needed
    # }
    return [index_sequences.get(idx, "UNKNOWN") for idx in indices]



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
            f_index_seq, r_index_seq = convert_indices_to_sequence([combo['f_idx'], combo['r_idx']])
            row = [
                full_sample_name,               # sample
                base_row["path"] + "/" + barcode_name,  # fastq_dir
                base_row["min"],
                base_row["max"],
                f_index_seq,                 # fwd_index
                combo['fwd_primer'],
                r_index_seq,                 # rev_index
                combo['rev_primer']
            ]
            rows.append(row)
    return rows

def main(fastq_dir=None, output_csv=None):
    # Configuration
    script_dir = Path(__file__).parent
    template_file = script_dir / "18SV4-9_index.tsv"
    if fastq_dir is None:
        fastq_dir = "demo/fastq"
    if output_csv is None:
        output_csv = "samplesheet/sampleSheet.target.csv"
    delimiter = ","  # CSV delimiter

    # Get barcode range from user (with auto-detection)
    start_num, end_num = get_barcode_range_from_user(fastq_dir)

    # Header labels for pipeline-compatible sample sheet
    header = [
        "sample", "fastq_dir", "min_len", "max_len",
        "fwd_index", "fwd_primer", "rev_index", "rev_primer"
    ]

    # Fixed data - adjust as needed
    base_row = {
        "date": "CHANGE_DATE_PROJECT",
        "desc": "CHANGE_DESCRIPTION",
        "flow_cell": "FLO-FLG114",
        "kit": "SQK-NBD114.96",
        "caller": "8",
        "path": fastq_dir,
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
        # パスが "-" なら sys.stdout を使い、それ以外なら open でファイルを開く
        with (sys.stdout if output_csv == "-" else open(output_csv, mode='w', encoding='utf-8', newline='')) as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(header)
            writer.writerows(all_data)

        # 標準出力（画面）にデータを出した場合は、以下のメッセージを表示させないか、stderrに出す
        if output_csv != "-":
            print(f"\nSuccess: '{output_csv}' created.")
            print(f"Details: {len(all_data)} rows generated.")
            # ...その他のメッセージ...
        else:
            # 画面にデータを出した場合でも、詳細をログとして出したいなら stderr を使う
            print(f"Details: {len(all_data)} rows generated.", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)

if __name__ == "__main__":
    fastq_dir = sys.argv[1] if len(sys.argv) > 1 else "demo/fastq"
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "samplesheet/sampleSheet.target.csv"
    main(fastq_dir, output_csv)