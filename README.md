# python_file_splitter
Splits files into whatever number of sub files keeping column names for faster procesing

split_tsv.py
```py


#!/usr/bin/env python3
"""
split_tsv.py

Split a large TSV file into N parts for job arrays, preserving the header in each part.
Outputs go to a user-specified directory (default: 'splitted_input/'). Filenames are
`<prefix><NN>.tsv` where numbering starts at **01** and goes up to **N**.

Usage examples:
  python split_tsv.py --input rna_single_cell_cluster.tsv --parts 20
  python split_tsv.py --input rna_single_cell_cluster.tsv --parts 50 --outdir chunks --prefix chunk_

Notes:
- Splits rows approximately evenly across parts (block-wise).
- Header is written to each output file.
- Numeric suffix is zero-padded to the width needed (e.g., 01..20 for 20 parts, 001..050 for 50 parts).
- Prints the number of lines (including header) for each output file.
"""
import argparse
import math
from pathlib import Path

def split_tsv(input_path: str, parts: int = 20, outdir: str = 'splitted_input', prefix: str = 'part_'):
    if parts <= 0:
        raise SystemExit("--parts must be a positive integer")

    inp = Path(input_path)
    if not inp.exists():
        raise SystemExit(f"Input file not found: {inp}")

    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # First pass: read header and count data lines (exclude header)
    with inp.open('r', encoding='utf-8') as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file or missing header")
        total_rows = sum(1 for _ in f)

    # Compute rows per part (block-wise split)
    rows_per_part = math.ceil(total_rows / parts)

    # Determine zero-padding width for filenames (starting at 1)
    pad = max(2, len(str(parts)))

    # Second pass: write rows into blocks, each file gets the header
    with inp.open('r', encoding='utf-8') as f:
        header = f.readline()
        part_idx = 1  # start numbering at 1
        rows_written_in_part = 0

        def out_name(i):
            return out_dir / f"{prefix}{i:0{pad}d}.tsv"

        out_file = out_name(part_idx).open('w', encoding='utf-8')
        out_file.write(header)
        lines_in_current_file = 1  # count header line

        # Track per-file line counts
        per_file_counts = []  # list of (filename, line_count)

        for line in f:
            if rows_written_in_part >= rows_per_part and part_idx < parts:
                # finalize current file
                out_file.close()
                per_file_counts.append((out_name(part_idx).name, lines_in_current_file))
                # open next file
                part_idx += 1
                rows_written_in_part = 0
                out_file = out_name(part_idx).open('w', encoding='utf-8')
                out_file.write(header)
                lines_in_current_file = 1  # reset count with header
            out_file.write(line)
            rows_written_in_part += 1
            lines_in_current_file += 1

        out_file.close()
        per_file_counts.append((out_name(part_idx).name, lines_in_current_file))

    # Print summary with line counts
    print(f"Split completed: {parts} parts written to '{out_dir}'.")
    for fname, nlines in per_file_counts:
        print(f"{fname}: {nlines} lines")

def main():
    ap = argparse.ArgumentParser(description='Split a TSV into N parts for job arrays, keeping header in each.')
    ap.add_argument('--input', '-i', required=True, help='Input TSV file path')
    ap.add_argument('--parts', '-n', type=int, default=20, help='Number of parts to split into (default: 20)')
    ap.add_argument('--outdir', '-o', default='splitted_input', help='Output directory (default: splitted_input)')
    ap.add_argument('--prefix', default='part_', help='Output filename prefix without numeric suffix (default: part_)')
    args = ap.parse_args()
    split_tsv(args.input, parts=args.parts, outdir=args.outdir, prefix=args.prefix)

if __name__ == '__main__':
    main()

```
run it as 
```
# 50 parts, folder 'chunks', files chunk_001.tsv..chunk_050.tsv
python split_tsv.py --input rna_single_cell_cluster.tsv --parts 20 --outdir splitted_input --prefix rna_single_cell_cluster
```
