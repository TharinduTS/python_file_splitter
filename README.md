# python_file_splitter
Splits files into whatever number of sub files keeping column names for faster procesing

Group-aware splitting (use --group-by "Gene,Cell type,Tissue" or "Cell type,Tissue"): keeps each unit intact, so within-tissue aggregation and the bootstrap performed in recommend_min_k() are not starved of clusters or mixed across chunks. Partial groups in different files can bias effective_clusters, CV, relative Δ, and rank ρ; grouping prevents that.
Global recompute after merging filtered clusters: Your enrichment is defined as a ratio to other cell types per gene; any chunk-wise run sees only a subset of cell types/genes and can inflate/deflate denominators. Merging all *_filtered_clusters.tsv and running enrichment once globally restores the correct background across all cell types.
Header verification & manifest: ensures consistent schema and lets you audit exactly which groups landed in each part.

split_tsv_robust.py
```py

#!/usr/bin/env python3
"""
Robust TSV splitter for single-cell cluster data (and similar large TSVs).

Focus: safeguards so splitting does **not** affect downstream enrichment calculations.

Key features:
- Header preserved in every part.
- Two splitting modes:
  1) block (raw even rows per part),
  2) group-aware (keeps complete groups together in one part), with greedy balancing.
- Recommended group key for your enrichment workflow: **Gene, Cell type, Tissue**
  (or at minimum **Cell type, Tissue**) so bootstrap and aggregation remain intact within parts.
- Optional gzip input/output; auto-detect by filename suffix.
- Zero-padded numeric suffix in filenames (01..N, or wider as needed).
- Safeguards: caps parts if total rows < parts; verifies required columns; writes a manifest JSON
  with per-part row counts and group coverage; optional dry-run.
- Merge subcommand to concatenate outputs (keeping one header) and schema checks.

Usage examples:

Split by groups to preserve units relevant to enrichment (recommended):
  python split_tsv_robust.py split \
    --input rna_single_cell_cluster.tsv \
    --parts 20 \
    --group-by "Gene,Cell type,Tissue" \
    --outdir $SCRATCH/splitted_input --prefix part_ --gzip-out

Block-wise split (simple):
  python split_tsv_robust.py split --input rna_single_cell_cluster.tsv --parts 20

Dry-run (plan only, no files):
  python split_tsv_robust.py split --input rna_single_cell_cluster.tsv --parts 20 --group-by "Cell type,Tissue" --dry-run

Merge chunk outputs back together (one header):
  python split_tsv_robust.py merge \
    --pattern "$SCRATCH/enrich_parts/adjusted_*_final_enrichment.tsv" \
    --output merged_final_enrichment.tsv

Two-pass safeguard (recommended end-to-end):
  1) Run array jobs on parts to produce per-chunk "*_filtered_clusters.tsv".
  2) Merge all filtered cluster chunks: `split_tsv_robust.py merge --pattern '.../*_filtered_clusters.tsv' --output merged_filtered_clusters.tsv`
  3) Run **global** enrichment once on the merged filtered clusters to avoid partial baselines:
     `python min_clusters_and_enrichment.py --clusters merged_filtered_clusters.tsv --out-prefix global_adjusted ...`

Required columns (for your enrichment script):
  Gene, Gene name, Tissue, Cluster, Cell type, Read count, nCPM
"""
import argparse
import gzip
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# -------------------------
# Utilities
# -------------------------
def open_maybe_gzip(path: Path, mode: str = 'rt', encoding: Optional[str] = 'utf-8'):
    """Open plain text or gzip based on suffix. Accepts text ('rt','wt') or binary modes."""
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, mode)
    else:
        if 't' in mode:
            return path.open(mode, encoding=encoding)
        return path.open(mode)

def parse_group_by(group_by: Optional[str]) -> List[str]:
    if not group_by:
        return []
    cols = [c.strip() for c in group_by.split(',') if c.strip()]
    return cols

REQUIRED_COLS = ["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nCPM"]

def check_required_columns(header: str, required_cols: List[str]) -> None:
    cols = header.rstrip('\n').split('\t')
    missing = [c for c in required_cols if c not in cols]
    if missing:
        raise SystemExit(f"Input header missing required columns: {missing}. Found: {cols}")

def zero_pad_width(parts: int) -> int:
    return max(2, len(str(parts)))

def make_out_name(outdir: Path, prefix: str, idx: int, pad: int, gzip_out: bool) -> Path:
    suffix = '.tsv.gz' if gzip_out else '.tsv'
    return outdir / f"{prefix}{idx:0{pad}d}{suffix}"

# -------------------------
# Block-wise split
# -------------------------
def split_block(input_path: Path, parts: int, outdir: Path, prefix: str, gzip_out: bool,
                rows_per_part_override: Optional[int], dry_run: bool) -> Dict:
    # First pass: read header + count rows
    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file or missing header")
        total_rows = sum(1 for _ in f)
    if total_rows == 0:
        raise SystemExit("No data rows found (only header)")

    # Safeguard: cap parts to total_rows or compute from override
    if rows_per_part_override and rows_per_part_override > 0:
        parts = math.ceil(total_rows / rows_per_part_override)
    parts = min(parts, total_rows)

    rows_per_part = math.ceil(total_rows / parts)
    pad = zero_pad_width(parts)

    plan = {
        'mode': 'block',
        'parts': parts,
        'rows_per_part': rows_per_part,
        'total_rows': total_rows,
        'outdir': str(outdir),
        'prefix': prefix,
        'files': []
    }

    if dry_run:
        counts = [rows_per_part] * (parts - 1) + [total_rows - rows_per_part * (parts - 1)]
        for i, cnt in enumerate(counts, start=1):
            plan['files'].append({'name': make_out_name(outdir, prefix, i, pad, gzip_out).name,
                                  'rows_including_header': cnt + 1})
        return plan

    outdir.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        part_idx = 1
        rows_written_in_part = 0
        out_file_path = make_out_name(outdir, prefix, part_idx, pad, gzip_out)
        out_file = open_maybe_gzip(out_file_path, 'wt')
        out_file.write(header)
        lines_in_current_file = 1
        per_file_counts = []

        for line in f:
            if rows_written_in_part >= rows_per_part and part_idx < parts:
                out_file.close()
                per_file_counts.append((out_file_path.name, lines_in_current_file))
                part_idx += 1
                rows_written_in_part = 0
                out_file_path = make_out_name(outdir, prefix, part_idx, pad, gzip_out)
                out_file = open_maybe_gzip(out_file_path, 'wt')
                out_file.write(header)
                lines_in_current_file = 1
            out_file.write(line)
            rows_written_in_part += 1
            lines_in_current_file += 1
        out_file.close()
        per_file_counts.append((out_file_path.name, lines_in_current_file))

    plan['files'] = [{'name': n, 'rows_including_header': c} for n, c in per_file_counts]
    return plan

# -------------------------
# Group-aware split
# -------------------------
def split_groups(input_path: Path, parts: int, outdir: Path, prefix: str, gzip_out: bool,
                 group_cols: List[str], required_cols: List[str], dry_run: bool) -> Dict:
    # First pass: header + group sizes
    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file or missing header")
        # Verify required columns exist (so downstream enrichment won't fail)
        check_required_columns(header, required_cols)
        cols = header.rstrip('\n').split('\t')
        col_idx = {c: i for i, c in enumerate(cols)}
        for gc in group_cols:
            if gc not in col_idx:
                raise SystemExit(f"--group-by column not found in header: {gc}. Available: {cols}")
        total_rows = 0
        group_counts: Dict[Tuple, int] = {}
        for line in f:
            total_rows += 1
            fields = line.rstrip('\n').split('\t')
            key = tuple(fields[col_idx[c]] for c in group_cols)
            group_counts[key] = group_counts.get(key, 0) + 1
    if total_rows == 0:
        raise SystemExit("No data rows found (only header)")

    # Safeguard: cap parts to min(total_rows, number_of_groups)
    num_groups = len(group_counts)
    parts = min(parts, max(1, num_groups))

    # Greedy bin packing: largest groups first to balance rows across parts
    sorted_groups = sorted(group_counts.items(), key=lambda kv: kv[1], reverse=True)
    bins = [{'rows': 0, 'groups': []} for _ in range(parts)]
    for key, gsize in sorted_groups:
        target = min(range(parts), key=lambda i: bins[i]['rows'])
        bins[target]['rows'] += gsize
        bins[target]['groups'].append(key)

    pad = zero_pad_width(parts)

    plan = {
        'mode': 'group',
        'parts': parts,
        'total_rows': total_rows,
        'num_groups': num_groups,
        'group_cols': group_cols,
        'outdir': str(outdir),
        'prefix': prefix,
        'files': [{'name': make_out_name(outdir, prefix, i+1, pad, gzip_out).name,
                   'assigned_groups': len(bins[i]['groups']),
                   'planned_rows_excluding_header': bins[i]['rows']}
                  for i in range(parts)]
    }

    if dry_run:
        return plan

    outdir.mkdir(parents=True, exist_ok=True)

    # Second pass: route each line to its assigned bin
    writers = []
    for i in range(parts):
        p = make_out_name(outdir, prefix, i+1, pad, gzip_out)
        fh = open_maybe_gzip(p, 'wt')
        writers.append({'path': p, 'fh': fh, 'count': 0})

    group_to_bin: Dict[Tuple, int] = {}
    for i in range(parts):
        for g in bins[i]['groups']:
            group_to_bin[g] = i

    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        for w in writers:
            w['fh'].write(header)
            w['count'] = 1
        cols = header.rstrip('\n').split('\t')
        col_idx = {c: i for i, c in enumerate(cols)}
        for line in f:
            fields = line.rstrip('\n').split('\t')
            key = tuple(fields[col_idx[c]] for c in group_cols)
            bin_idx = group_to_bin.get(key)
            if bin_idx is None:
                raise SystemExit(f"Internal error: group {key} not assigned to any bin")
            w = writers[bin_idx]
            w['fh'].write(line)
            w['count'] += 1

    # Close writers and collect counts (include assigned_groups)
    per_file_infos = []
    for i, w in enumerate(writers):
        w['fh'].close()
        per_file_infos.append({
            'name': w['path'].name,
            'rows_including_header': w['count'],
            'assigned_groups': len(bins[i]['groups'])
        })

    plan['files'] = per_file_infos
    return plan

# -------------------------
# Merge helper (concatenate with single header)
# -------------------------
def merge_files(pattern: str, output: Path, gzip_out: bool = False) -> Dict:
    import glob
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No files matched pattern: {pattern}")

    def read_header(fp: str) -> str:
        p = Path(fp)
        with open_maybe_gzip(p, 'rt') as f:
            hdr = f.readline()
            if not hdr:
                raise SystemExit(f"Empty file or missing header: {fp}")
            return hdr

    main_header = read_header(files[0])

    # Verify all headers match exactly
    for fp in files[1:]:
        hdr = read_header(fp)
        if hdr.rstrip('\n') != main_header.rstrip('\n'):
            raise SystemExit(f"Header mismatch between {files[0]} and {fp}\n{main_header}\n!=\n{hdr}")

    out_path = output
    if gzip_out and not str(out_path).endswith('.gz'):
        out_path = Path(str(out_path) + '.gz')

    with open_maybe_gzip(out_path, 'wt') as out:
        out.write(main_header)
        total_lines = 1
        for fp in files:
            p = Path(fp)
            with open_maybe_gzip(p, 'rt') as f:
                _ = f.readline()  # skip header
                for line in f:
                    out.write(line)
                    total_lines += 1

    return {'output': str(out_path), 'lines_including_header': total_lines, 'merged_files': files}

# -------------------------
# CLI
# -------------------------
def main():
    ap = argparse.ArgumentParser(description='Robust TSV splitter/merger with safeguards for enrichment workflows')
    sub = ap.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('split', help='Split a large TSV into parts')
    sp.add_argument('--input', '-i', required=True, help='Input TSV file path (.tsv or .tsv.gz)')
    sp.add_argument('--parts', '-n', type=int, default=20, help='Number of parts to split into (default: 20)')
    sp.add_argument('--rows-per-part', type=int, default=None, help='Override rows per part (block mode only)')
    sp.add_argument('--outdir', '-o', default='splitted_input', help='Output directory (default: splitted_input)')
    sp.add_argument('--prefix', default='part_', help='Output filename prefix (default: part_)')
    sp.add_argument('--gzip-out', action='store_true', help='Write gzip-compressed parts (.tsv.gz)')
    sp.add_argument('--group-by', default='', help='Comma-separated column names to keep rows grouped (e.g., "Gene,Cell type,Tissue")')
    sp.add_argument('--require-cols', default=','.join(REQUIRED_COLS), help='Comma-separated required columns to verify in header')
    sp.add_argument('--manifest', default='split_manifest.json', help='Where to write JSON manifest of the split plan/results')
    sp.add_argument('--dry-run', action='store_true', help='Plan only; do not write files')

    mp = sub.add_parser('merge', help='Concatenate TSVs (one header)')
    mp.add_argument('--pattern', required=True, help='Glob pattern for input files to merge')
    mp.add_argument('--output', required=True, help='Output merged file path (.tsv or .tsv.gz)')
    mp.add_argument('--gzip-out', action='store_true', help='Write gzip-compressed output (adds .gz if missing)')

    args = ap.parse_args()

    if args.cmd == 'split':
        input_path = Path(args.input)
        if not input_path.exists():
            raise SystemExit(f"Input file not found: {input_path}")
        outdir = Path(args.outdir)
        required_cols = [c.strip() for c in args.require_cols.split(',') if c.strip()]
        group_cols = parse_group_by(args.group_by)

        # Choose mode
        if group_cols:
            plan = split_groups(input_path, parts=args.parts, outdir=outdir, prefix=args.prefix,
                                gzip_out=args.gzip_out, group_cols=group_cols,
                                required_cols=required_cols, dry_run=args.dry_run)
        else:
            plan = split_block(input_path, parts=args.parts, outdir=outdir, prefix=args.prefix,
                               gzip_out=args.gzip_out, rows_per_part_override=args.rows_per_part, dry_run=args.dry_run)

        # Write manifest
        manifest_path = Path(args.manifest)
        data = {
            'input': str(input_path),
            'mode': plan['mode'],
            'parts': plan['parts'],
            'outdir': plan['outdir'],
            'prefix': plan['prefix'],
            'total_rows': plan.get('total_rows'),
            'num_groups': plan.get('num_groups'),
            'group_cols': plan.get('group_cols'),
            'files': plan['files']
        }
        with manifest_path.open('w', encoding='utf-8') as mf:
            json.dump(data, mf, indent=2)
        print(f"Split {'planned' if args.dry_run else 'completed'}: {plan['parts']} parts to '{outdir}'. Manifest: {manifest_path}")
        for info in plan['files']:
            name = info['name']
            rows = info.get('rows_including_header', info.get('planned_rows_excluding_header', 0) + 1)
            if plan['mode'] == 'group' and 'assigned_groups' in info:
                print(f"{name}: {rows} lines (groups: {info['assigned_groups']})")
            else:
                print(f"{name}: {rows} lines")

    elif args.cmd == 'merge':
        output = Path(args.output)
        res = merge_files(args.pattern, output, gzip_out=args.gzip_out)
        print(f"Merged {len(res['merged_files'])} files into: {res['output']} ({res['lines_including_header']} lines)")

if __name__ == '__main__':
    main()

```
run it as 
```
python split_tsv_robust.py split \
  --input rna_single_cell_cluster.tsv \
  --parts 20 \
  --group-by "Gene,Cell type,Tissue" \
  --outdir splitted_input --prefix part_ --gzip-out \
  --manifest split_manifest.json
```
Then I am running the following enrichment calculation script for all parts

min_clusters_and_enrichment.py
```py
#!/usr/bin/env python3
"""
min_clusters_and_enrichment.py

Purpose:
  - Use ONLY rna_single_cell_cluster.tsv (columns: Gene, Gene name, Tissue, Cluster, Cell type, Read count, nCPM).
  - Recommend the minimum number of clusters per (Cell type, Tissue) for stable enrichment (weighted by Read count).
  - Apply these replication adjustments (exclude unstable pairs) and recompute proper enrichment values.

Outputs:
  - min_clusters_per_celltype.tsv : recommendations per (Cell type, Tissue) with diagnostics.
  - filtered_clusters.tsv         : cluster-level rows retained for enrichment after applying recommendations.
  - final_enrichment.tsv          : gene × cell type table with adjusted enrichment.

CLI example:
  python min_clusters_and_enrichment.py \
    --clusters rna_single_cell_cluster.tsv \
    --out-prefix adjusted \
    --B 100 --cv 0.20 --reldelta 0.10 --rho 0.90 --min-k 2 \
    --tissue-weighting weighted --pseudocount 0.01
"""

# --- Method notes (short) ---
# Weighted within‑tissue aggregation:
#   Combine cluster nCPM using weights (Read count) so large clusters contribute more than tiny ones.
#   Formula: weighted_mean = Σ(nCPM_i * w_i) / Σ(w_i). Prevents small/noisy clusters skewing the cell type profile.
#
# Weighted bootstrap (example):
#   Sampling clusters with probability ∝ Read count. If clusters have Read counts [1000, 200, 50],
#   then p = [0.8, 0.16, 0.04]. For k=2, most samples include the large cluster + one smaller —
#   resamples reflect the true population.
#
# CV of enrichment (target ≤ 0.20):
#   CV = SD(enrichment) / Mean(enrichment). Values ≤ 0.20 indicate low relative variability → stable.
#   Example: values [10,12,8,10] → mean=10, SD≈1.63 → CV≈0.163 (stable).
#
# Median relative change vs baseline (target ≤ 0.10):
#   median(|E_boot - E_base| / |E_base|). ≤ 10% means resampled enrichment stays close to baseline.
#   Example: base=[10,12,8], boot=[9,11,8] → changes=[0.10,0.083,0] → median=0.083 (stable).
#
# Spearman rank correlation vs baseline (target ≥ 0.90):
#   Checks if gene rank order is preserved. ρ ≥ 0.90 ⇒ very similar ranking.
#   Example: ranks differ slightly → ρ≈0.95 (stable).
# ---

import argparse
import numpy as np
import pandas as pd

# ---------------------------
# Utilities
# ---------------------------

def spearman_corr(x, y):
    """Spearman correlation without SciPy: rank then Pearson."""
    xr = pd.Series(x).rank(method="average").values
    yr = pd.Series(y).rank(method="average").values
    if len(xr) < 2:
        return np.nan
    c = np.corrcoef(xr, yr)
    return c[0, 1]


def effective_clusters(weights):
    """Effective number of clusters given weights (e.g., Read count)."""
    w = np.asarray(weights, dtype=float)
    ss = (w ** 2).sum()
    s = w.sum()
    return (s ** 2 / ss) if ss > 0 else 0.0


# ---------------------------
# Aggregation & Enrichment (HPA-consistent)
# ---------------------------

# Weighted within‑tissue aggregation using Read count as weight

def aggregate_within_tissue(df, expr_col="nCPM", weight_col="Read count"):
    """
    Weighted mean per (Gene, Gene name, Cell type, Tissue) using Read count as weight.
    Also computes clusters_used_tissue and weight_sum_tissue for diagnostics.
    """
    df = df.copy()
    df["_w"] = pd.to_numeric(df[weight_col], errors="coerce").fillna(1.0)
    df["_val"] = pd.to_numeric(df[expr_col], errors="coerce")

    gcols = ["Gene", "Gene name", "Cell type", "Tissue"]

    # Vectorized weighted sum per group (no groupby.apply)
    df["__prod"] = df["_w"] * df["_val"]
    num = df.groupby(gcols)["__prod"].sum()
    den = df.groupby(gcols)["_w"].sum()

    out = (num / den).reset_index(name="avg_nCPM")

    # diagnostics
    out["clusters_used_tissue"] = df.groupby(gcols)["Cluster"].nunique().values
    out["weight_sum_tissue"] = den.values

    # cleanup temp column
    df.drop(columns="__prod", inplace=True)

    return out


def integrate_across_tissues(df_ct, how="weighted"):
    """
    Final cell type profile by averaging tissue-level profiles.
    - 'weighted': weighted mean by weight_sum_tissue
    - 'unweighted': simple mean
    Also returns diagnostics (#tissues, total clusters across tissues, total weight).
    """
    df = df_ct.copy()
    gcols = ["Gene", "Gene name", "Cell type"]
    if how == "weighted" and ("weight_sum_tissue" in df.columns):
        # Vectorized weighted sum across tissues (no groupby.apply)
        df["__prod_ct"] = df["avg_nCPM"] * df["weight_sum_tissue"]
        num = df.groupby(gcols)["__prod_ct"].sum()
        den = df.groupby(gcols)["weight_sum_tissue"].sum()
        out = (num / den).reset_index(name="avg_nCPM")
        df.drop(columns="__prod_ct", inplace=True)
    else:
        out = (df.groupby(gcols, as_index=False)
                 .agg(avg_nCPM=("avg_nCPM", "mean")))

    diag = (df.groupby(gcols, as_index=False)
              .agg(datasets_used=("Tissue", "nunique"),
                   clusters_used=("clusters_used_tissue", "sum"),
                   total_weight=("weight_sum_tissue", "sum")))
    out = out.merge(diag, on=gcols, how="left")
    return out


def add_enrichment(agg_df, gene_col="Gene", value_col="avg_nCPM",
                   out_col="Enrichment score",
                   min_background=1e-3, min_expression=0.0,
                   pseudocount=None):
    """
    Enrichment = value / mean(value in other cell types of the same gene), with safeguards.
    """
    df = agg_df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    sums = df.groupby(gene_col)[value_col].transform("sum")
    counts = df.groupby(gene_col)[value_col].transform("count")
    denom_counts = counts - 1
    avg_other = (sums - df[value_col]) / denom_counts
    avg_other = avg_other.mask(denom_counts <= 0, np.nan)

    if pseudocount is not None:
        avg_other = avg_other + pseudocount
        # Optionally stabilize numerator: df[value_col] = df[value_col] + pseudocount

    denom = np.maximum(avg_other, min_background)
    numer = df[value_col].where(df[value_col] >= min_expression, np.nan)

    df[out_col] = np.divide(numer, denom,
                            out=np.full(df.shape[0], np.nan),
                            where=(denom > 0))
    df.loc[avg_other.isna(), out_col] = np.nan
    return df


# ---------------------------
# Baseline build
# ---------------------------

def build_baseline(cluster_df, expr_col="nCPM", tissue_weighting="weighted", pseudocount=None):
    """
    Build baseline aggregated profiles and enrichment from full dataset.
    Returns:
      within_tissue (DataFrame), across_tissue (DataFrame), baseline_enrich (DataFrame), base_map (dict)
    """
    within_tissue = aggregate_within_tissue(cluster_df, expr_col=expr_col, weight_col="Read count")
    across_tissue = integrate_across_tissues(within_tissue, how=tissue_weighting)
    baseline_enrich = add_enrichment(across_tissue, pseudocount=pseudocount)

    baseline_enrich["key"] = baseline_enrich["Gene"].astype(str) + "||" + baseline_enrich["Cell type"].astype(str)
    base_map = dict(zip(baseline_enrich["key"], baseline_enrich["Enrichment score"]))
    return within_tissue, across_tissue, baseline_enrich, base_map


# ---------------------------
# Leave-One-Tissue-Out (dataset guard)
# ---------------------------

def lodo_stability(cluster_df, base_map, expr_col="nCPM", tissue_weighting="weighted", pseudocount=None):
    """
    Compute per-cell-type max relative change in enrichment when leaving one Tissue out.
    Returns dict: Cell type -> max_relative_change (lower is better).
    """
    df = cluster_df.copy()
    tissues = df["Tissue"].dropna().unique().tolist()
    ct_max_rel = {}

    if not tissues:
        return ct_max_rel

    for t in tissues:
        sub = df[df["Tissue"] != t]
        if sub.empty:
            continue
        within = aggregate_within_tissue(sub, expr_col=expr_col, weight_col="Read count")
        across = integrate_across_tissues(within, how=tissue_weighting)
        enr = add_enrichment(across, pseudocount=pseudocount)
        enr["key"] = enr["Gene"].astype(str) + "||" + enr["Cell type"].astype(str)
        enr["base"] = enr["key"].map(base_map)
        enr = enr.dropna(subset=["Enrichment score", "base"])
        if enr.empty:
            continue

        rel = np.abs(enr["Enrichment score"] - enr["base"]) / (np.abs(enr["base"]) + 1e-9)

        for ct, grp in enr.groupby("Cell type"):
            mx = rel.loc[grp.index].max()
            prev = ct_max_rel.get(ct, 0.0)
            ct_max_rel[ct] = max(prev, mx)

    return ct_max_rel


# ---------------------------
# Weighted bootstrap for min k
# ---------------------------

# Weighted bootstrap stability to recommend minimal clusters k

def recommend_min_k(cluster_df,
                    B=100,
                    cv_thresh=0.20,
                    rel_thresh=0.10,
                    rank_r_thresh=0.90,
                    min_k_rule=2,
                    expr_col="nCPM",
                    tissue_weighting="weighted",
                    pseudocount=None):
    """
    For each (Cell type, Tissue), recommend minimal number of clusters k that yields stable enrichment.
    Sampling is within the (Cell type, Tissue) cluster set. The rest of the dataset stays intact.
    """
    df = cluster_df.copy()
    within_base, across_base, enr_base, base_map = build_baseline(df, expr_col=expr_col,
                                                                  tissue_weighting=tissue_weighting,
                                                                  pseudocount=pseudocount)

    # weights: Read count per (Cell type, Tissue, Cluster)
    wdf = (df.groupby(["Cell type", "Tissue", "Cluster"])['Read count']
             .sum().reset_index())
    weight_map = {(r['Cell type'], r['Tissue'], r['Cluster']): r['Read count'] for _, r in wdf.iterrows()}

    # dataset/tissue counts per cell type
    ct_tissues_used = df.groupby("Cell type")["Tissue"].nunique().to_dict()

    # LODO guard
    lodo = lodo_stability(df, base_map, expr_col=expr_col,
                          tissue_weighting=tissue_weighting,
                          pseudocount=pseudocount)

    results = []

    # Iterate per (Cell type, Tissue)
    for (ct, tissue), clusters in df.groupby(["Cell type", "Tissue"])["Cluster"].unique().items():
        clusters = list(clusters)
        weights = [weight_map.get((ct, tissue, c), 1.0) for c in clusters]
        neff = effective_clusters(weights)
        k_max = len(clusters)
        ds_used = ct_tissues_used.get(ct, np.nan)
        lodo_max = lodo.get(ct, np.nan)

        # Pre-filters
        if k_max < min_k_rule or neff < 2:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": np.nan,
                "reason": f"Insufficient replication (k<{min_k_rule} or neff<2)",
                "median_cv": np.nan,
                "median_rel_delta": np.nan,
                "median_rank_corr": np.nan
            })
            continue

        # Dataset guard: prefer ≥2 tissues or LODO max change ≤ 0.15
        dataset_guard = (not np.isnan(ds_used)) and (ds_used >= 2)
        lodo_guard = (not np.isnan(lodo_max)) and (lodo_max <= 0.15)

        # Candidate k
        k_choice = None
        med_cv = med_rel = med_r = np.nan

        # Weighted sampling without replacement within this (ct, tissue)
        w = np.asarray(weights, dtype=float)
        p = (w / w.sum()) if w.sum() > 0 else None

        for k in range(min_k_rule, k_max + 1):
            cvs, rels, rhos = [], [], []

            for b in range(B):
                idxs = np.random.choice(np.arange(k_max), size=k, replace=False, p=p)
                chosen = [clusters[i] for i in idxs]

                # Build a bootstrapped dataset:
                sub_ct_tissue = df[(df["Cell type"] == ct) & (df["Tissue"] == tissue) & (df["Cluster"].isin(chosen))]
                other = df[~((df["Cell type"] == ct) & (df["Tissue"] == tissue))]
                boot_df = pd.concat([other, sub_ct_tissue], ignore_index=True)

                within_b = aggregate_within_tissue(boot_df, expr_col=expr_col, weight_col="Read count")
                across_b = integrate_across_tissues(within_b, how=tissue_weighting)
                enr_b = add_enrichment(across_b, pseudocount=pseudocount)

                # Align with baseline for this cell type only
                e_ct = enr_b[enr_b["Cell type"] == ct].copy()
                e_ct["key"] = e_ct["Gene"].astype(str) + "||" + e_ct["Cell type"].astype(str)
                e_ct["base"] = e_ct["key"].map(base_map)
                e_ct = e_ct.dropna(subset=["Enrichment score", "base"])
                if e_ct.empty:
                    continue

                vals = e_ct["Enrichment score"].values
                mu = np.mean(vals)
                sd = np.std(vals, ddof=1)
                cvs.append(sd / mu if mu > 0 else np.inf)

                rel = np.median(np.abs(vals - e_ct["base"].values) / (np.abs(e_ct["base"].values) + 1e-9))
                rels.append(rel)

                rhos.append(spearman_corr(vals, e_ct["base"].values))

            med_cv = np.median(cvs) if len(cvs) else np.inf
            med_rel = np.median(rels) if len(rels) else np.inf
            med_r = np.median(rhos) if len(rhos) else 0.0

            if (med_cv <= cv_thresh) and (med_rel <= rel_thresh) and (med_r >= rank_r_thresh):
                # If dataset guards fail, you may choose to require one higher k for conservatism
                k_choice = k if (dataset_guard or lodo_guard) else min(k_max, k + 1)
                break

        if k_choice is None:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": k_max,
                "reason": "No k met thresholds; using max_k (flag for caution)",
                "median_cv": (None if np.isinf(med_cv) else round(float(med_cv), 3)),
                "median_rel_delta": (None if np.isinf(med_rel) else round(float(med_rel), 3)),
                "median_rank_corr": (None if np.isnan(med_r) else round(float(med_r), 3))
            })
        else:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": k_choice,
                "reason": "Meets stability thresholds" + ("" if (dataset_guard or lodo_guard) else " (dataset guard raised k by +1)"),
                "median_cv": round(float(med_cv), 3),
                "median_rel_delta": round(float(med_rel), 3),
                "median_rank_corr": round(float(med_r), 3)
            })

    return pd.DataFrame(results).sort_values(["Cell type", "Tissue"])


# ---------------------------
# Apply adjustments and recompute enrichment
# ---------------------------

# Apply replication filters and recompute final enrichment

def apply_adjustments_and_enrich(cluster_df, recs_df,
                                 expr_col="nCPM",
                                 tissue_weighting="weighted",
                                 pseudocount=None):
    """
    Filter cluster_df to include only (Cell type, Tissue) pairs that meet replication:
      - available_clusters >= recommended_min_k
      - effective_clusters >= 2
      - recommended_min_k is not NaN
    Then recompute within-tissue aggregation, across-tissue integration, and enrichment.
    Returns: filtered_clusters, final_enrichment
    """
    # Merge recommendations on (Cell type, Tissue)
    key_cols = ["Cell type", "Tissue"]
    keep = recs_df.dropna(subset=["recommended_min_k"]).copy()
    keep = keep[(keep["available_clusters"] >= keep["recommended_min_k"]) & (keep["effective_clusters"] >= 2)]

    if keep.empty:
        # No pairs meet replication; return empty results
        return cluster_df.iloc[0:0].copy(), pd.DataFrame(columns=["Gene","Gene name","Cell type","avg_nCPM","datasets_used","clusters_used","total_weight","Enrichment score"]) 

    allowed = keep[key_cols].drop_duplicates()
    filt = cluster_df.merge(allowed, on=key_cols, how="inner")

    # Recompute enrichment on filtered clusters
    within = aggregate_within_tissue(filt, expr_col=expr_col, weight_col="Read count")
    across = integrate_across_tissues(within, how=tissue_weighting)
    final_enr = add_enrichment(across, pseudocount=pseudocount)

    return filt, final_enr


# ---------------------------
# CLI
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Recommend minimum clusters and recompute enrichment using only rna_single_cell_cluster.tsv")
    ap.add_argument("--clusters", required=True, help="Input TSV: rna_single_cell_cluster.tsv")
    ap.add_argument("--out-prefix", default="adjusted", help="Output file prefix")
    ap.add_argument("--B", type=int, default=100, help="Bootstrap replicates")
    ap.add_argument("--cv", type=float, default=0.20, help="CV threshold")
    ap.add_argument("--reldelta", type=float, default=0.10, help="Relative Δ threshold")
    ap.add_argument("--rho", type=float, default=0.90, help="Spearman rank correlation threshold")
    ap.add_argument("--min-k", type=int, default=2, help="Minimum k to consider")
    ap.add_argument("--tissue-weighting", choices=["weighted","unweighted"], default="weighted", help="Across-tissue integration weighting")
    ap.add_argument("--pseudocount", type=float, default=0.01, help="Optional pseudocount for enrichment denominator")

    args = ap.parse_args()

    # Load clusters file
    df = pd.read_csv(args.clusters, sep='\t')
    required = ["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nCPM"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise SystemExit(f"Input is missing required columns: {miss}")

    # Recommend min k
    recs = recommend_min_k(df,
                           B=args.B,
                           cv_thresh=args.cv,
                           rel_thresh=args.reldelta,
                           rank_r_thresh=args.rho,
                           min_k_rule=args.min_k,
                           expr_col="nCPM",
                           tissue_weighting=args.tissue_weighting,
                           pseudocount=args.pseudocount)

    # Apply adjustments & recompute enrichment
    filt_clusters, final_enrichment = apply_adjustments_and_enrich(df, recs,
                                                                  expr_col="nCPM",
                                                                  tissue_weighting=args.tissue_weighting,
                                                                  pseudocount=args.pseudocount)

    # Save outputs
    recs_out = f"{args.out_prefix}_min_clusters_per_celltype.tsv"
    filt_out = f"{args.out_prefix}_filtered_clusters.tsv"
    enr_out  = f"{args.out_prefix}_final_enrichment.tsv"

    recs.to_csv(recs_out, sep='\t', index=False)
    filt_clusters.to_csv(filt_out, sep='\t', index=False)
    final_enrichment.to_csv(enr_out, sep='\t', index=False)

    # Console summary
    n_pairs_total = df.groupby(["Cell type","Tissue"]).ngroups
    n_pairs_kept  = filt_clusters.groupby(["Cell type","Tissue"]).ngroups if not filt_clusters.empty else 0
    print(f"Saved: {recs_out}\nSaved: {filt_out}\nSaved: {enr_out}")
    print(f"Pairs total: {n_pairs_total}, kept after adjustments: {n_pairs_kept}")


if __name__ == "__main__":
    main()
```
Run this as a job array
```
#!/bin/bash
#SBATCH --job-name=enrich_array
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-20%5
#SBATCH --output=logs/enrich_%A_%a.out
#SBATCH --error=logs/enrich_%A_%a.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# ---------------------------
# Configurable parameters (defaults target the current directory)
# Override at submit-time if needed, e.g.:
#   PARTS=50 INDIR=./chunks PREFIX=chunk_ OUTDIR=./parts_out sbatch enrich_array.sbatch
# ---------------------------
PARTS=${PARTS:-20}                # total number of chunks (must match your split)
INDIR=${INDIR:-./splitted_input}  # directory holding part_XX.tsv(.gz) relative to CWD
PREFIX=${PREFIX:-part_}           # chunk filename prefix (e.g., part_01.tsv.gz)
OUTDIR=${OUTDIR:-./enrich_parts}  # where per-chunk outputs will be written (relative to CWD)

# ---------------------------
# Environment
# ---------------------------
module load gcc arrow
module load python

python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate

mkdir -p "$OUTDIR" logs

# ---------------------------
# Map array index -> chunk file
# ---------------------------
PADDED=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

# Prefer .tsv.gz, fall back to .tsv if .gz doesn’t exist
INPUT_GZ="${INDIR}/${PREFIX}${PADDED}.tsv.gz"
INPUT_TSV="${INDIR}/${PREFIX}${PADDED}.tsv"
if [ -s "$INPUT_GZ" ]; then
  INPUT="$INPUT_GZ"
elif [ -s "$INPUT_TSV" ]; then
  INPUT="$INPUT_TSV"
else
  echo "[array] Missing input for index ${SLURM_ARRAY_TASK_ID}:"
  echo "  tried: $INPUT_GZ and $INPUT_TSV" >&2
  exit 2
fi

PREFIX_OUT="${OUTDIR}/adjusted_${PADDED}"

echo "[array] task ${SLURM_ARRAY_TASK_ID}/${PARTS}"
echo "[array] input:      $INPUT"
echo "[array] out-prefix: $PREFIX_OUT"

# ---------------------------
# Run enrichment script per chunk
# ---------------------------
python min_clusters_and_enrichment.py \
  --clusters "$INPUT" \
  --out-prefix "$PREFIX_OUT" \
  --B 100 \
  --cv 0.20 \
  --reldelta 0.10 \
  --rho 0.90 \
  --min-k 2 \
  --tissue-weighting weighted \
  --pseudocount 0.01
```
















merge
```
python split_tsv_robust.py merge \
  --pattern "$SCRATCH/enrich_parts/adjusted_*_filtered_clusters.tsv" \
  --output $SCRATCH/merged_filtered_clusters.tsv
```
global enrichment recompute on the merged filtered clusters (to avoid any chunk-wise denominator bias):
```

python min_clusters_and_enrichment.py \
  --clusters $SCRATCH/merged_filtered_clusters.tsv \
  --out-prefix $SCRATCH/global_adjusted \
  --B 100 --cv 0.20 --reldelta 0.10 --rho 0.90 --min-k 2 \
  --tissue-weighting weighted --pseudocount 0.01
```
