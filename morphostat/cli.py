"""
morphostat.cli
==============
"One-click" command line interface.

    morphostat run --input <folder> --group-col moa --control DMSO \
                   --batch-col Plate --source cellprofiler --outdir results
"""
from __future__ import annotations

import argparse
import sys

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="morphostat",
        description="Automated morphological post-processing for FilamentSensor / "
                    "ImageJ / CellProfiler outputs: clean -> normalize -> stats -> figures.")
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser(
        "run", help="run the full pipeline on a folder of CSVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Read all input tables, detect features, clean, normalise to controls, "
                    "test every feature, and write REPORT.md, tables and figures.")
    r.add_argument("--input", required=True, help="folder of input files, or a single file")
    r.add_argument("--outdir", default="morphostat_results",
                   help="output folder (created if absent; files with the same names are overwritten)")
    r.add_argument("--group-col", required=True, help="name of the column holding the experimental group")
    r.add_argument("--control", required=True,
                   help="value in the group column that marks control rows (case-sensitive)")
    r.add_argument("--batch-col", default=None,
                   help="plate/replicate/day column; enables per-batch normalisation to that batch's controls")
    r.add_argument("--source", default="generic",
                   choices=["generic", "cellprofiler", "filamentsensor"],
                   help="input preset; 'filamentsensor' with --fs-sample-col summarises per-filament rows per image")
    r.add_argument("--meta-cols", nargs="*", default=[],
                   help="additional columns that are NOT features (always list numeric metadata such as concentration)")
    r.add_argument("--features", nargs="*", default=None,
                   help="analyse only these columns (default: detect numeric measurement columns automatically)")
    r.add_argument("--panel", nargs="*", default=None,
                   help="features shown in the violin panel (default: the 6 with the largest effects)")
    r.add_argument("--normalization", default="robust_z_to_control",
                   choices=["robust_z_to_control", "standard", "none"],
                   help="robust_z_to_control: (x - control median) / (1.4826 x control MAD), per batch; "
                        "standard: z-score with mean and SD of all rows; none: raw values")
    r.add_argument("--impute", default="median", choices=["median", "mean", "drop", "none"],
                   help="how remaining missing values are handled (drop = remove rows with any missing value)")
    r.add_argument("--nan-frac-thresh", type=float, default=0.5,
                   help="features with a larger fraction of missing values are removed")
    r.add_argument("--alpha", type=float, default=0.05,
                   help="significance threshold for Benjamini-Hochberg-adjusted p-values (FDR)")
    r.add_argument("--fs-sample-col", default=None,
                   help="FilamentSensor only: column identifying the image/sample of each filament")
    r.add_argument("--pattern", default="*.csv", help="file-name pattern read from the input folder, e.g. '*.tsv'")
    r.add_argument("--recursive", action="store_true", help="also read files in sub-folders")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        try:
            res = run_pipeline(
                input_path=args.input, outdir=args.outdir,
                group_col=args.group_col, control=args.control, source=args.source,
                batch_col=args.batch_col, meta_cols=args.meta_cols,
                explicit_features=args.features, panel_features=args.panel,
                normalization=args.normalization, impute=args.impute,
                nan_frac_thresh=args.nan_frac_thresh, alpha=args.alpha,
                fs_sample_col=args.fs_sample_col, pattern=args.pattern,
                recursive=args.recursive,
            )
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f"\n[ERROR] {exc}", file=sys.stderr)
            return 2
        m = res["manifest"]
        print(f"\n[OK] {m['n_samples']} samples x {m['n_features']} features "
              f"-> {m['n_significant_pairwise']} significant hits (FDR<{m['alpha']})")
        print(f"     outputs written to: {args.outdir}/")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
