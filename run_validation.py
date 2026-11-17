"""
run_validation.py
=================
Reproduce, end to end, the validation reported in the Bio-protocol article
"MorphoStat: An Automated Pipeline for Normalization, Assumption-Aware
Statistics, and Publication-Ready Visualization of Cell Morphology Profiles".

What it does
------------
1. Obtains the BBBC021 per-well morphological profiles (MCF-7 cells; Ljosa et al.
   2013) from the Broad Institute ``cytominergallery`` repository, unless a local
   copy is supplied with --extdata.
2. Builds one CSV per plate with ``prepare_bbbc021.py`` (632 wells, 55 plates).
3. Writes two analysis subsets and runs the MorphoStat pipeline on each:
     validation_results/two_group/     DMSO vs Actin disruptors      (345 wells)
     validation_results/multi_group/   DMSO + 4 cytoskeletal classes (450 wells)
4. Performs the independent raw-data direction check on cell area
   (matched-plate medians, no normalisation, no statistics).
5. Records the software versions and machine details used for the run.

This script is the supported way to regenerate every figure, table and
supplementary file of the article (Figure 7, direction_check.txt and
environment.txt are produced only by this script).

Usage (identical on Linux, macOS and Windows)
---------------------------------------------
    python run_validation.py
    python run_validation.py --extdata path/to/cytominergallery/inst/extdata
    python run_validation.py --outdir my_validation

Requires git on the PATH for the automatic download. Without git, download
https://github.com/cytomining/cytominergallery/archive/refs/heads/master.zip,
unzip it, and pass its inst/extdata folder with --extdata.
"""
from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import morphostat  # noqa: E402
from morphostat import viz as mviz  # noqa: E402
from morphostat import run_pipeline  # noqa: E402
from prepare_bbbc021 import main as prepare_bbbc021  # noqa: E402

GALLERY_URL = "https://github.com/cytomining/cytominergallery.git"
CYTOSKELETAL = ["DMSO", "Actin disruptors", "Microtubule stabilizers",
                "Microtubule destabilizers", "Aurora kinase inhibitors"]
# Metadata columns written by prepare_bbbc021.py. "concentration" is numeric, so it
# must be declared as metadata; otherwise it would be analysed as a feature.
META = ["Well", "compound", "concentration"]
# The eight CellProfiler AreaShape descriptors discussed in the article (Figure 1, Table 3).
SHAPE = ["Cells_AreaShape_Area", "Cells_AreaShape_Perimeter", "Cells_AreaShape_FormFactor",
         "Cells_AreaShape_Solidity", "Cells_AreaShape_Compactness",
         "Cells_AreaShape_MajorAxisLength", "Cells_AreaShape_MinorAxisLength",
         "Cells_AreaShape_Eccentricity"]
COMMON = dict(group_col="moa", control="DMSO", batch_col="Plate",
              source="cellprofiler",
              normalization="robust_z_to_control", impute="median",
              nan_frac_thresh=0.5, alpha=0.05)


def get_extdata(extdata: str | None, workdir: str) -> str:
    if extdata:
        return extdata
    clone = os.path.join(workdir, "cytominergallery")
    if not os.path.isdir(clone):
        if shutil.which("git") is None:
            sys.exit("git was not found. Download the cytominergallery zip (see the "
                     "header of this script) and re-run with --extdata.")
        print(f"[1/5] Downloading BBBC021 profiles (about 10 MB) from {GALLERY_URL} ...")
        subprocess.run(["git", "clone", "--depth", "1", GALLERY_URL, clone], check=True)
    return os.path.join(clone, "inst", "extdata")


def write_subset(src_dir: str, dst_dir: str, groups: list[str]) -> None:
    os.makedirs(dst_dir, exist_ok=True)
    for f in sorted(os.listdir(src_dir)):
        if f.endswith(".csv"):
            d = pd.read_csv(os.path.join(src_dir, f))
            d = d[d["moa"].isin(groups)]
            if len(d):
                d.to_csv(os.path.join(dst_dir, f), index=False)


def raw_direction_check(input_dir: str, feature: str = "Cells_AreaShape_Area") -> str:
    frames = [pd.read_csv(os.path.join(input_dir, f))
              for f in sorted(os.listdir(input_dir)) if f.endswith(".csv")]
    d = pd.concat(frames, ignore_index=True)
    lines, diffs = [], []
    for plate, g in d.groupby("Plate"):
        t = g.loc[g["moa"] == "Actin disruptors", feature]
        c = g.loc[g["moa"] == "DMSO", feature]
        if len(t) and len(c):
            diff = t.median() - c.median()
            diffs.append(diff)
            lines.append(f"  {plate:<12s} treated={t.median():9.1f}  DMSO={c.median():9.1f}  diff={diff:9.1f}")
    lower = sum(x < 0 for x in diffs)
    head = (f"Raw-data direction check ({feature}, un-normalised, matched plates)\n"
            f"  plates with both groups: {len(diffs)}; treated median lower on {lower}/{len(diffs)}\n"
            f"  median within-plate difference: {pd.Series(diffs).median():.0f} px^2\n")
    return head + "\n".join(lines) + "\n"


def environment() -> str:
    import matplotlib, numpy, scipy, seaborn, sklearn
    try:
        import statsmodels
        sm = statsmodels.__version__
    except Exception:  # pragma: no cover
        sm = "not installed"
    return "\n".join([
        f"MorphoStat   {morphostat.__version__}",
        f"Python       {platform.python_version()}",
        f"NumPy        {numpy.__version__}",
        f"pandas       {pd.__version__}",
        f"SciPy        {scipy.__version__}",
        f"scikit-learn {sklearn.__version__}",
        f"statsmodels  {sm}",
        f"Matplotlib   {matplotlib.__version__}",
        f"seaborn      {seaborn.__version__}",
        f"OS           {platform.platform()}",
        f"CPU cores    {os.cpu_count()}",
    ]) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Reproduce the MorphoStat BBBC021 validation.")
    ap.add_argument("--extdata", default=None,
                    help="existing cytominergallery/inst/extdata folder (skips download)")
    ap.add_argument("--outdir", default="validation_results", help="output folder")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    extdata = get_extdata(a.extdata, a.outdir)

    print("[2/5] Building per-plate input tables ...")
    all_dir = os.path.join(a.outdir, "input_all_plates")
    prepare_bbbc021(extdata, all_dir)

    two_in = os.path.join(a.outdir, "input_two_group")
    multi_in = os.path.join(a.outdir, "input_multi_group")
    write_subset(all_dir, two_in, ["DMSO", "Actin disruptors"])
    write_subset(all_dir, multi_in, CYTOSKELETAL)

    print("[3/5] Two-group analysis: DMSO vs Actin disruptors ...")
    r2 = run_pipeline(two_in, os.path.join(a.outdir, "two_group"), meta_cols=META,
                      panel_features=SHAPE, **COMMON)
    print("[4/5] Multi-group analysis: DMSO + four cytoskeletal classes ...")
    rm = run_pipeline(multi_in, os.path.join(a.outdir, "multi_group"), meta_cols=META,
                      panel_features=SHAPE, **COMMON)
    # Figure 4: Spearman correlation among the eight shape descriptors (all 450 wells)
    mviz.correlation_heatmap(rm["normalized"], SHAPE,
                             os.path.join(a.outdir, "multi_group", "figures"),
                             name="shape_correlation_spearman", method="spearman",
                             title="Shape descriptor correlation (Spearman)")

    print("[5/5] Raw-data direction check and environment record ...")
    check = raw_direction_check(two_in)
    with open(os.path.join(a.outdir, "direction_check.txt"), "w") as fh:
        fh.write(check)
    with open(os.path.join(a.outdir, "environment.txt"), "w") as fh:
        fh.write(environment())

    pw, om = r2["pairwise"], rm["omnibus"]
    print("\n================ VALIDATION SUMMARY ================")
    print(f"Two-group : {r2['manifest']['n_samples']} wells x {r2['manifest']['n_features']} features; "
          f"{int(pw['significant'].sum())} significant at FDR < 0.05; "
          f"tests used: {pw['test'].value_counts().to_dict()}")
    print(f"Multi-group: {rm['manifest']['n_samples']} wells; omnibus significant: "
          f"{int(om['significant'].sum())} of {len(om)}; tests used: {om['test'].value_counts().to_dict()}")
    print(check)
    print("Article Figure 4 = two_group/figures/violin_panel.*")
    print("Article Figure 5 = two_group/figures/volcano_Actin_disruptors.*")
    print("Article Figure 6 = multi_group/figures/pca.*")
    print("Article Figure 7 = multi_group/figures/shape_correlation_spearman.*")
    print("Article Table 5 and File S2 = two_group/tables/pairwise_vs_control.csv")
    print("Article File S3 = multi_group/tables/omnibus_across_groups.csv")
    print(f"All outputs written to: {os.path.abspath(a.outdir)}")


if __name__ == "__main__":
    main()
