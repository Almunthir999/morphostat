"""
morphostat.viz
==============
Publication-ready, *editable* vector figures (SVG + PDF).

All text is kept as real text in the SVG (``svg.fonttype='none'``) so figures
remain fully editable in Illustrator / Inkscape / Affinity. Every helper saves
both a .svg and a .pdf and returns the matplotlib Figure.
"""
from __future__ import annotations

import os
from typing import Sequence

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
from scipy import stats as _sps
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform


def set_style():
    sns.set_theme(context="paper", style="ticks")
    plt.rcParams.update({
        "svg.fonttype": "none",          # keep text editable in vector output
        "pdf.fonttype": 42,              # TrueType, editable in PDF
        "font.family": "sans-serif",
        # Inter Medium where available (slightly heavier, very legible), else Arial / DejaVu
        "font.sans-serif": ["Inter", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.weight": 500, "axes.labelweight": 500,
        "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.labelsize": 12, "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
        "legend.fontsize": 10.5, "legend.title_fontsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 1.2, "xtick.major.width": 1.1, "ytick.major.width": 1.1,
        "xtick.major.size": 4, "ytick.major.size": 4,
        "text.color": "#111111", "axes.labelcolor": "#111111",
        "xtick.color": "#111111", "ytick.color": "#111111",
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


# ---------------------------------------------------------------------------
# Readable feature names (full names are always kept in the output tables)
_SHAPE_NAMES = {
    "Area": "Area", "Perimeter": "Perimeter", "FormFactor": "Circularity (form factor)",
    "Solidity": "Solidity", "Compactness": "Compactness", "Eccentricity": "Eccentricity",
    "MajorAxisLength": "Major-axis length", "MinorAxisLength": "Minor-axis length",
    "Orientation": "Orientation angle", "Extent": "Extent",
    "MaxFeretDiameter": "Max Feret diameter", "MinFeretDiameter": "Min Feret diameter",
    "MaximumRadius": "Maximum radius", "MeanRadius": "Mean radius", "MedianRadius": "Median radius",
}
_COMPARTMENTS = {"Cells": "Cell", "Nuclei": "Nuc", "Cytoplasm": "Cyto"}


def pretty_name(feature: str, max_len: int = 42) -> str:
    """Shorten CellProfiler-style names for plotting, e.g.
    Cells_AreaShape_FormFactor -> 'Circularity (form factor)',
    Nuclei_Intensity_StdIntensity_CorrTub -> 'Nuc_Intensity_StdIntensity [Tub]'.
    Names that do not follow the CellProfiler pattern are returned unchanged."""
    f = str(feature)
    parts = f.split("_")
    if len(parts) >= 3 and parts[0] in _COMPARTMENTS and parts[1] == "AreaShape" and len(parts) == 3:
        nice = _SHAPE_NAMES.get(parts[2])
        if nice:
            return nice if parts[0] == "Cells" else f"{_COMPARTMENTS[parts[0]]} {nice.lower()}"
    if parts and parts[0] in _COMPARTMENTS:
        parts[0] = _COMPARTMENTS[parts[0]]
        if len(parts) >= 3 and parts[1] == "Intensity":
            parts[2] = parts[2].replace("Intensity", "") or parts[2]
        ci = next((i for i, t in enumerate(parts) if t.startswith("Corr") and len(t) > 4 and i > 1), None)
        if ci is not None:
            chan = parts.pop(ci)[4:]
            f = "_".join(parts) + f" [{chan}]"
        else:
            f = "_".join(parts)
    return f if len(f) <= max_len else f[: max_len - 1] + "\u2026"


def _save(fig, outdir: str, name: str) -> dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    paths = {}
    for ext in ("svg", "pdf", "png"):   # svg+pdf are the editable deliverables; png is a preview
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p)
        paths[ext] = p
    return paths


def _palette(groups: Sequence[str], control: str | None = None) -> dict:
    groups = list(dict.fromkeys(groups))
    base = sns.color_palette("colorblind", n_colors=max(len(groups), 3))
    pal = {g: base[i % len(base)] for i, g in enumerate(groups)}
    if control in pal:
        pal[control] = (0.45, 0.45, 0.45)  # grey control
    return pal


# ---------------------------------------------------------------------------
def violin_panel(
    df: pd.DataFrame, features: Sequence[str], group_col: str,
    outdir: str, name: str = "violin_panel", control: str | None = None,
    order: Sequence[str] | None = None, ncols: int = 3, title: str | None = None,
    stats: pd.DataFrame | None = None, ylabel: str = "Robust z-score (vs control)",
) -> dict[str, str]:
    """Grid of violin + box + jittered points for selected features by group.
    If ``stats`` (the pairwise table) is given, the adjusted p-value of each
    treatment-vs-control comparison is printed in the panel."""
    set_style()
    features = [f for f in features if f in df.columns]
    n = len(features); ncols = min(ncols, n); nrows = int(np.ceil(n / ncols))
    if order is None:
        order = list(df[group_col].dropna().unique())
        if control in order:
            order = [control] + [g for g in order if g != control]
    pal = _palette(order, control)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.75 * ncols, 2.65 * nrows), squeeze=False)
    for i, feat in enumerate(features):
        ax = axes[i // ncols][i % ncols]
        sns.violinplot(data=df, x=group_col, y=feat, order=order, hue=group_col,
                       palette=pal, inner=None, cut=0, linewidth=0.9, legend=False, ax=ax)
        for coll in ax.collections:
            coll.set_alpha(0.6)
        sns.boxplot(data=df, x=group_col, y=feat, order=order, width=0.16,
                    showcaps=True, boxprops={"facecolor": "white", "zorder": 3},
                    whiskerprops={"linewidth": 0.9}, fliersize=0, linewidth=0.9, ax=ax)
        sns.stripplot(data=df, x=group_col, y=feat, order=order, color="black",
                      size=2.0, alpha=0.4, jitter=0.18, ax=ax)
        ax.set_title(pretty_name(feat, 30), fontsize=11.5, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(ylabel if i % ncols == 0 else "", fontsize=10.5)
        rot = 0 if len(order) <= 2 else 30
        ax.tick_params(axis="x", rotation=rot, labelsize=10)
        if rot:
            for lbl in ax.get_xticklabels():
                lbl.set_ha("right")
        if stats is not None and len(order) == 2:
            row = stats[(stats["feature"] == feat) & (stats["group"] == order[1])]
            if len(row) and pd.notna(row["fdr_bh"].iloc[0]):
                q = float(row["fdr_bh"].iloc[0])
                y0, y1 = ax.get_ylim()
                ax.set_ylim(y0, y1 + 0.20 * (y1 - y0))      # headroom for the p-value
                ax.text(0.98, 0.98, "adj. p = " + _fmt_p(q), transform=ax.transAxes,
                        ha="right", va="top", fontsize=10, color="#111111")
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    if title:
        fig.suptitle(title, y=1.01, fontsize=13.5, fontweight="bold")
    fig.tight_layout()
    paths = _save(fig, outdir, name)
    plt.close(fig)
    return paths


def _fmt_p(q: float) -> str:
    if q >= 0.001:
        return f"{q:.3f}"
    m, e = f"{q:.1e}".split("e")
    return f"{m}\u00d710$^{{{int(e)}}}$"


# ---------------------------------------------------------------------------
def pca_plot(
    df: pd.DataFrame, feature_cols: Sequence[str], group_col: str,
    outdir: str, name: str = "pca", control: str | None = None,
    order: Sequence[str] | None = None, draw_ellipses: bool = True,
    title: str | None = None,
) -> tuple[dict[str, str], pd.DataFrame]:
    """PCA scatter of the full morphological signature, coloured by group."""
    set_style()
    X = df[list(feature_cols)].to_numpy(float)
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(X)
    ev = pca.explained_variance_ratio_ * 100
    scores = df[[group_col]].copy()
    scores["PC1"], scores["PC2"] = pcs[:, 0], pcs[:, 1]

    if order is None:
        order = list(df[group_col].dropna().unique())
        if control in order:
            order = [control] + [g for g in order if g != control]
    pal = _palette(order, control)

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    ax.axhline(0, color="#dddddd", lw=0.8, zorder=0); ax.axvline(0, color="#dddddd", lw=0.8, zorder=0)
    for g in order:
        sub = scores[scores[group_col] == g]
        ax.scatter(sub["PC1"], sub["PC2"], s=34, alpha=0.85, label=str(g),
                   color=pal[g], edgecolor="white", linewidth=0.5, zorder=3)
        if draw_ellipses and len(sub) >= 3:
            _confidence_ellipse(sub["PC1"].to_numpy(), sub["PC2"].to_numpy(), ax,
                                n_std=2.0, edgecolor=pal[g], facecolor="none",
                                linewidth=0.9, alpha=0.55, zorder=2)
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}% variance)")
    ax.legend(title=group_col, frameon=False, loc="best", markerscale=1.1,
              handletextpad=0.3, borderaxespad=0.4)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    paths = _save(fig, outdir, name)
    plt.close(fig)
    return paths, scores


def _confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    if x.size < 3:
        return
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * n_std * np.sqrt(np.maximum(vals, 0))
    e = Ellipse((np.mean(x), np.mean(y)), width=w, height=h, angle=theta, **kwargs)
    ax.add_patch(e)


# ---------------------------------------------------------------------------
def correlation_heatmap(
    df: pd.DataFrame, feature_cols: Sequence[str], outdir: str,
    name: str = "correlation_heatmap", method: str = "pearson",
    cluster: bool = True, max_features: int = 40, title: str | None = None,
) -> dict[str, str]:
    """Clustered feature-feature correlation heatmap with readable feature names."""
    set_style()
    feats = list(feature_cols)
    if len(feats) > max_features:  # keep the most variable features for readability
        var = df[feats].var().sort_values(ascending=False)
        feats = var.index[:max_features].tolist()
    corr = df[feats].corr(method=method)
    if cluster and len(feats) > 2:
        d = (1 - corr.abs()).to_numpy(dtype=float).copy()
        d = (d + d.T) / 2.0           # enforce exact symmetry
        np.fill_diagonal(d, 0.0)
        Z = linkage(squareform(d, checks=False), method="average")
        idx = leaves_list(Z)
        corr = corr.iloc[idx, idx]
    labels = [pretty_name(f, 40) for f in corr.columns]
    corr.index, corr.columns = labels, labels
    n = len(feats)
    small = n <= 15
    side = 5.6 if small else max(7.0, 0.24 * n)
    fig, ax = plt.subplots(figsize=(side + 1.6, side))
    sns.heatmap(corr, cmap="RdBu_r", center=0, vmin=-1, vmax=1, square=True,
                linewidths=0.4 if small else 0, linecolor="white",
                cbar_kws={"shrink": 0.7, "label": f"{method.capitalize()} correlation"},
                xticklabels=True, yticklabels=True, ax=ax)
    ax.tick_params(labelsize=11 if small else 7.5, length=2)
    plt.setp(ax.get_xticklabels(), rotation=90)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10.5)
    cbar.set_label(f"{method.capitalize()} correlation", fontsize=12)
    ax.set_title(title or f"Feature correlation structure ({method.capitalize()})", pad=10)
    fig.tight_layout()
    paths = _save(fig, outdir, name)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------
def effect_volcano(
    pairwise: pd.DataFrame, group: str, outdir: str,
    name: str | None = None, alpha: float = 0.05, top_n_labels: int = 12,
    control: str | None = None,
) -> dict[str, str]:
    """Volcano-style summary: effect size vs -log10(adjusted p) for one treatment group.
    The strongest significant features are labelled outside the plot (higher in
    treatment on the right, lower on the left) with leader lines."""
    set_style()
    sub = pairwise[pairwise["group"] == group].copy()
    sub = sub.dropna(subset=["effect", "fdr_bh"])
    sub["nlog10fdr"] = -np.log10(sub["fdr_bh"].clip(lower=1e-300))
    UP, DOWN, NS = "#C8394B", "#3C64C0", "#C9CBD3"
    fig = plt.figure(figsize=(9.4, 5.6))
    ax = fig.add_axes([0.335, 0.17, 0.33, 0.72])
    sig = sub["fdr_bh"] < alpha
    up = sig & (sub["effect"] > 0); dn = sig & (sub["effect"] < 0)
    ax.scatter(sub.loc[~sig, "effect"], sub.loc[~sig, "nlog10fdr"], s=26, color=NS,
               edgecolor="white", linewidth=0.3, label="not significant", zorder=2)
    ax.scatter(sub.loc[up, "effect"], sub.loc[up, "nlog10fdr"], s=30, color=UP,
               edgecolor="white", linewidth=0.3, label="higher in treatment", zorder=3)
    ax.scatter(sub.loc[dn, "effect"], sub.loc[dn, "nlog10fdr"], s=30, color=DOWN,
               edgecolor="white", linewidth=0.3, label="lower in treatment", zorder=3)
    ax.axhline(-np.log10(alpha), color="#777777", ls="--", lw=0.9, zorder=1)
    ax.axvline(0, color="#dddddd", lw=0.9, zorder=0)
    et = sub["effect_type"].dropna().iloc[0] if sub["effect_type"].notna().any() else "effect"
    et_nice = {"cliffs_delta": "Cliff's \u03b4", "hedges_g": "Hedges' g"}.get(et, et)
    ref = control if control else "control"
    ax.set_xlabel(f"Effect size ({et_nice}): {group} relative to {ref}")
    ax.set_ylabel("\u2212log$_{10}$ (adjusted p)")
    ax.set_title(f"Differential morphology: {group} vs {ref}")
    lim = max(1.05, float(sub["effect"].abs().max()) * 1.08)
    ax.set_xlim(-lim, lim)
    ymax = float(sub["nlog10fdr"].max())
    ax.set_ylim(-0.3, ymax * 1.06)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
              handletextpad=0.2, columnspacing=1.2)

    # ---- labels outside the axes, one column per side, with leader lines ----
    sig_df = sub.loc[sig].drop_duplicates("feature")
    sig_df = sig_df.reindex(sig_df["effect"].abs().sort_values(ascending=False).index)
    k_side = max(1, top_n_labels // 2)
    for side, x_txt, ha, col in [(1, 1.04, "left", UP), (-1, -0.20, "right", DOWN)]:
        hits = sig_df[np.sign(sig_df["effect"]) == side].head(k_side)
        hits = hits.sort_values("nlog10fdr", ascending=False)
        for k, (_, r) in enumerate(hits.iterrows()):
            y_txt = 0.97 - k * (0.80 / max(k_side - 1, 1))
            ax.annotate(pretty_name(r["feature"], 48), xy=(r["effect"], r["nlog10fdr"]),
                        xycoords="data", xytext=(x_txt, y_txt), textcoords="axes fraction",
                        ha=ha, va="center", fontsize=10.5, color="#111111",
                        arrowprops=dict(arrowstyle="-", lw=0.8, color=col, alpha=0.6,
                                        shrinkA=1, shrinkB=2),
                        annotation_clip=False)
    ax.yaxis.set_label_coords(-0.11, 0.5)
    ax.yaxis.label.set_bbox(dict(facecolor="white", edgecolor="none", pad=1.5))
    ax.yaxis.label.set_zorder(5)
    paths = _save(fig, outdir, name or f"volcano_{_safe(group)}")
    plt.close(fig)
    return paths


def _short(s: str, n: int = 26) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "\u2026"


def _safe(s: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(s))
