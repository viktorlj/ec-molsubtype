#!/usr/bin/env python3
"""Generate publication-grade figures for ec-molsubtype validation.

Produces:
  figures/fig1_validation_overview.pdf  — 4-panel: subtypes, TMB, FGA, I-index
  figures/fig2_tp53_cin.pdf             — TP53 vs CIN by histology
  figures/fig3_msi_proxy.pdf            — MSI proxy calibration
  figures/fig4_cross_panel.pdf          — Cross-panel concordance
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------- Publication style ----------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.0,
    "pdf.fonttype": 42,  # TrueType fonts in PDF
    "ps.fonttype": 42,
})

# Subtype colors — distinct, colorblind-safe
SUBTYPE_COLORS = {
    "POLEmut": "#009E73",  # bluish green
    "MMRd":    "#0072B2",  # blue
    "p53abn":  "#D55E00",  # vermillion
    "NSMP":    "#999999",  # gray
}
SUBTYPE_ORDER = ["POLEmut", "MMRd", "p53abn", "NSMP"]

# Panel display names (shortened)
PANEL_NAMES = {
    "MSK-IMPACT468": "MSK-468",
    "MSK-IMPACT505": "MSK-505",
    "UCSF-IDTV5-TO": "UCSF",
    "DFCI-ONCOPANEL-3.1": "DFCI-3.1",
    "DFCI-ONCOPANEL-3": "DFCI-3",
}
PANEL_ORDER = ["MSK-468", "MSK-505", "UCSF", "DFCI-3.1", "DFCI-3"]

TCGA_REF = {"POLEmut": 7, "MMRd": 28, "p53abn": 26, "NSMP": 39}

OUTDIR = Path("figures")
OUTDIR.mkdir(exist_ok=True)


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv("data/genie/ec_multi_panel_classifications.tsv", sep="\t")
    df["panel"] = df["SEQ_ASSAY_ID"].map(PANEL_NAMES)
    tp53 = pd.read_csv("data/genie/ec_msk_tp53_cin.tsv", sep="\t")
    return df, tp53


def panel_label(ax, label, x=-0.12, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


# ============================================================
# Figure 1: Validation overview (4-panel)
# ============================================================
def fig1_validation_overview(df: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.5))

    # --- A: Subtype distribution by panel ---
    ax = axes[0, 0]
    panel_label(ax, "a")

    # Compute percentages per panel
    panels = []
    for panel in PANEL_ORDER:
        sub = df[df["panel"] == panel]
        n = len(sub)
        if n == 0:
            continue
        row = {"panel": f"{panel}\n(n={n})"}
        for st in SUBTYPE_ORDER:
            row[st] = 100 * (sub["subtype"] == st).sum() / n
        panels.append(row)

    # Add TCGA reference
    panels.append({"panel": "TCGA\nref.", **{st: TCGA_REF[st] for st in SUBTYPE_ORDER}})

    pdf = pd.DataFrame(panels)
    x = np.arange(len(pdf))
    width = 0.18
    for i, st in enumerate(SUBTYPE_ORDER):
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, pdf[st], width, label=st,
                      color=SUBTYPE_COLORS[st], edgecolor="white", linewidth=0.3)
        # Add percentage labels on TCGA bars
        if True:
            for j, (xi, val) in enumerate(zip(x + offset, pdf[st])):
                if j == len(pdf) - 1 and val > 3:  # TCGA row
                    ax.text(xi, val + 1, f"{val:.0f}", ha="center", va="bottom",
                            fontsize=5, color=SUBTYPE_COLORS[st])

    ax.set_xticks(x)
    ax.set_xticklabels(pdf["panel"], fontsize=6)
    ax.set_ylabel("Proportion (%)")
    ax.set_ylim(0, 75)
    ax.legend(loc="upper right", frameon=False, ncol=2, fontsize=6)
    ax.set_title("Subtype distribution across panels")

    # --- B: TMB by subtype ---
    ax = axes[0, 1]
    panel_label(ax, "b")

    plot_df = df[df["tmb"].notna() & (df["tmb"] > 0)].copy()
    plot_df["subtype"] = pd.Categorical(plot_df["subtype"], categories=SUBTYPE_ORDER, ordered=True)

    sns.violinplot(data=plot_df, x="subtype", y="tmb", hue="subtype",
                   order=SUBTYPE_ORDER, palette=SUBTYPE_COLORS, legend=False,
                   inner=None, linewidth=0.5, cut=0, density_norm="width", ax=ax)
    sns.boxplot(data=plot_df, x="subtype", y="tmb", order=SUBTYPE_ORDER,
                color="white", width=0.15, linewidth=0.5,
                fliersize=0, ax=ax,
                boxprops=dict(zorder=2), medianprops=dict(color="black", linewidth=1))
    ax.set_yscale("symlog", linthresh=1)
    ax.set_yticks([0, 1, 5, 10, 50, 100, 500])
    ax.set_yticklabels(["0", "1", "5", "10", "50", "100", "500"])
    ax.set_ylabel("TMB (mutations/Mb)")
    ax.set_xlabel("")
    ax.set_title("Tumor mutational burden")

    # Reference lines
    ax.axhline(10, color="#888", ls=":", lw=0.5, zorder=0)
    ax.axhline(100, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(3.6, 10, "10", fontsize=5, color="#888", va="bottom")
    ax.text(3.6, 100, "100", fontsize=5, color="#888", va="bottom")

    # --- C: FGA by subtype ---
    ax = axes[1, 0]
    panel_label(ax, "c")

    plot_fga = df[df["fga"].notna()].copy()
    plot_fga["subtype"] = pd.Categorical(plot_fga["subtype"], categories=SUBTYPE_ORDER, ordered=True)

    sns.violinplot(data=plot_fga, x="subtype", y="fga", hue="subtype",
                   order=SUBTYPE_ORDER, palette=SUBTYPE_COLORS, legend=False,
                   inner=None, linewidth=0.5, cut=0, density_norm="width", ax=ax)
    sns.boxplot(data=plot_fga, x="subtype", y="fga", order=SUBTYPE_ORDER,
                color="white", width=0.15, linewidth=0.5,
                fliersize=0, ax=ax,
                boxprops=dict(zorder=2), medianprops=dict(color="black", linewidth=1))

    ax.set_ylabel("Fraction genome altered")
    ax.set_xlabel("")
    ax.set_title("Copy number burden")
    ax.set_ylim(-0.02, 1.0)
    ax.axhline(0.2, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(3.6, 0.2, "0.2", fontsize=5, color="#888", va="bottom")

    # --- D: Indel fraction by subtype ---
    ax = axes[1, 1]
    panel_label(ax, "d")

    plot_ii = df[df["i_index"].notna()].copy()
    plot_ii["subtype"] = pd.Categorical(plot_ii["subtype"], categories=SUBTYPE_ORDER, ordered=True)

    sns.violinplot(data=plot_ii, x="subtype", y="i_index", hue="subtype",
                   order=SUBTYPE_ORDER, palette=SUBTYPE_COLORS, legend=False,
                   inner=None, linewidth=0.5, cut=0, density_norm="width", ax=ax)
    sns.boxplot(data=plot_ii, x="subtype", y="i_index", order=SUBTYPE_ORDER,
                color="white", width=0.15, linewidth=0.5,
                fliersize=0, ax=ax,
                boxprops=dict(zorder=2), medianprops=dict(color="black", linewidth=1))

    ax.set_ylabel("Indel fraction (I-index)")
    ax.set_xlabel("")
    ax.set_title("Substitution spectrum: indel fraction")
    ax.set_ylim(-0.02, 0.8)
    ax.axhline(0.15, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(3.6, 0.15, "0.15", fontsize=5, color="#888", va="bottom")

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"fig1_validation_overview.{fmt}")
    plt.close(fig)
    print(f"Figure 1 saved: {OUTDIR / 'fig1_validation_overview.pdf'}")


# ============================================================
# Figure 2: TP53 and chromosomal instability
# ============================================================
def fig2_tp53_cin(tp53: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5))

    tp53_colors = {"wild_type": "#56B4E9", "mutated": "#D55E00"}

    # --- A: FGA by TP53 status ---
    ax = axes[0]
    panel_label(ax, "a")

    sns.violinplot(data=tp53, x="tp53_status", y="fga", hue="tp53_status",
                   order=["wild_type", "mutated"],
                   palette=tp53_colors, legend=False, inner=None, linewidth=0.5,
                   cut=0, density_norm="width", ax=ax)
    sns.boxplot(data=tp53, x="tp53_status", y="fga",
                order=["wild_type", "mutated"],
                color="white", width=0.2, linewidth=0.5,
                fliersize=0, ax=ax,
                boxprops=dict(zorder=2), medianprops=dict(color="black", linewidth=1))

    ax.set_ylabel("Fraction genome altered")
    ax.set_xlabel("")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["TP53 WT", "TP53 mut"])
    ax.set_title("TP53 status vs CIN")

    # Add median annotations
    for i, status in enumerate(["wild_type", "mutated"]):
        med = tp53[tp53["tp53_status"] == status]["fga"].median()
        ax.text(i, med + 0.03, f"{med:.3f}", ha="center", fontsize=6, style="italic")

    # --- B: FGA by histological subtype ---
    ax = axes[1]
    panel_label(ax, "b")

    hist_order = ["UEC", "USC", "UCS", "UCEC", "UMEC"]
    sub = tp53[tp53["ONCOTREE_CODE"].isin(hist_order)]

    sns.boxplot(data=sub, x="ONCOTREE_CODE", y="fga", hue="ONCOTREE_CODE",
                order=hist_order,
                palette="Set2", legend=False, linewidth=0.5, fliersize=1, ax=ax,
                medianprops=dict(color="black", linewidth=1))
    ax.set_ylabel("Fraction genome altered")
    ax.set_xlabel("")
    ax.set_title("CIN by histology")
    ax.set_ylim(-0.02, 1.0)

    # Add n and TP53 rate as x-tick sub-labels
    n_labels = []
    for code in hist_order:
        subset = sub[sub["ONCOTREE_CODE"] == code]
        tp53_pct = 100 * (subset["tp53_status"] == "mutated").sum() / len(subset) if len(subset) > 0 else 0
        n_labels.append(f"{code}\n({tp53_pct:.0f}% TP53)")
    ax.set_xticks(range(len(hist_order)))
    ax.set_xticklabels(n_labels, fontsize=5.5)

    # --- C: FGA by TP53 variant class ---
    ax = axes[2]
    panel_label(ax, "c")

    tp53_mut = tp53[tp53["tp53_status"] == "mutated"].copy()
    class_order = ["hotspot", "truncating", "missense", "other"]
    tp53_mut = tp53_mut[tp53_mut["tp53_variant_class"].isin(class_order)]

    sns.boxplot(data=tp53_mut, x="tp53_variant_class", y="fga", hue="tp53_variant_class",
                order=class_order,
                palette="Oranges_d", legend=False, linewidth=0.5, fliersize=1, ax=ax,
                medianprops=dict(color="black", linewidth=1))
    ax.set_ylabel("Fraction genome altered")
    ax.set_xlabel("")
    ax.set_title("CIN by TP53 variant class")
    ax.set_ylim(-0.02, 1.0)

    # Add n per group
    for i, cls in enumerate(class_order):
        n = (tp53_mut["tp53_variant_class"] == cls).sum()
        ax.text(i, 0.95, f"n={n}", ha="center", fontsize=5, color="#666")

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"fig2_tp53_cin.{fmt}")
    plt.close(fig)
    print(f"Figure 2 saved: {OUTDIR / 'fig2_tp53_cin.pdf'}")


# ============================================================
# Figure 3: MSI proxy performance
# ============================================================
def fig3_msi_proxy(df: pd.DataFrame, tp53: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5))

    # --- A: Pseudo MSI % distribution by subtype ---
    ax = axes[0]
    panel_label(ax, "a")

    # Only samples with some MSI signal for visibility
    for st in SUBTYPE_ORDER:
        sub = df[df["subtype"] == st]
        vals = sub["pseudo_msi_pct"].dropna()
        # Histogram-like KDE
        if len(vals) > 10 and vals.std() > 0:
            sns.kdeplot(vals, ax=ax, label=st, color=SUBTYPE_COLORS[st],
                        linewidth=1.2, clip=(0, 100))

    ax.axvline(9.0, color="#888", ls="--", lw=0.7)
    ax.text(10, ax.get_ylim()[1] * 0.9, "threshold\n(9%)", fontsize=5, color="#888")
    ax.set_xlabel("Pseudo MSI %")
    ax.set_ylabel("Density")
    ax.set_title("MSI proxy score by subtype")
    ax.legend(frameon=False, fontsize=6)
    ax.set_xlim(-2, 60)

    # --- B: cMS genes hit vs I-index (2D scatter) ---
    ax = axes[1]
    panel_label(ax, "b")

    for st in SUBTYPE_ORDER:
        sub = df[df["subtype"] == st]
        ax.scatter(sub["cms_high_spec"], sub["i_index"],
                   c=SUBTYPE_COLORS[st], label=st, s=3, alpha=0.3,
                   edgecolors="none", rasterized=True)

    ax.set_xlabel("High-specificity cMS genes hit")
    ax.set_ylabel("Indel fraction (I-index)")
    ax.set_title("MSI indicator genes vs indel fraction")
    ax.legend(frameon=False, fontsize=5, markerscale=3, loc="upper right")
    ax.set_xlim(-0.5, 7)
    ax.set_ylim(-0.02, 0.7)

    # --- C: TMB vs I-index colored by subtype ---
    ax = axes[2]
    panel_label(ax, "c")

    for st in SUBTYPE_ORDER:
        sub = df[(df["subtype"] == st) & (df["tmb"] > 0)]
        ax.scatter(sub["tmb"], sub["i_index"],
                   c=SUBTYPE_COLORS[st], label=st, s=3, alpha=0.3,
                   edgecolors="none", rasterized=True)

    ax.set_xscale("log")
    ax.set_xlabel("TMB (mutations/Mb)")
    ax.set_ylabel("Indel fraction (I-index)")
    ax.set_title("TMB vs indel fraction")
    ax.legend(frameon=False, fontsize=5, markerscale=3, loc="upper left")

    # Quadrant lines
    ax.axvline(100, color="#888", ls=":", lw=0.5, zorder=0)
    ax.axhline(0.15, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(120, 0.02, "POLEmut\nzone", fontsize=5, color="#888", ha="left")
    ax.text(3, 0.4, "MMRd\nzone", fontsize=5, color="#888", ha="center")

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"fig3_msi_proxy.{fmt}")
    plt.close(fig)
    print(f"Figure 3 saved: {OUTDIR / 'fig3_msi_proxy.pdf'}")


# ============================================================
# Figure 4: Cross-panel and UEC-specific comparison
# ============================================================
def fig4_cross_panel(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # --- A: All histologies ---
    ax = axes[0]
    panel_label(ax, "a")

    _stacked_bar(df, ax, "All endometrial cancer")

    # --- B: UEC only ---
    ax = axes[1]
    panel_label(ax, "b")

    uec = df[df["ONCOTREE_CODE"] == "UEC"]
    _stacked_bar(uec, ax, "Endometrioid (UEC) only")

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"fig4_cross_panel.{fmt}")
    plt.close(fig)
    print(f"Figure 4 saved: {OUTDIR / 'fig4_cross_panel.pdf'}")


def _stacked_bar(df: pd.DataFrame, ax, title: str):
    """Horizontal stacked bar chart of subtype proportions."""
    rows = []
    for panel in PANEL_ORDER:
        sub = df[df["panel"] == panel]
        n = len(sub)
        if n < 20:
            continue
        row = {"panel": f"{panel} (n={n})"}
        for st in SUBTYPE_ORDER:
            row[st] = 100 * (sub["subtype"] == st).sum() / n
        rows.append(row)
    rows.append({"panel": "TCGA ref.", **{st: TCGA_REF[st] for st in SUBTYPE_ORDER}})

    bdf = pd.DataFrame(rows)
    y = np.arange(len(bdf))

    left = np.zeros(len(bdf))
    for st in SUBTYPE_ORDER:
        ax.barh(y, bdf[st], left=left, color=SUBTYPE_COLORS[st],
                label=st, edgecolor="white", linewidth=0.3, height=0.6)
        # Add percentage text
        for i, (l, w) in enumerate(zip(left, bdf[st])):
            if w > 5:
                ax.text(l + w / 2, i, f"{w:.0f}", ha="center", va="center",
                        fontsize=5, color="white", fontweight="bold")
        left += bdf[st].values

    ax.set_yticks(y)
    ax.set_yticklabels(bdf["panel"], fontsize=6)
    ax.set_xlabel("Proportion (%)")
    ax.set_xlim(0, 100)
    ax.set_title(title)
    ax.legend(loc="lower right", frameon=False, fontsize=6, ncol=2)
    ax.invert_yaxis()


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    df, tp53 = load_data()
    fig1_validation_overview(df)
    fig2_tp53_cin(tp53)
    fig3_msi_proxy(df, tp53)
    fig4_cross_panel(df)
    print(f"\nAll figures saved to {OUTDIR}/")
