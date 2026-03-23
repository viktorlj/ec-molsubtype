#!/usr/bin/env python3
"""Generate publication-grade manuscript figures for ec-molsubtype.

Produces:
  figures/manuscript_fig1.pdf  — Algorithm schematic + validation overview (2x3)
  figures/manuscript_fig2.pdf  — Multi-panel generalizability (1x2)
  figures/manuscript_fig3.pdf  — Secondary evidence and subtype separation (1x3)
  figures/manuscript_sfig1.pdf — TP53 DBD missense impact (supplementary)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
import seaborn as sns

# ---------- Publication style (identical to generate_figures.py) ----------
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

# Subtype colors — distinct, colorblind-safe (identical to generate_figures.py)
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


def panel_label(ax, label, x=-0.12, y=1.06):
    """Add a bold lowercase panel label."""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load GENIE multi-panel, TP53/CIN, and TCGA validation data."""
    df = pd.read_csv("data/genie/ec_multi_panel_classifications.tsv", sep="\t")
    df["panel"] = df["SEQ_ASSAY_ID"].map(PANEL_NAMES)
    tp53 = pd.read_csv("data/genie/ec_msk_tp53_cin.tsv", sep="\t")
    tcga = pd.read_csv("data/tcga/tcga_validation_results.tsv", sep="\t")
    return df, tp53, tcga


# ============================================================
# Helper: draw confusion matrix heatmap
# ============================================================
def _draw_confusion_matrix(ax, cm, subtypes, title, n_total, accuracy):
    """Draw a confusion matrix heatmap on the given axes."""
    cm_pct = cm.astype(float)
    for i in range(4):
        row_total = cm[i].sum()
        if row_total > 0:
            cm_pct[i] = 100 * cm[i] / row_total

    im = ax.imshow(cm_pct, cmap="Blues", vmin=0, vmax=100, aspect="equal")

    for i in range(4):
        for j in range(4):
            val = cm[i, j]
            pct = cm_pct[i, j]
            if val > 0:
                color = "white" if pct > 60 else "black"
                ax.text(j, i, f"{val}\n({pct:.0f}%)", ha="center", va="center",
                        fontsize=5.5, color=color,
                        fontweight="bold" if i == j else "normal")

    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(subtypes, fontsize=5.5, rotation=30, ha="right")
    ax.set_yticklabels(subtypes, fontsize=5.5)
    ax.set_xlabel("Predicted", fontsize=6)
    ax.set_ylabel("Ground truth", fontsize=6)
    ax.set_title(f"{title}\n(n={n_total}, acc={accuracy:.1f}%)", fontsize=7)


# ============================================================
# Helper: stacked horizontal bar for cross-panel comparison
# ============================================================
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
# Figure 1: Algorithm schematic + validation overview (2x3)
# ============================================================
def manuscript_fig1(df: pd.DataFrame, tcga: pd.DataFrame):
    fig = plt.figure(figsize=(7.0, 5.0))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.40)

    subtypes = SUBTYPE_ORDER

    # --- (a) Algorithm flowchart ---
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "a", x=-0.15, y=1.08)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Classification algorithm", fontsize=7, pad=8)

    # Flow positions: decision diamonds on left, outcome boxes on right
    steps = [
        {"y": 8.8, "test": "POLE EDM?", "subtype": "POLEmut", "color": SUBTYPE_COLORS["POLEmut"]},
        {"y": 6.6, "test": "MMR\ndeficient?", "subtype": "MMRd", "color": SUBTYPE_COLORS["MMRd"]},
        {"y": 4.4, "test": "TP53\npathogenic?", "subtype": "p53abn", "color": SUBTYPE_COLORS["p53abn"]},
        {"y": 2.2, "test": None, "subtype": "NSMP", "color": SUBTYPE_COLORS["NSMP"]},
    ]

    diamond_x = 2.8
    box_x = 7.5
    diamond_w = 2.4
    diamond_h = 1.5
    box_w = 2.8
    box_h = 0.9

    for i, step in enumerate(steps):
        cy = step["y"]

        if step["test"] is not None:
            # Diamond (decision)
            diamond = mpatches.FancyBboxPatch(
                (diamond_x - diamond_w / 2, cy - diamond_h / 2),
                diamond_w, diamond_h,
                boxstyle="round,pad=0.05",
                facecolor="#f0f0f0", edgecolor="#333333", linewidth=0.7,
                transform=ax.transData
            )
            # Actually draw a rotated square (diamond shape)
            dx = diamond_w / 2
            dy = diamond_h / 2
            diamond_pts = np.array([
                [diamond_x, cy + dy],
                [diamond_x + dx, cy],
                [diamond_x, cy - dy],
                [diamond_x - dx, cy],
            ])
            diamond_patch = mpatches.Polygon(
                diamond_pts, closed=True,
                facecolor="#f7f7f7", edgecolor="#333333", linewidth=0.7
            )
            ax.add_patch(diamond_patch)
            ax.text(diamond_x, cy, step["test"], ha="center", va="center",
                    fontsize=5, fontweight="bold", color="#333")

            # "Yes" arrow to subtype box
            ax.annotate("", xy=(box_x - box_w / 2 - 0.1, cy),
                        xytext=(diamond_x + dx + 0.1, cy),
                        arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.7))
            ax.text((diamond_x + dx + box_x - box_w / 2) / 2, cy + 0.2,
                    "Yes", ha="center", fontsize=5, color="#009E73", fontweight="bold")

            # "No" arrow going down
            if i < len(steps) - 1:
                next_y = steps[i + 1]["y"]
                target_y = next_y + (diamond_h / 2 if steps[i + 1]["test"] is not None else box_h / 2)
                ax.annotate("", xy=(diamond_x, target_y + 0.1),
                            xytext=(diamond_x, cy - dy - 0.1),
                            arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.7))
                ax.text(diamond_x + 0.3, (cy - dy + target_y) / 2,
                        "No", ha="left", fontsize=5, color="#888")
        else:
            # Final box (NSMP — no decision, just the default)
            # Draw an arrow arriving from above pointing to the box directly
            pass

        # Subtype outcome box
        box = FancyBboxPatch(
            (box_x - box_w / 2, cy - box_h / 2),
            box_w, box_h,
            boxstyle="round,pad=0.15",
            facecolor=step["color"], edgecolor="white", linewidth=0.5,
            alpha=0.85
        )
        ax.add_patch(box)
        ax.text(box_x, cy, step["subtype"], ha="center", va="center",
                fontsize=6, fontweight="bold", color="white")

        # For NSMP, draw an arrow from the diamond_x position to the box
        if step["test"] is None:
            ax.annotate("", xy=(box_x - box_w / 2 - 0.1, cy),
                        xytext=(diamond_x + 0.1, cy),
                        arrowprops=dict(arrowstyle="-|>", color="#333", lw=0.7))

    # Step labels on left
    for i, step in enumerate(steps):
        if step["test"] is not None:
            ax.text(0.3, step["y"], f"Step {i + 1}", ha="left", va="center",
                    fontsize=5, color="#666", style="italic")
        else:
            ax.text(0.3, step["y"], "Default", ha="left", va="center",
                    fontsize=5, color="#666", style="italic")

    # --- (b) TCGA confusion matrix ---
    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "b")

    cm_tcga = np.zeros((4, 4), dtype=int)
    for i, gt in enumerate(subtypes):
        for j, pred in enumerate(subtypes):
            cm_tcga[i, j] = ((tcga["gt_subtype"] == gt) & (tcga["pred_subtype"] == pred)).sum()

    accuracy_tcga = 100 * tcga["match"].mean()
    _draw_confusion_matrix(ax, cm_tcga, subtypes, "TCGA validation",
                           len(tcga), accuracy_tcga)

    # --- (c) CPTAC confusion matrix (hardcoded) ---
    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "c")

    cm_cptac = np.array([
        [6,  0,  0, 1],
        [0, 25,  0, 0],
        [0,  0, 14, 6],
        [0,  0,  2, 41],
    ], dtype=int)
    n_cptac = cm_cptac.sum()
    acc_cptac = 100 * np.trace(cm_cptac) / n_cptac
    _draw_confusion_matrix(ax, cm_cptac, subtypes, "CPTAC validation",
                           n_cptac, acc_cptac)

    # --- (d) Per-subtype F1 grouped bar (TCGA + CPTAC) ---
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "d")

    def compute_f1_metrics(cm):
        """Compute precision, recall, F1 per subtype from confusion matrix."""
        metrics = []
        for i, st in enumerate(subtypes):
            tp = cm[i, i]
            gt_n = cm[i, :].sum()
            pred_n = cm[:, i].sum()
            recall = tp / gt_n if gt_n > 0 else 0
            precision = tp / pred_n if pred_n > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            metrics.append({"subtype": st, "F1": 100 * f1, "Recall": 100 * recall,
                            "Precision": 100 * precision})
        return pd.DataFrame(metrics)

    tcga_metrics = compute_f1_metrics(cm_tcga)
    cptac_metrics = compute_f1_metrics(cm_cptac)

    x = np.arange(4)
    w = 0.35
    bars_tcga = ax.bar(x - w / 2, tcga_metrics["F1"], w, label="TCGA",
                       color=[SUBTYPE_COLORS[s] for s in subtypes],
                       edgecolor="white", linewidth=0.3, alpha=0.9)
    bars_cptac = ax.bar(x + w / 2, cptac_metrics["F1"], w, label="CPTAC",
                        color=[SUBTYPE_COLORS[s] for s in subtypes],
                        edgecolor="white", linewidth=0.3, alpha=0.5,
                        hatch="//")

    for i in range(4):
        ax.text(x[i] - w / 2, tcga_metrics["F1"].iloc[i] + 1.5,
                f"{tcga_metrics['F1'].iloc[i]:.0f}", ha="center", fontsize=5, va="bottom")
        ax.text(x[i] + w / 2, cptac_metrics["F1"].iloc[i] + 1.5,
                f"{cptac_metrics['F1'].iloc[i]:.0f}", ha="center", fontsize=5, va="bottom")

    ax.set_xticks(x)
    ax.set_xticklabels(subtypes, fontsize=6)
    ax.set_ylabel("F1 score (%)")
    ax.set_ylim(0, 115)
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.set_title("Per-subtype F1 scores", fontsize=7)

    # --- (e) TMB by subtype violin (from GENIE data) ---
    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "e")

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
    ax.set_ylabel("TMB (mut/Mb)")
    ax.set_xlabel("")
    ax.set_title("Tumor mutational burden", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.axhline(10, color="#888", ls=":", lw=0.5, zorder=0)
    ax.axhline(100, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(3.6, 10, "10", fontsize=5, color="#888", va="bottom")
    ax.text(3.6, 100, "100", fontsize=5, color="#888", va="bottom")

    # --- (f) FGA by subtype violin (from GENIE data) ---
    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "f")

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
    ax.set_title("Copy number burden", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_ylim(-0.02, 1.0)
    ax.axhline(0.2, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(3.6, 0.2, "0.2", fontsize=5, color="#888", va="bottom")

    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"manuscript_fig1.{fmt}")
    plt.close(fig)
    print(f"Figure 1 saved: {OUTDIR / 'manuscript_fig1.pdf'}")


# ============================================================
# Figure 2: Multi-panel generalizability (1x2)
# ============================================================
def manuscript_fig2(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # --- (a) All EC histologies ---
    ax = axes[0]
    panel_label(ax, "a")
    _stacked_bar(df, ax, "All endometrial cancer")

    # --- (b) UEC only ---
    ax = axes[1]
    panel_label(ax, "b")
    uec = df[df["ONCOTREE_CODE"] == "UEC"]
    _stacked_bar(uec, ax, "Endometrioid (UEC) only")

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"manuscript_fig2.{fmt}")
    plt.close(fig)
    print(f"Figure 2 saved: {OUTDIR / 'manuscript_fig2.pdf'}")


# ============================================================
# Figure 3: Secondary evidence and subtype separation (1x3)
# ============================================================
def manuscript_fig3(df: pd.DataFrame, tp53: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5))

    # --- (a) TMB vs indel fraction scatter with zone annotations ---
    ax = axes[0]
    panel_label(ax, "a")

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

    # Quadrant lines and zone annotations
    ax.axvline(100, color="#888", ls=":", lw=0.5, zorder=0)
    ax.axhline(0.15, color="#888", ls=":", lw=0.5, zorder=0)
    ax.text(120, 0.02, "POLEmut\nzone", fontsize=5, color=SUBTYPE_COLORS["POLEmut"],
            ha="left", fontweight="bold")
    ax.text(3, 0.4, "MMRd\nzone", fontsize=5, color=SUBTYPE_COLORS["MMRd"],
            ha="center", fontweight="bold")

    # --- (b) TP53 status vs FGA violin ---
    ax = axes[1]
    panel_label(ax, "b")

    tp53_colors = {"wild_type": "#56B4E9", "mutated": "#D55E00"}

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

    # --- (c) MSI proxy (pseudo_msi_pct) by subtype strip/jitter ---
    ax = axes[2]
    panel_label(ax, "c")

    plot_msi = df[df["pseudo_msi_pct"].notna()].copy()
    plot_msi["subtype"] = pd.Categorical(plot_msi["subtype"],
                                         categories=SUBTYPE_ORDER, ordered=True)

    sns.stripplot(data=plot_msi, x="subtype", y="pseudo_msi_pct", hue="subtype",
                  order=SUBTYPE_ORDER, palette=SUBTYPE_COLORS, legend=False,
                  size=1.5, alpha=0.25, jitter=0.3, ax=ax, rasterized=True)

    # Add boxplot overlay for summary
    sns.boxplot(data=plot_msi, x="subtype", y="pseudo_msi_pct",
                order=SUBTYPE_ORDER,
                color="white", width=0.3, linewidth=0.5,
                fliersize=0, ax=ax,
                boxprops=dict(zorder=2, alpha=0.7),
                medianprops=dict(color="black", linewidth=1),
                whiskerprops=dict(linewidth=0.5),
                capprops=dict(linewidth=0.5))

    ax.set_ylabel("Pseudo MSI %")
    ax.set_xlabel("")
    ax.set_title("MSI proxy score by subtype")
    ax.axhline(9.0, color="#888", ls="--", lw=0.7, zorder=0)
    ax.text(3.6, 10, "threshold (9%)", fontsize=5, color="#888", va="bottom")
    ax.set_ylim(-2, 65)

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"manuscript_fig3.{fmt}")
    plt.close(fig)
    print(f"Figure 3 saved: {OUTDIR / 'manuscript_fig3.pdf'}")


# ============================================================
# Supplementary Figure 1: TP53 DBD missense impact
# ============================================================
def manuscript_sfig1(tcga: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(5.5, 2.8))

    subtypes = SUBTYPE_ORDER

    # --- (a) Before vs after comparison ---
    ax = axes[0]
    panel_label(ax, "a")

    # Before values (from validation run before DBD fix)
    before = {"POLEmut": 93.9, "MMRd": 83.8, "p53abn": 41.7, "NSMP": 95.9}
    # After values (current)
    after = {}
    for st in subtypes:
        gt_n = (tcga["gt_subtype"] == st).sum()
        tp = ((tcga["gt_subtype"] == st) & (tcga["pred_subtype"] == st)).sum()
        after[st] = 100 * tp / gt_n if gt_n > 0 else 0

    x = np.arange(4)
    w = 0.35
    ax.bar(x - w / 2, [before[s] for s in subtypes], w, label="Before (hotspots only)",
           color="#bbb", edgecolor="white", linewidth=0.3)
    bars_after = ax.bar(x + w / 2, [after[s] for s in subtypes], w,
                        label="After (all DBD missense)",
                        color=[SUBTYPE_COLORS[s] for s in subtypes],
                        edgecolor="white", linewidth=0.3)

    for i, st in enumerate(subtypes):
        delta = after[st] - before[st]
        if abs(delta) > 1:
            ax.annotate(f"+{delta:.0f}%", xy=(i + w / 2, after[st] + 1),
                        fontsize=6, ha="center", va="bottom",
                        color="#D55E00", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(subtypes, fontsize=6)
    ax.set_ylabel("Recall (%)")
    ax.set_ylim(0, 110)
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    ax.set_title("Impact of TP53 DBD missense\nreclassification on recall", fontsize=7)

    # --- (b) Misclassification waterfall ---
    ax = axes[1]
    panel_label(ax, "b")

    errors = (tcga[~tcga["match"]]
              .groupby(["gt_subtype", "pred_subtype"])
              .size()
              .reset_index(name="count"))
    errors = errors.sort_values("count", ascending=True)
    errors["label"] = errors["gt_subtype"] + " \u2192 " + errors["pred_subtype"]

    y = np.arange(len(errors))
    colors = [SUBTYPE_COLORS.get(row["gt_subtype"], "#999") for _, row in errors.iterrows()]
    ax.barh(y, errors["count"], color=colors, edgecolor="white",
            linewidth=0.3, height=0.6)

    for i, (_, row) in enumerate(errors.iterrows()):
        ax.text(row["count"] + 0.5, i, str(row["count"]), va="center", fontsize=6)

    ax.set_yticks(y)
    ax.set_yticklabels(errors["label"], fontsize=6)
    ax.set_xlabel("Count")
    n_errors = len(tcga[~tcga["match"]])
    ax.set_title(f"Remaining misclassifications\n({n_errors}/{len(tcga)} samples)", fontsize=7)

    plt.tight_layout()
    for fmt in ["pdf", "png"]:
        fig.savefig(OUTDIR / f"manuscript_sfig1.{fmt}")
    plt.close(fig)
    print(f"Supplementary Figure 1 saved: {OUTDIR / 'manuscript_sfig1.pdf'}")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    df, tp53, tcga = load_data()
    manuscript_fig1(df, tcga)
    manuscript_fig2(df)
    manuscript_fig3(df, tp53)
    manuscript_sfig1(tcga)
    print(f"\nAll manuscript figures saved to {OUTDIR}/")
