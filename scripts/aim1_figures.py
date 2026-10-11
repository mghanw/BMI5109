"""Figures for Phase 2 Aim 1 (reports/figures/phase2_aim1_*.{png,svg}).

Usage:
    python scripts/aim1_figures.py [--run-name phase2_aim1]

Reads the aggregate tables in results/summary/ and, for the ECDF panel only, the cell-level scores
kept in outputs/<run_name>/scores.csv.gz (never committed). Colours: reference categorical slots
1-2 (DMSO blue, Vem orange) with circle / square markers as a secondary encoding.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results" / "summary"
FIGURES = ROOT / "reports" / "figures"
FUFA = "Fufa2019_AZD_up_SOX10_bound_pigmentation"
COLORS = {"DMSO": "#2a78d6", "Vem": "#eb6834"}
MARKERS = {"DMSO": "o", "Vem": "s"}
MUTED, INK, GRID = "#a9a8a2", "#0b0b0b", "#d9d8d4"


def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "legend.frameon": False, "axes.edgecolor": "#52514e",
        "xtick.color": "#52514e", "ytick.color": "#52514e", "text.color": INK,
    })


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(FIGURES / f"{name}.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)


def forest(mean, disp):
    gene = mean[mean.kind.eq("gene") & mean.readout.eq(FUFA)]
    lvr = disp[disp.kind.eq("gene") & disp.readout.eq(FUFA) & disp.statistic.eq("log_var_ratio")]
    order = gene[gene.condition.eq("DMSO")].sort_values("estimate").target.tolist()
    y = {t: i for i, t in enumerate(order)}
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 230 / 25.4), sharey=True)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.085, wspace=0.08)
    for ax, data, title, xlabel in [
            (axes[0], gene, "Mean shift", "Knockdown − NTC mean\n(DMSO NTC SD)"),
            (axes[1], lvr, "Dispersion change", "log(variance knockdown / variance NTC)")]:
        for condition, offset in (("DMSO", -0.18), ("Vem", 0.18)):
            d = data[data.condition.eq(condition)]
            ys = d.target.map(y).to_numpy() + offset
            sig = d.q_value.lt(0.05).to_numpy()
            ax.hlines(ys, d.ci_low, d.ci_high, color=COLORS[condition], linewidth=0.8, alpha=0.9)
            ax.scatter(d.estimate[sig], ys[sig], color=COLORS[condition], marker=MARKERS[condition], s=12,
                       zorder=3, label=f"{condition} (BH q < 0.05)")
            ax.scatter(d.estimate[~sig], ys[~sig], facecolor="white", edgecolor=COLORS[condition],
                       marker=MARKERS[condition], s=12, linewidth=0.8, zorder=3, label=f"{condition} (q ≥ 0.05)")
        ax.axvline(0, color=MUTED, linewidth=0.7, linestyle="--", zorder=0)
        ax.set_title(title, loc="left")
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
    axes[0].set_yticks(range(len(order)), order)
    axes[0].set_ylim(len(order) - 0.5, -0.5)
    axes[0].tick_params(axis="y", length=0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", ncol=4, bbox_to_anchor=(0.105, 0.968))
    fig.suptitle("Fufa differentiation score: knockdown vs NTC, 56 eligible targets", x=0.11, y=0.99,
                 ha="left", fontsize=9)
    fig.text(0.11, 0.012, "Lines: 95% stratified cell-bootstrap intervals (5,000 draws). Targets sorted by the "
             "DMSO mean shift.", fontsize=6.5, color="#52514e")
    save(fig, "phase2_aim1_forest_fufa")


def scatter(highlights, mean, disp):
    gene = mean[mean.kind.eq("gene") & mean.readout.eq(FUFA)][["target", "condition", "estimate", "ci_low",
                                                               "ci_high"]]
    lvr = disp[disp.kind.eq("gene") & disp.readout.eq(FUFA) & disp.statistic.eq("log_var_ratio")][
        ["target", "condition", "estimate", "q_value"]]
    data = gene.merge(lvr, on=["target", "condition"], suffixes=("_mean", "_lvr"))
    calls = highlights[highlights.readout.eq(FUFA) & highlights.primary_call.eq(True)]
    called = set(zip(calls.target, calls.condition))
    smallest = set(zip(*[data.sort_values("q_value").groupby("condition").head(3)[c]
                         for c in ("target", "condition")]))
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 85 / 25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.76, bottom=0.15, wspace=0.08)
    for ax, condition in zip(axes, ("DMSO", "Vem")):
        d = data[data.condition.eq(condition)]
        sig = d.q_value.lt(0.05)
        ax.axhline(0, color=MUTED, linewidth=0.7, linestyle="--", zorder=0)
        ax.axvline(0, color=MUTED, linewidth=0.7, linestyle="--", zorder=0)
        ax.scatter(d.estimate_mean[~sig], d.estimate_lvr[~sig], s=14, color=MUTED, edgecolor="white",
                   linewidth=0.5, label="log variance ratio q ≥ 0.05")
        ax.scatter(d.estimate_mean[sig], d.estimate_lvr[sig], s=16, color=COLORS[condition],
                   marker=MARKERS[condition], edgecolor="white", linewidth=0.5, label="log variance ratio q < 0.05")
        labelled = [row for row in d.sort_values("estimate_mean").itertuples()
                    if (row.target, condition) in called or (row.target, condition) in smallest]
        for k, row in enumerate(labelled):
            ax.annotate(row.target, (row.estimate_mean, row.estimate_lvr), xytext=(3, 2 if k % 2 else -8),
                        textcoords="offset points", fontsize=6, color=INK,
                        fontweight="bold" if (row.target, condition) in called else "normal")
        ax.set_title(condition, loc="left")
        ax.set_xlabel("Mean shift, knockdown − NTC (DMSO NTC SD)")
    axes[0].set_ylabel("log variance ratio")
    fig.suptitle("Mean shift vs dispersion change in the Fufa score", x=0.08, y=0.98, ha="left", fontsize=9)
    fig.text(0.08, 0.905, "Bold label: primary call (mean CI includes 0 and log variance ratio BH q < 0.05 "
             "within condition). Plain labels: three smallest log variance ratio q per condition.",
             fontsize=6.5, color="#52514e")
    handles = [plt.Line2D([], [], linestyle="", marker="o", color=MUTED, label="log variance ratio q ≥ 0.05")]
    handles += [plt.Line2D([], [], linestyle="", marker=MARKERS[c], color=COLORS[c],
                           label=f"{c}: log variance ratio q < 0.05") for c in ("DMSO", "Vem")]
    fig.legend(handles=handles, loc="upper left", ncol=3, bbox_to_anchor=(0.075, 0.885), fontsize=6.5)
    save(fig, "phase2_aim1_mean_vs_dispersion")


def ecdf(highlights, scores_path):
    if not scores_path.exists():
        print(f"skipping ECDF: {scores_path} not found")
        return
    h = highlights[highlights.readout.eq(FUFA)].copy()
    calls = h[h.primary_call.eq(True)]
    title = "Primary calls"
    if calls.empty:
        calls = h.sort_values("log_var_ratio_q").groupby("condition").head(3)
        title = "No primary calls; three smallest log-variance-ratio q per condition"
    scores = pd.read_csv(scores_path, usecols=["condition", "target", FUFA])
    panels = list(calls.itertuples())
    ncol = min(4, len(panels))
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(183 / 25.4, (48 * nrow + 20) / 25.4), squeeze=False,
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.98, top=1 - 14 / (48 * nrow + 20), bottom=11 / (48 * nrow + 20),
                        wspace=0.08, hspace=0.45)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    for ax, row in zip(axes.ravel(), panels):
        sub = scores[scores.condition.eq(row.condition)]
        for group, color, lw in (("NTC", MUTED, 1.4), (row.target, COLORS[row.condition], 1.4)):
            v = np.sort(sub.loc[sub.target.eq(group), FUFA].to_numpy())
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=color, linewidth=lw,
                    label=f"{group} (n={len(v)})")
        ax.set_title(f"{row.target}, {row.condition}", loc="left")
        ax.text(0.02, 0.97, f"Δ mean {row.mean_estimate:+.2f}\nlog VR {row.log_var_ratio_estimate:+.2f}\n"
                f"q = {row.log_var_ratio_q:.2g}", transform=ax.transAxes, va="top", fontsize=6, color="#52514e")
        ax.legend(loc="lower right", fontsize=6, frameon=True, facecolor="white", edgecolor="none", framealpha=0.9)
        ax.grid(color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
    for ax in axes[-1]:
        ax.set_xlabel("Fufa score (DMSO NTC SD)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative fraction")
    fig.suptitle(f"Fufa score distributions vs NTC. {title}", x=0.07, y=0.995, ha="left", fontsize=9)
    save(fig, "phase2_aim1_ecdf_fufa")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-name", default="phase2_aim1")
    args = parser.parse_args(argv)
    style()
    mean = pd.read_csv(SUMMARY / "phase2_mean_contrasts.csv")
    disp = pd.read_csv(SUMMARY / "phase2_dispersion_contrasts.csv")
    highlights = pd.read_csv(SUMMARY / "phase2_highlight_targets.csv")
    forest(mean, disp)
    scatter(highlights, mean, disp)
    ecdf(highlights, ROOT / "outputs" / args.run_name / "scores.csv.gz")
    print(f"wrote figures to {FIGURES}")


if __name__ == "__main__":
    sys.exit(main())
