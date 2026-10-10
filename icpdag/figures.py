"""Effect-size figure and Chinese report for the explicitly scoped ATE analysis."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


YAP = "Kanai2024_YAP_TAZ_up"
DIFFERENTIATION = "Fufa2019_AZD_up_SOX10_bound_pigmentation"


def effect_row(effects, kind, target, outcome, condition, comparator=None):
    selected = effects.loc[effects.kind.eq(kind) & effects.target.eq(target)
                           & effects.outcome.eq(outcome) & effects.condition.eq(condition)]
    if comparator is not None:
        selected = selected.loc[selected.comparator.eq(comparator)]
    if len(selected) != 1:
        raise ValueError("Report contrast must identify exactly one estimate")
    return selected.iloc[0]


def formatted(row):
    return f"{row.estimate:+.3f} [{row.ci_low:+.3f}, {row.ci_high:+.3f}]"


def draw_figure(root, effects, targets, bootstrap_repeats):
    """183-mm quantitative grid; all seven targets, no significance filtering."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "legend.frameon": False,
    })
    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 100 / 25.4))
    fig.subplots_adjust(left=0.105, right=0.985, top=0.80, bottom=0.24, wspace=0.63)
    fig.suptitle("Average intervention contrasts in observed cells", x=0.105, y=0.98, ha="left", fontsize=10)
    source_rows = []
    colors = {"DMSO": "#476A83", "Vem": "#B46D3B"}
    markers = {"DMSO": "o", "Vem": "s"}
    for ax, outcome, title in zip(axes[:2], [YAP, DIFFERENTIATION],
                                  ["YAP/TAZ program", "Differentiation program"]):
        for condition, offset in [("DMSO", -0.14), ("Vem", 0.14)]:
            data = pd.DataFrame([effect_row(effects, "gene", t, outcome, condition) for t in targets])
            # Percentile intervals are drawn directly, rather than assuming that
            # the point estimate must lie inside a bootstrap interval.
            ys = np.arange(len(targets)) + offset
            ax.hlines(ys, data.ci_low, data.ci_high, color=colors[condition], linewidth=1)
            ax.scatter(data.estimate, ys, color=colors[condition], marker=markers[condition], s=15,
                       label=condition, zorder=3)
            source_rows.append(data.assign(panel=title))
        ax.set_title(title, loc="left", pad=9)
        ax.set_xlabel("Guide assignment − NTC\n(DMSO NTC SD)")
    interaction = pd.DataFrame([effect_row(effects, "interaction", t, DIFFERENTIATION,
                                           "Vem_minus_DMSO") for t in targets])
    axes[2].hlines(np.arange(len(targets)), interaction.ci_low, interaction.ci_high, color="#625A70", linewidth=1)
    axes[2].scatter(interaction.estimate, np.arange(len(targets)), color="#625A70", marker="D", s=15, zorder=3)
    axes[2].set_title("Differentiation interaction", loc="left", pad=9)
    axes[2].set_xlabel("Gene × drug contrast\n(DMSO NTC SD)")
    source_rows.append(interaction.assign(panel="Differentiation interaction"))
    for letter, ax in zip("abc", axes):
        ax.axvline(0, linestyle="--", linewidth=0.7, color="#999999", zorder=0)
        ax.set_yticks(np.arange(len(targets)), targets)
        ax.set_ylim(len(targets)-0.5, -0.5)
        ax.tick_params(length=3, width=0.7)
        ax.text(-0.25, 1.075, letter, transform=ax.transAxes, fontsize=9, weight="bold")
        ax.margins(x=0.13)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.99, 0.944), ncol=2, fontsize=7)
    fig.text(0.105, 0.085, f"Points: observed mean contrasts. Lines: 95% cell-bootstrap intervals ({bootstrap_repeats:,} draws).", fontsize=7)
    fig.text(0.105, 0.04, "Frozen scores; intervals cover this cell pool only. Independent screen uncertainty is unavailable.", fontsize=7)
    fig.savefig(figures / "average-effects.svg")
    fig.savefig(figures / "average-effects.pdf")
    fig.savefig(figures / "average-effects.png", dpi=600)
    plt.close(fig)
    pd.concat(source_rows, ignore_index=True).to_csv(figures / "source_data.csv", index=False)
