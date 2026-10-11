"""Figure for Phase 2b: power of the primary dispersion call (reports/figures/phase2b_power_curves.{png,svg}).

Usage:
    python scripts/phase2b_figures.py

Reads results/summary/phase2b_power.csv and phase2_dispersion_contrasts.csv (aggregate tables only).
One panel per condition: power vs true log variance ratio. Scale-model lines by pseudo-target size
(ordinal blue ramp, steps 250 / 450 / 650 of the reference palette), solid for alpha 0.05 and dashed
for the Bonferroni level 0.05/56; mixture scenarios at the median size as open markers.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from aim1_figures import FUFA, GRID, INK, MUTED, SUMMARY, save, style  # noqa: E402

RAMP = {0.1: "#86b6ef", 0.5: "#2a78d6", 0.9: "#104281"}
ALPHAS = {"a05": ("-", "o", "alpha 0.05"), "bonferroni": ("--", "s", "alpha 0.05/56")}


def power_curves(power, disp):
    max_obs = (disp[disp.statistic.eq("log_var_ratio") & disp.kind.eq("gene") & disp.readout.eq(FUFA)]
               .assign(a=lambda d: d.estimate.abs()).groupby("condition").a.max())
    conditions = list(dict.fromkeys(power.condition))
    fig, axes = plt.subplots(1, len(conditions), figsize=(7.0, 3.0), sharey=True)
    for ax, condition in zip(axes, conditions):
        sub = power[power.condition.eq(condition)]
        sizes = []
        for q, color in RAMP.items():
            g = sub[sub.n_quantile.eq(q) & sub.model.isin(["no_injection", "scale"])].sort_values("true_lvr")
            sizes.append(str(int(g.n_target.iloc[0])))
            for label, (ls, _, _) in ALPHAS.items():
                ax.plot(g.true_lvr, g[f"power_lvr_{label}"], ls=ls, lw=2 if ls == "-" else 1.5, color=color,
                        solid_capstyle="round")
        mix = sub[sub.model.eq("mixture") & sub.n_quantile.eq(0.5)]
        for label, (_, marker, _) in ALPHAS.items():
            ax.scatter(mix.true_lvr, mix[f"power_lvr_{label}"], s=16, marker=marker, facecolors="white",
                       edgecolors=RAMP[0.5], linewidths=1.0, zorder=3)
        ax.axhline(0.8, color=MUTED, lw=0.8, ls=":")
        ax.axhline(0.05, color=GRID, lw=0.8)
        if condition in max_obs:
            ax.axvline(max_obs[condition], color=MUTED, lw=0.8, ls="-.")
        ax.set_title(f"{condition}  (n = {' / '.join(sizes)} cells)", loc="left")
        ax.set_xlabel("True log variance ratio (pseudo-target vs NTC)")
        ax.set_xlim(0, 1.6)
        ax.set_ylim(0, 1.02)
        ax.grid(axis="y", color=GRID, lw=0.5)
    axes[0].set_ylabel("Power (log variance ratio test)")
    handles = [Line2D([], [], color=INK, ls="-", lw=2, label="Scale model, alpha 0.05"),
               Line2D([], [], color=INK, ls="--", lw=1.5, label="Scale model, alpha 0.05/56"),
               Line2D([], [], color=RAMP[0.5], marker="o", ls="", mfc="white", label="Mixture, median n, 0.05"),
               Line2D([], [], color=RAMP[0.5], marker="s", ls="", mfc="white", label="Mixture, median n, 0.05/56")]
    handles += [Line2D([], [], color=c, lw=2, label=f"n = p{int(q * 100)} of target sizes") for q, c in RAMP.items()]
    handles += [Line2D([], [], color=MUTED, lw=0.8, ls="-.", label="Largest observed |lvr| (Phase 2)"),
                Line2D([], [], color=MUTED, lw=0.8, ls=":", label="80% power")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Power of the Phase 2 primary dispersion test (Fufa score, NTC pseudo-targets)", x=0.01,
                 ha="left", fontsize=8)
    fig.tight_layout(rect=(0, 0.17, 1, 0.95))
    save(fig, "phase2b_power_curves")


def main():
    style()
    power = pd.read_csv(SUMMARY / "phase2b_power.csv")
    disp = pd.read_csv(SUMMARY / "phase2_dispersion_contrasts.csv")
    power_curves(power, disp)


if __name__ == "__main__":
    main()
