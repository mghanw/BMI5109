"""Split-half reliability of the program scores and power of the primary dispersion call.

Usage:
    python scripts/power_analysis.py [--config config/power_analysis.json] [--splits S] [--reps R]
                                     [--bootstrap B] [--run-name NAME]

NTC cells only (no target data). Scores are recomputed with icpdag.program_scores exactly as in
the Aim 1 analysis and checked against outputs/aim1/scores.csv.gz. Step 2: for each readout, 200 seeded
random gene halves are scored with the same background procedure; the half-score correlation across
DMSO NTC and Vem NTC cells gives the Spearman-Brown reliability. Step 3: pseudo-targets of NTC cells
(sizes = 10th / 50th / 90th percentile of the eligible-target sizes) receive mean-preserving
heterogeneity and are tested with the Aim 1 primary call (icpdag/power.py).
Any reduced count (--splits / --reps / --bootstrap) is a smoke test: outputs go to outputs/<run>/ only.
"""
import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import aim1  # noqa: E402
from icpdag import dispersion as ds  # noqa: E402
from icpdag import gene_sets as gs  # noqa: E402
from icpdag import power as pw  # noqa: E402
from icpdag.program_scores import program_scores  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

LOG = logging.getLogger("power_analysis")
NTC = "NTC"


# ---------------------------------------------------------------- data
def load_ntc(acfg):
    """NTC cells only: log1p CP10k of the scoring genes, gene annotation and cell metadata."""
    counts, var_all, obs = aim1.load_data(acfg)
    ntc = obs.target.eq(NTC).to_numpy()
    keep = gs.scoring_genes(var_all)
    var = var_all[keep]
    obs = obs[ntc]
    X = aim1.log_cp10k(counts[ntc][:, keep], obs.total_counts.to_numpy(float))
    meta = obs[["condition", "target"]].copy()
    return X, var, meta


def regression_check(scores, path, tolerance):
    old = pd.read_csv(path, index_col=0)
    old = old.loc[scores.index, scores.columns]
    gap = float(np.abs(old.to_numpy() - scores.to_numpy()).max())
    if not gap < tolerance:
        raise ValueError(f"Recomputed NTC scores differ from the Aim 1 analysis (max abs diff {gap:.3g})")
    LOG.info("regression check: NTC scores match the Aim 1 analysis (max abs diff %.2e)", gap)
    return gap


# ---------------------------------------------------------------- step 2: reliability
def reliability(X, var, meta, programs, library_symbols, conditions, splits, seed, progress):
    symbols = var.symbol.to_numpy()
    lookup = set(symbols)
    annotation = var[["symbol", "gene_type"]]
    all_program_genes = sorted({g for info in programs.values() for g in info["genes"]})
    masks = {c: meta.condition.eq(c).to_numpy() for c in conditions}
    rows = []
    for readout, info in programs.items():
        genes = [g for g in info["genes"] if g in lookup]
        for s in range(splits):
            rng = ds.stratum_rng(seed, "split_half", readout, s)
            half_a, half_b = pw.split_halves(genes, rng)
            halves = {"half_a": {"genes": half_a}, "half_b": {"genes": half_b}}
            cfg = {"network_genes": all_program_genes, "gene_readouts": sorted(library_symbols),
                   "seed": int(rng.integers(2 ** 31))}
            adjusted, raw, _ = program_scores(X, symbols, annotation, meta, halves, cfg)
            for variant, frame in (("background_adjusted", adjusted), ("raw_no_background", raw)):
                for c, mask in masks.items():
                    a, b = frame.half_a.to_numpy()[mask], frame.half_b.to_numpy()[mask]
                    rows.append(dict(readout=readout, condition=c, variant=variant, split=s,
                                     n_genes=len(genes), n_half_a=len(half_a), n_half_b=len(half_b),
                                     n_cells=int(mask.sum()), pearson_r=float(np.corrcoef(a, b)[0, 1]),
                                     spearman_r=float(stats.spearmanr(a, b).statistic)))
        progress(f"reliability {readout}")
    per_split = pd.DataFrame(rows)
    per_split["reliability"] = pw.spearman_brown(per_split.pearson_r)
    keys = ["readout", "condition", "variant", "n_genes", "n_half_a", "n_half_b", "n_cells"]
    out = per_split.groupby(keys, sort=False).agg(
        n_splits=("split", "size"),
        r_median=("pearson_r", "median"), r_min=("pearson_r", "min"),
        r_p05=("pearson_r", lambda v: v.quantile(0.05)), r_p95=("pearson_r", lambda v: v.quantile(0.95)),
        r_max=("pearson_r", "max"), spearman_r_median=("spearman_r", "median"),
        reliability_median=("reliability", "median"), reliability_min=("reliability", "min"),
        reliability_p05=("reliability", lambda v: v.quantile(0.05)),
        reliability_p95=("reliability", lambda v: v.quantile(0.95)),
        reliability_max=("reliability", "max")).reset_index()
    out.insert(1, "readout_role", out.readout.map(gs.ROLES))
    return out


# ---------------------------------------------------------------- step 3: power
def target_sizes(path, readout, quantiles, conditions):
    d = pd.read_csv(path)
    d = d[d.statistic.eq("log_var_ratio") & d.kind.eq("gene") & d.readout.eq(readout)]
    out = {}
    for c in conditions:
        n = d[d.condition.eq(c)].n_target.to_numpy(float)
        out[c] = [(q, int(np.round(np.quantile(n, q)))) for q in quantiles]
    return out


def power_tables(scores, meta, sizes, pc, seed, progress):
    alphas = {"a05": pc["alpha"], "bonferroni": pc["alpha"] / pc["n_targets_tested"]}
    scenarios = pw.make_scenarios(pc["mixture_pi"], pc["mixture_d_sd"],
                                  pc["scale_true_lvr"] + pc["scale_true_lvr_extension"])
    floor = 2 / (pc["bootstrap"] + 1)
    reachable = {k: a for k, a in alphas.items() if floor < a}
    for label in set(alphas) - set(reachable):
        LOG.warning("B=%d cannot reach alpha %s (%.1e): its power is NaN (resolution-limited smoke run)",
                    pc["bootstrap"], label, alphas[label])
    frames = []
    for condition, ns in sizes.items():
        ntc = scores[pc["readout"]].to_numpy()[meta.condition.eq(condition).to_numpy()]
        for q, n in ns:
            res = pw.simulate(ntc, n, scenarios, pc["reps"], pc["bootstrap"], seed, condition, reachable.values())
            summary = pw.summarize(res, scenarios, alphas)
            for label in set(alphas) - set(reachable):
                summary.loc[:, summary.columns.str.endswith(f"_{label}")] = np.nan
            summary.insert(0, "condition", condition)
            summary.insert(1, "n_quantile", q)
            summary.insert(2, "n_target", n)
            summary.insert(3, "n_control", len(ntc) - n)
            summary.insert(4, "ntc_sd", float(np.std(ntc, ddof=1)))
            summary["bootstrap"] = pc["bootstrap"]
            frames.append(summary)
            progress(f"power {condition} n={n}")
    table = pd.concat(frames, ignore_index=True)
    table["scale_grid_extension"] = table.model.eq("scale") & table.true_lvr.isin(pc["scale_true_lvr_extension"])
    return table, alphas


def mde_table(power, alphas, target, rel):
    """MDE per condition x n x alpha x test; signal scale only when the split p05 reliability is > 0.

    rel: Fufa background-adjusted reliability rows indexed by condition (median, p05, p95 over splits).
    """
    rows = []
    for (c, q, n), g in power.groupby(["condition", "n_quantile", "n_target"]):
        curve = g[g.model.isin(["no_injection", "scale"])]
        rho = rel.loc[c]
        estimable = rho.reliability_p05 > 0
        for label in alphas:
            for kind in ("lvr", "call"):
                mde = pw.min_detectable(curve.true_lvr, curve[f"power_{kind}_{label}"], target)
                signal = {f"mde_lvr_signal_scale{suffix}": float(pw.signal_scale_lvr(mde, rho[col])) if estimable
                          else np.nan
                          for suffix, col in (("", "reliability_median"), ("_at_rel_p95", "reliability_p95"),
                                              ("_at_rel_p05", "reliability_p05"))}
                rows.append(dict(condition=c, n_quantile=q, n_target=n, alpha=label, test=kind, mde_lvr=mde,
                                 reliability_median=rho.reliability_median, reliability_p05=rho.reliability_p05,
                                 reliability_p95=rho.reliability_p95, signal_scale_estimable=bool(estimable),
                                 **signal))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- summary
def fmt(x, digits=2):
    return "n/a" if pd.isna(x) else f"{x:.{digits}f}"


def write_summary(path, rel, power, mde, alphas, pc, gap, max_obs, run):
    adj = rel[rel.variant.eq("background_adjusted")]
    raw = rel[rel.variant.eq("raw_no_background")].set_index(["readout", "condition"])
    lines = ["# Score reliability and power of the primary dispersion call", "",
             f"Run `outputs/{run}` (scripts/power_analysis.py, config/power_analysis.json). NTC cells only; no target "
             f"data were used. Recomputed NTC scores match the Aim 1 analysis (max abs diff {gap:.1e}).", "",
             "## Split-half reliability (NTC cells)",
             f"Program genes split at random into halves ({int(adj.n_splits.iloc[0])} seeded splits); each half "
             "is scored with the same expression-matched background procedure. Reliability = Spearman-Brown "
             "2r/(1+r) of the half-score correlation: the share of score variance reproducible across gene "
             "halves (biology plus any shared cell-level technical factor), the rest being gene-sampling noise.",
             "", "| Readout | Role | Genes | Condition | Reliability median [min, max] | Raw (no background) median |",
             "|---|---|---:|---|---|---:|"]
    for r in adj.itertuples():
        rr = raw.loc[(r.readout, r.condition), "reliability_median"]
        lines.append(f"| {r.readout} | {r.readout_role} | {r.n_genes} | {r.condition} | "
                     f"{fmt(r.reliability_median)} [{fmt(r.reliability_min)}, {fmt(r.reliability_max)}] | {fmt(rr)} |")
    fufa = adj[adj.readout.eq(pc["readout"])].set_index("condition")
    fraw = raw.xs(pc["readout"], level="readout")
    rho = fufa.reliability_median
    lines += ["", f"**Fufa (primary).** Background-adjusted reliability is "
              + " and ".join(f"{fmt(rho[c])} ({c})" for c in rho.index)
              + "; without the background it is "
              + " and ".join(f"{fmt(fraw.loc[c, 'reliability_median'])} ({c})" for c in rho.index)
              + ". The two gene halves of the raw score share a cell-level component, but the expression-matched "
              "background removes it, and what remains does not replicate across gene halves. In NTC cells, "
              "nearly all between-cell variance of the Fufa score as defined is therefore gene-sampling noise. "
              "Population means still carry signal (Aim 1 mean contrasts); single-cell dispersion of this "
              "score carries little."]
    bonf = alphas["bonferroni"]
    lines += ["", "## Power of the primary dispersion call (Fufa)",
              f"Pseudo-targets of n NTC cells vs the remaining NTC cells; {pc['reps']} replicates per scenario, "
              f"B = {pc['bootstrap']} bootstrap draws (Aim 1: 5000). Test: log variance ratio bootstrap p "
              f"< alpha (\"lvr\"), and the full primary call that also requires the mean CI to include 0 "
              f"(\"call\"). alpha 0.05 is an optimistic bound for BH across 56 targets; 0.05/56 = {bonf:.1e} "
              "(Bonferroni) is the bound for a single true hit.", ""]
    null = power[power.model.eq("no_injection")]
    lines.append("Null rejection rates (lvr test, alpha 0.05): " + ", ".join(
        f"{r.condition} n={r.n_target}: {fmt(r.power_lvr_a05, 3)}" for r in null.itertuples())
        + f" (mean {fmt(null.power_lvr_a05.mean(), 3)}; Monte Carlo SE about "
        f"{fmt(np.sqrt(0.05 * 0.95 / pc['reps']), 3)} each). The percentile bootstrap test is close to nominal.")
    lines += ["", "### Minimal detectable true log variance ratio at 80% power (scale model)",
              "Linear interpolation on the scale-model grid (0, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5). Signal scale: "
              "the implied log variance ratio of the reproducible signal, log(1 + (e^MDE - 1)/reliability), "
              "assuming equal noise variance in target and NTC cells; median reliability, with the range at the "
              "split p95 and p05 reliabilities. It is not estimable when the split p05 reliability is <= 0.", "",
              "| Condition | n (quantile) | MDE, alpha 0.05 | MDE, Bonferroni | Full call, Bonferroni | "
              "Bonferroni MDE, signal scale |", "|---|---|---:|---:|---:|---|"]
    m = mde.set_index(["condition", "n_target", "alpha", "test"])
    for (c, q, n), _ in mde.groupby(["condition", "n_quantile", "n_target"], sort=False):
        b = m.loc[(c, n, "bonferroni", "lvr")]
        signal = (f"{fmt(b.mde_lvr_signal_scale)} [{fmt(b.mde_lvr_signal_scale_at_rel_p95)}, "
                  f"{fmt(b.mde_lvr_signal_scale_at_rel_p05)}]" if b.signal_scale_estimable
                  else "not estimable (reliability ~ 0)")
        lines.append(f"| {c} | {n} (p{int(round(q * 100))}) | {fmt(m.loc[(c, n, 'a05', 'lvr'), 'mde_lvr'])} | "
                     f"{fmt(b.mde_lvr)} | {fmt(m.loc[(c, n, 'bonferroni', 'call'), 'mde_lvr'])} | {signal} |")
    lines += ["", "### Mixture scenarios on the same scale (median n)",
              "A fraction pi of cells shifted by +d and pi by -d (d in NTC SD units); true log variance ratio "
              "log(1 + 2 pi d^2). Cells: true lvr, then power of the lvr test at alpha 0.05 / Bonferroni.", ""]
    med = power[power.model.eq("mixture") & power.n_quantile.eq(0.5)]
    d_values = sorted(med.d_sd.unique())
    lines += ["| Condition | pi | " + " | ".join(f"d = {d:g}" for d in d_values) + " |",
              "|---|---:|" + "---|" * len(d_values)]
    for (c, pi), g in med.groupby(["condition", "pi"], sort=False):
        g = g.set_index("d_sd")
        cells = [f"{fmt(g.loc[d, 'true_lvr'])}: {fmt(g.loc[d, 'power_lvr_a05'])} / "
                 f"{fmt(g.loc[d, 'power_lvr_bonferroni'])}" for d in d_values]
        lines.append(f"| {c} | {pi:g} | " + " | ".join(cells) + " |")
    med_mde = mde[mde.n_quantile.eq(0.5) & mde.test.eq("lvr") & mde.alpha.eq("bonferroni")]
    lines += ["", "## How to report the Aim 1 null",
              f"Largest observed |Fufa log variance ratio| in Aim 1: "
              + ", ".join(f"{c} {fmt(v)}" for c, v in max_obs.items()) + ".", "",
              "Suggested wording: \"No knockdown changed the between-cell dispersion of the Fufa score without a "
              "mean shift. At the median target size the design had 80% power, after multiplicity control, only "
              "for log variance ratios of about "
              + " / ".join(f"{fmt(r.mde_lvr)} ({r.condition})" for r in med_mde.itertuples())
              + ", and the score's split-half reliability in NTC cells is only "
              + " / ".join(f"{fmt(rho[c])} ({c})" for c in rho.index)
              + ". The null result therefore means that heterogeneity of the Fufa response is not measurable "
              "with this single-cell score; it is not evidence that the response is homogeneous.\"",
              "", "Figure: `reports/figures/power_curves.png`. Tables: `score_reliability.csv`, "
              "`power_simulation.csv`."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------- main
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "power_analysis.json"))
    parser.add_argument("--splits", type=int)
    parser.add_argument("--reps", type=int)
    parser.add_argument("--bootstrap", type=int)
    parser.add_argument("--run-name")
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    with open(ROOT / cfg["aim1_config"], encoding="utf-8") as f:
        acfg = json.load(f)
    pc = cfg["power"]
    smoke = any(v is not None for v in (args.splits, args.reps, args.bootstrap))
    splits = args.splits or cfg["reliability"]["splits"]
    pc["reps"] = args.reps or pc["reps"]
    pc["bootstrap"] = args.bootstrap or pc["bootstrap"]
    run = args.run_name or ("smoke_power" if smoke else "power_analysis")
    run_dir = ROOT / acfg["paths"]["outputs_dir"] / run
    run_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(run_dir / "log.txt", mode="w", encoding="utf-8")])
    with open(run_dir / "config.json", "w", encoding="utf-8") as f:
        json.dump(dict(cfg, splits=splits, run_name=run, smoke=smoke, argv=sys.argv[1:]), f, indent=2)
    t0 = time.time()
    progress = lambda msg: LOG.info("%s (%.0fs)", msg, time.time() - t0)
    out_dir = run_dir if smoke else ROOT / acfg["paths"]["summary_dir"]
    seed = cfg["seed"]
    conditions = tuple(acfg["conditions"])

    X, var, meta = load_ntc(acfg)
    progress(f"loaded {X.shape[0]} NTC cells x {X.shape[1]} scoring genes")
    eligible = sorted(load_split(ROOT / acfg["paths"]["target_split"]).target)
    library_symbols, _, _ = aim1.library_genes(acfg, var, eligible)
    _, mapping = aim1.load_programs(ROOT / acfg["paths"]["programs"])
    programs = gs.programs_for_scoring(mapping)
    scores, _, _ = aim1.score_cells(X, var, meta, programs, library_symbols, acfg["seed"])
    gap = regression_check(scores, ROOT / cfg["aim1_scores"], cfg["regression_tolerance"])

    rel = reliability(X, var, meta, programs, library_symbols, conditions, splits, seed, progress)
    aim1.assert_no_cell_ids(rel)
    rel.to_csv(out_dir / "score_reliability.csv", index=False)

    sizes = target_sizes(ROOT / cfg["aim1_dispersion"], pc["readout"], pc["n_quantiles"], conditions)
    LOG.info("pseudo-target sizes: %s", sizes)
    power, alphas = power_tables(scores, meta, sizes, pc, seed, progress)
    aim1.assert_no_cell_ids(power)
    power.to_csv(out_dir / "power_simulation.csv", index=False)

    fufa = rel[rel.readout.eq(pc["readout"]) & rel.variant.eq("background_adjusted")].set_index("condition")
    mde = mde_table(power, alphas, pc["target_power"], fufa)
    mde.to_csv(run_dir / "mde.csv", index=False)
    disp = pd.read_csv(ROOT / cfg["aim1_dispersion"])
    disp = disp[disp.statistic.eq("log_var_ratio") & disp.kind.eq("gene") & disp.readout.eq(pc["readout"])]
    max_obs = disp.assign(a=disp.estimate.abs()).groupby("condition").a.max().to_dict()
    write_summary(out_dir / "power_summary.md", rel, power, mde, alphas, pc, gap, max_obs, run)
    progress("done")


if __name__ == "__main__":
    main()
