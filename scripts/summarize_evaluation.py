"""Summary and figures of the held-out evaluation (reads the tables of scripts/evaluate_heldout.py).

Usage:
    python scripts/summarize_evaluation.py            # results/summary/evaluation_summary.md, reports/figures/
    python scripts/summarize_evaluation.py --smoke    # outputs/evaluation_smoke/ only

Writes results/summary/evaluation_summary.md, reports/figures/evaluation_paired_z.png and
reports/figures/evaluation_r_vs_noise_ceiling.png. The post hoc magnitude paragraph (not pre-registered,
no test) is computed from the stored predictions and outputs/evaluation/observed_changes.npz; its
medians are also written to outputs/evaluation/post_hoc_downstream_magnitude.csv. Nothing is printed
except the written paths.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag import heldout_evaluation as ev  # noqa: E402
from icpdag.sem_inputs import panel_nodes, target_nodes  # noqa: E402
from icpdag.sem_prediction import MODELS  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

PRIMARY, REFERENCE = "cell_specific_cscn_sem", "shared_cscn_sem"
LABEL = {"cell_specific_cscn_sem": "Cell-specific CSCN-SEM", "shared_cscn_sem": "Shared CSCN-SEM",
         "cluster_level_sem": "Cluster-level SEM", "local_association_network_sem": "Local association SEM",
         "unweighted_propagation_ablation": "Unweighted propagation", "zero_downstream_baseline": "Zero downstream"}
COLORS = {PRIMARY: "#2F6690", REFERENCE: "#B46D3B"}


def f(x, digits=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{digits}f}"


def p_fmt(p):
    if not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"


def e(x):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2e}"


def test_line(row, back=True, sci=False):
    if back:
        effect = (f"HL {f(row.hl)} on the z scale (tanh {f(row.hl_back)}; 95% bootstrap CI "
                  f"[{f(row.hl_ci_low_back)}, {f(row.hl_ci_high_back)}])")
    elif sci:
        effect = f"HL {e(row.hl)} (95% bootstrap CI [{e(row.hl_ci_low)}, {e(row.hl_ci_high)}])"
    else:
        effect = f"HL {f(row.hl, 4)} (95% bootstrap CI [{f(row.hl_ci_low, 4)}, {f(row.hl_ci_high, 4)}])"
    return (f"n = {int(row.n_targets)}, NaN pairs dropped {int(row.n_nan_dropped)}, zero differences dropped "
            f"{int(row.n_zero_dropped)}, n_eff = {int(row.n_effective)}, {row.method}; W = {f(row.statistic, 1)}, "
            f"p = {p_fmt(row.p_value)}; {effect}")


def medians_table(prim, section, extra_col=None):
    sub = prim[prim.section.eq(section)]
    lines = []
    keys = [extra_col] if extra_col else []
    head = "| " + " | ".join([*(k.capitalize() for k in keys), "Model", "Condition", "Median r (NaN excluded, n)",
                               "Median r (r = 0, n)", "Median RMSE"]) + " |"
    lines += [head, "|" + "---|" * (len(keys) + 5)]
    for key, s in sub.groupby([*keys, "model", "condition"], sort=False):
        key = key if isinstance(key, tuple) else (key,)
        get = lambda v, m: s[s.version.eq(v)][m]  # noqa: E731
        main, zero, err = get("nan_excluded", "median_r"), get("r_zero", "median_r"), get("all", "median_rmse")
        n_main, n_zero = get("nan_excluded", "n_targets"), get("r_zero", "n_targets")
        cell = lambda v, n: f"{f(v.iloc[0])} ({int(n.iloc[0])})" if len(v) else "not computed"  # noqa: E731
        name = LABEL[key[len(keys)]]
        lead = [str(int(k)) if isinstance(k, float) and k.is_integer() else str(k) for k in key[:len(keys)]]
        lines.append("| " + " | ".join([*lead, name, key[-1], cell(main, n_main), cell(zero, n_zero),
                                         f(err.iloc[0], 4)]) + " |")
    return lines


def downstream_magnitude(predictions, observed, conditions, targets, nodes, target_node, models=MODELS):
    """Post hoc, per (model, target, condition, direction): RMS of the predicted downstream change over the
    evaluation genes, its ratio to |delta| and to the RMS of the observed downstream change.

    observed: array (condition, target, direction, gene) in the order of conditions and targets.
    """
    rows = []
    for ci, condition in enumerate(conditions):
        p = predictions[condition]
        for ti, target in enumerate(targets):
            t = p["targets"].index(target)
            mask = ev.evaluation_mask(nodes, target_node[target])
            for c in range(2):
                o = np.asarray(observed[ci, ti, c, mask], float)
                rms_obs, delta = float(np.sqrt(np.mean(o ** 2))), float(abs(p["delta"][t, c]))
                for i, model in enumerate(models):
                    rms = float(np.sqrt(np.mean(np.asarray(p["scaled"][i, t, c, mask], float) ** 2)))
                    rows.append(dict(model=model, target=target, condition=condition, direction=c, rms_predicted=rms,
                                     abs_delta=delta, rms_observed=rms_obs, ratio_to_delta=rms / delta,
                                     ratio_to_observed=rms / rms_obs))
    return pd.DataFrame(rows)


def magnitude_summary(units, fufa):
    """Per model medians over targets, conditions and directions, with the median |Fufa KD - NTC| shifts."""
    cols = ["rms_predicted", "ratio_to_delta", "ratio_to_observed", "abs_delta", "rms_observed"]
    out = units.groupby("model", sort=False)[cols].median()
    out["n_units"] = units.groupby("model", sort=False).size()
    kd = fufa[fufa.section.eq("target") & fufa.contrast.astype(str).str.startswith("KD_minus_NTC")
              & fufa.direction.astype(str).isin(["0", "1"])]
    out["fufa_abs_predicted"] = kd.groupby("model")["predicted"].apply(lambda s: s.abs().median())
    out["fufa_abs_observed"] = kd.groupby("model")["observed"].apply(lambda s: s.abs().median())
    return out.reset_index()


def post_hoc_lines(mag):
    mag = mag.set_index("model")
    fitted = mag.drop(index="zero_downstream_baseline", errors="ignore")
    worst = float(fitted.ratio_to_observed.max())
    first = mag.iloc[0]
    L = ["## Post hoc description (not pre-registered)", "",
         "Descriptive only; no test. Computed from the stored predictions and the observed change vectors over the "
         "249 evaluation genes: for each model, target, condition and direction, the root mean square (RMS) of the "
         "predicted downstream change, relative to the knockdown size |delta| and to the RMS of the observed "
         f"downstream change; medians over the {int(first.n_units)} target x condition x direction units. For scale, "
         f"the median |delta| is {f(float(first.abs_delta))} and the median observed RMS is "
         f"{f(float(first.rms_observed))} (log1p CP10k); the median absolute observed Fufa KD - NTC shift is "
         f"{f(float(first.fufa_abs_observed))} (DMSO NTC SD).", "",
         "| Model | Median predicted RMS | Ratio to abs(delta) | Ratio to observed RMS | "
         "Median abs. predicted Fufa KD - NTC shift |",
         "|---|---|---|---|---|"]
    for model, r in mag.iterrows():
        L.append(f"| {LABEL[model]} | {e(r.rms_predicted)} | {e(r.ratio_to_delta)} | {e(r.ratio_to_observed)} | "
                 f"{e(r.fufa_abs_predicted)} |")
    L.append("")
    if worst < 0.05:
        L.append(f"- The predicted downstream changes of every model are much smaller than the observed downstream "
                 f"changes (median ratio to the observed RMS at most {worst:.1%}), and the predicted Fufa shifts are "
                 "much smaller than the observed ones. Each model's RMSE is therefore essentially the RMS of the "
                 "observed change, which is the RMSE of the zero-downstream baseline, and each model's Fufa error is "
                 "essentially the observed Fufa shift. This explains why every model is close to the zero baseline "
                 "on RMSE and on the Fufa endpoint. Pearson r does not depend on the scale of the prediction, so the "
                 "small magnitude does not by itself explain the near-zero correlations.")
    else:
        L.append(f"- The largest median ratio of predicted to observed downstream RMS is {worst:.1%}.")
    return L


def write_summary(tables, run, out_path, lam, magnitude=None):
    m, prim, sec, fufa, dist = (tables[k] for k in ("metrics", "primary", "secondary", "fufa", "distribution"))
    tests = prim[prim.section.eq("test")].set_index("analysis")
    pr = tests.loc["primary_decision"]
    cov = prim[prim.section.eq("coverage")]
    flag = bool(prim[prim.section.eq("nan_flag")].nan_flag.iloc[0])
    L = ["# Held-out evaluation of the CSCN-SEM knockdown predictions", ""]
    L += [f"Run: {run.get('date', 'NA')}, code at commit {run.get('head', 'NA')[:7] if run.get('head') else 'NA'}; "
          "rules in `config/evaluation_run.json` (frozen before the run), endpoints in "
          "`config/evaluation_endpoints.json` (authoritative).", ""]
    L += ["## Deviations", ""]
    if run.get("deviation"):
        L += [f"- {run['deviation']} (this run). Earlier run(s): "
              + "; ".join(f"{h.get('date', 'NA')} at {str(h.get('head', 'NA'))[:7]}" for h in run.get("history", []))
              + ". Both results are documented here."]
    else:
        L += ["- None: the evaluation was run once with the committed code."]
    L += ["", "## Design", "",
          "- Observed response: target - NTC mean of log1p CP10k (whole-transcriptome totals) over the 249 panel genes "
          "other than the perturbed gene, on the evaluation cells: the prediction calibrated on calib_fold c "
          "(delta from the held-out target's own gene on calib_fold c) is evaluated on calib_fold 1 - c, for both c.",
          "- Predicted means average over all NTC reference cells of the condition (both calib folds), while the "
          "delta uses only calib_fold c cells; NTC reference cells carry no target response.",
          "- Per target: Pearson r and RMSE per condition and direction; directions and then conditions combined on "
          "the Fisher z scale (atanh, |r| clipped at 0.999999999999); RMSE averaged.",
          "- Amendment (2026-10-21, before model fitting): in the primary decision test and every paired test on z, "
          "a constant prediction gets r = 0 (z = 0), so all 56 targets enter; the pre-registered complete-pairs "
          "analysis is supporting. Wilcoxon zero differences are dropped (zero_method 'wilcox') and counted.",
          "- Sign convention: differences are cell-specific minus reference; a positive z difference or a negative "
          "RMSE difference favours the cell-specific model. Hodges-Lehmann (HL) estimates include zero differences; "
          "CIs are percentile bootstraps over targets (5,000 draws, seed 20261011).", ""]
    L += ["## Primary decision test (cell-specific vs shared CSCN-SEM, amended rule)", "",
          f"- {test_line(pr)}.",
          f"- Decision at alpha {f(pr.alpha, 2)}: H0 (no difference in held-out Fisher z) is "
          f"{'rejected' if pr.reject_h0 else 'not rejected'}.",
          "- Median back-transformed r over the paired targets: cell-specific "
          f"{f(np.tanh(pr.median_model))}, shared {f(np.tanh(pr.median_reference))}.", ""]
    cp = tests.loc["complete_pairs"]
    L += ["## Supporting analyses", "",
          "### Coverage (share of targets with a non-constant prediction) and accuracy (complete pairs)", "",
          "| Model | Non-constant DMSO | Non-constant Vem | Non-constant in both | NaN pooled z (share) | > 10% flag |",
          "|---|---|---|---|---|---|"]
    for model, s in cov.groupby("model", sort=False):
        s = s.set_index("condition")
        both = s.loc["both"]
        nan = f"{int(both.n_nan)} ({f(both.nan_share, 2)})" if np.isfinite(both.get("n_nan", np.nan)) else "r not computed"
        fl = ("yes" if both.nan_flag else "no") if isinstance(both.get("nan_flag"), (bool, np.bool_)) else "-"
        L.append(f"| {LABEL[model]} | {int(s.loc['DMSO'].n_nonconstant)}/56 | {int(s.loc['Vem'].n_nonconstant)}/56 | "
                 f"{int(both.n_nonconstant)}/56 | {nan} | {fl} |")
    L += ["", f"- Complete pairs (pre-registered NaN rule, supporting): {test_line(cp)}.",
          "- Coverage and accuracy are separate questions: the amended primary test credits a model for predicting "
          "a non-constant change for more targets (constant predictions score r = 0), whereas the complete-pairs "
          "test compares accuracy only on targets where both models predict a non-constant change.",
          f"- NaN flag (> 10% of targets with undefined pooled z under the pre-registered rule, any model with r): "
          f"{'FLAGGED' if flag else 'not flagged'}.", ""]
    L += ["### Per-condition and RMSE tests", ""]
    for c in ("DMSO", "Vem"):
        L.append(f"- {c} (amended z): {test_line(tests.loc[f'per_condition_{c}'])}.")
    L += [f"- RMSE (mean over conditions): {test_line(tests.loc['rmse'], back=False, sci=True)}.", ""]
    L += ["### Medians over targets", ""] + medians_table(prim, "median") + [""]
    L += ["### Per-fold medians (pooled over conditions)", ""] + medians_table(prim, "fold_median", "fold") + [""]
    L += ["### By connectivity stratum (descriptive; connected = an edge in >= 10% of NTC cells in both conditions)", ""]
    L += medians_table(prim, "stratum_median", "stratum") + [""]
    joint = m[m.direction.eq("combined") & m.condition.isin(["DMSO", "Vem"])]
    L += ["### Split-half noise ceiling (context only; not part of any test)", ""]
    for c in ("DMSO", "Vem"):
        s = joint[joint.model.eq(PRIMARY) & joint.condition.eq(c)]
        q = s.noise_ceiling_r.quantile([0.25, 0.5, 0.75])
        above = []
        for model in (PRIMARY, REFERENCE):
            t = joint[joint.model.eq(model) & joint.condition.eq(c)]
            t = t[~t.constant_prediction.astype(bool)]
            above.append(f"{LABEL[model]} {int((t.r > t.noise_ceiling_r).sum())} of {len(t)}")
        L.append(f"- {c}: median noise ceiling {f(q[0.5])} (IQR {f(q[0.25])} to {f(q[0.75])}); targets with a ceiling "
                 f"< 0: {int((s.noise_ceiling_r < 0).sum())} of {len(s)}.")
        L.append(f"  - Descriptive, not part of any test: targets whose held-out r exceeds their noise ceiling, among "
                 f"targets with a non-constant prediction: {'; '.join(above)}.")
    L += ["", "## Secondary comparisons (cell-specific vs each comparator; descriptive, Holm within each metric family)", "",
          "| Comparator | Metric | Version | n_eff (zeros dropped) | Method | p | Holm p | HL [95% CI] |",
          "|---|---|---|---|---|---|---|---|"]
    for row in sec.itertuples():
        back = row.metric == "z"
        hl = (f"{f(row.hl_back)} [{f(row.hl_ci_low_back)}, {f(row.hl_ci_high_back)}] (tanh)" if back
              else f"{e(row.hl)} [{e(row.hl_ci_low)}, {e(row.hl_ci_high)}]")
        L.append(f"| {LABEL[row.reference]} | {row.metric} | {row.version} ({row.role}) | {int(row.n_effective)} "
                 f"({int(row.n_zero_dropped)}) | {row.method} | {p_fmt(row.p_value)} | {p_fmt(row.p_holm)} | {hl} |")
    L += ["", "## Secondary endpoint: Fufa score mean-effect error", "",
          "Observed: Aim 1 per-cell scores, target - NTC on the evaluation fold; predicted: mean stored delta-scaled "
          "Fufa shift over included reference cells (genes outside the panel change by 0 for every model). Error "
          "|predicted - observed| per direction, averaged over directions.", "",
          "| Contrast | " + " | ".join(LABEL[k] for k in LABEL) + " |", "|---|" + "---|" * len(LABEL)]
    med = fufa[fufa.section.eq("median")]
    for contrast, s in med.groupby("contrast", sort=False):
        s = s.set_index("model").median_abs_error
        L.append(f"| {contrast} | " + " | ".join(f(s[k]) for k in LABEL) + " |")
    L += [""]
    for row in fufa[fufa.section.eq("test")].itertuples():
        L.append(f"- {row.contrast}, cell-specific - shared absolute error: {test_line(row, back=False, sci=True)} "
                 "(no multiplicity correction pre-specified).")
    L += ["", "## Tertiary endpoints (expected low power)", "",
          "Aim 1 found 0 primary dispersion calls and the Fufa score has low single-cell reliability, so these "
          "endpoints have low power (results/summary/power_summary.md). Predicted distribution: evaluation-fold NTC "
          "scores plus each cell's own stored shift (shared and zero-downstream models: a constant shift, so their "
          "predicted dispersion change is exactly 0); the predicted statistic compares the shifted with the same "
          "unshifted cells.", "",
          "| Model | Median abs. error: log var ratio | IQR difference | Centred W1 | Over-dispersion rate | Coverage |",
          "|---|---|---|---|---|---|"]
    both = dist[dist.condition.eq("both")]
    for model in LABEL:
        s = both[both.model.eq(model)]
        g = lambda st: s[s.section.eq("median") & s.statistic.eq(st)].median_abs_error.iloc[0]  # noqa: E731
        cal = s[s.section.eq("calibration")].iloc[0]
        L.append(f"| {LABEL[model]} | {f(g('log_var_ratio'))} | {f(g('iqr_diff'))} | {f(g('centered_wasserstein'))} | "
                 f"{f(cal.over_dispersion_rate)} | {f(cal.coverage)} |")
    drop = dist[dist.section.eq("dropped_cells")]
    L += ["", "- Over-dispersion: predicted log variance ratio above the upper bound of the observed 95% bootstrap "
          "interval (evaluation cells), averaged over directions per target x condition row (112 rows); nominal "
          "<= 0.025. Aim 1 all-cell intervals are listed alongside in evaluation_distribution.csv (columns only).",
          "- Evaluation NTC cells without a defined shift (dropped, summed over targets): "
          + "; ".join(f"{LABEL[r.model]} {r.condition} c={int(float(r.direction))}: {int(r.n_dropped)}"
                      for r in drop[drop.n_dropped > 0].itertuples()) if (drop.n_dropped > 0).any()
          else "- No evaluation NTC cell was dropped for an undefined shift."]
    if magnitude is not None:
        L += [""] + post_hoc_lines(magnitude)
    L += ["", "## Interpretation notes and limitations", "",
          f"- The cell-specific vs shared contrast is mainly structural: the NTC-only cross-validation selected "
          f"lambda {lam['DMSO']} (DMSO) and {lam['Vem']} (Vem), the largest grid value, and pure pooled coefficients "
          "(lambda = infinity, reference only) fitted marginally better, so local coefficients add little beyond "
          "per-cell graph structure.",
          "- The shared and cluster graphs use a population Spearman partial-correlation PC search, whereas the "
          "cell-specific graphs use the local CSCN count test: the comparison differs in CI-test family as well as "
          "in locality. Coverage matching equalises the number of targets with an edge, not the number of edges.",
          "- Weak target connectivity: connectivity follows NTC detection of the target gene; many targets have few "
          "or no downstream edges, and 40 targets are constant under the shared model in at least one condition.",
          "- Low score reliability (Fufa split-half 0.01 / 0.12) limits the distribution-level endpoints.",
          "- 70 DMSO cells (3.5%) are excluded from the cell-specific and ablation pseudobulk means because their "
          "stored graphs contain a cycle among fixed edges.",
          "- Fallback completions (relaxed rule) may introduce v-structures not supported by the data; the local "
          "association model has no separating sets and no collider step.", ""]
    out_path.write_text("\n".join(L), encoding="utf-8")


def post_hoc_magnitude(base, fufa):
    """Load the verified stored predictions and the observed vectors, then summarise the magnitudes."""
    run_cfg = json.loads((ROOT / "config" / "evaluation_run.json").read_text(encoding="utf-8"))
    sem = json.loads((ROOT / "config" / "sem_prediction.json").read_text(encoding="utf-8"))
    split = load_split(ROOT / sem["inputs"]["target_split"])
    panel = json.loads((ROOT / sem["inputs"]["panel"]).read_text(encoding="utf-8"))
    _, nodes = panel_nodes(panel)
    preds = ev.verify_predictions(ROOT / "outputs" / "sem_prediction" / "predictions",
                                  run_cfg["inputs"]["prediction_config_sha256"], split, nodes, MODELS,
                                  run_cfg["conditions"])
    with np.load(base / "observed_changes.npz") as z:
        if [str(n) for n in z["nodes"]] != list(nodes):
            raise ValueError("Observed change vectors use another gene order")
        conditions = [str(c) for c in z["conditions"]]
        targets = [str(t) for t in z["targets"]]
        observed = z["observed"]
    units = downstream_magnitude(preds, observed, conditions, targets, nodes, target_nodes(panel))
    return magnitude_summary(units, fufa)


def sem_lambda(condition):
    path = ROOT / "outputs" / "sem_prediction" / "lambda_cv" / f"{condition}.json"
    return json.loads(path.read_text(encoding="utf-8"))["lambda_selected"]


def style():
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 7,
                         "axes.titlesize": 8, "axes.labelsize": 7, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 0.7, "legend.frameon": False})


def paired_z_figure(m, path):
    pooled = m[m.condition.eq("pooled") & m.direction.eq("combined")]
    a = pooled[pooled.model.eq(PRIMARY)].set_index("target")
    b = pooled[pooled.model.eq(REFERENCE)].set_index("target").reindex(a.index)
    fig, ax = plt.subplots(figsize=(90 / 25.4, 85 / 25.4))
    cats = [("Both non-constant", ~a.constant_prediction & ~b.constant_prediction, "o", "#2F6690", True),
            ("Shared constant (r = 0)", ~a.constant_prediction & b.constant_prediction, "s", "#B46D3B", False),
            ("Cell-specific constant (r = 0)", a.constant_prediction & ~b.constant_prediction, "^", "#625A70", False),
            ("Both constant (r = 0)", a.constant_prediction & b.constant_prediction, "x", "#777777", True)]
    for label, sel, marker, color, filled in cats:
        if sel.any():
            kw = dict(color=color) if filled else dict(facecolors="none", edgecolors=color)
            ax.scatter(b.z_amended[sel], a.z_amended[sel], marker=marker, s=16, linewidths=0.8, label=label, **kw)
    lim = np.nanmax(np.abs(np.r_[a.z_amended, b.z_amended, 0.1])) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "--", color="#999999", linewidth=0.7, zorder=0)
    ax.axhline(0, color="#cccccc", linewidth=0.5, zorder=0)
    ax.axvline(0, color="#cccccc", linewidth=0.5, zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Shared CSCN-SEM: Fisher z (mean of DMSO and Vem)")
    ax.set_ylabel("Cell-specific CSCN-SEM: Fisher z")
    ax.set_title("Held-out targets: paired Fisher z (amended rule)", loc="left")
    ax.legend(loc="upper left", fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def ceiling_figure(m, path):
    fig, axes = plt.subplots(1, 2, figsize=(160 / 25.4, 75 / 25.4), sharex=True, sharey=True)
    joint = m[m.direction.eq("combined")]
    for ax, condition in zip(axes, ("DMSO", "Vem")):
        for model, marker in ((PRIMARY, "o"), (REFERENCE, "s")):
            s = joint[joint.model.eq(model) & joint.condition.eq(condition)]
            const = s.constant_prediction.astype(bool)
            ax.scatter(s.noise_ceiling_r[~const], s.r_amended[~const], marker=marker, s=14, color=COLORS[model],
                       label=f"{LABEL[model]}")
            ax.scatter(s.noise_ceiling_r[const], s.r_amended[const], marker=marker, s=14, facecolors="none",
                       edgecolors=COLORS[model], linewidths=0.8, label=f"{LABEL[model]}, constant (r = 0)")
        ax.plot([-1, 1], [-1, 1], "--", color="#999999", linewidth=0.7, zorder=0)
        ax.axhline(0, color="#cccccc", linewidth=0.5, zorder=0)
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-1.02, 1.02)
        ax.set_title(condition, loc="left")
        ax.set_xlabel("Split-half noise ceiling r (Spearman-Brown)")
    axes[0].set_ylabel("Held-out prediction r (directions combined)")
    axes[1].legend(loc="lower right", fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(argv)
    base = ROOT / "outputs" / ("evaluation_smoke" if args.smoke else "evaluation")
    tables_dir = base / "tables" if args.smoke else ROOT / "results" / "summary"
    md = base / "evaluation_summary.md" if args.smoke else tables_dir / "evaluation_summary.md"
    figures = base if args.smoke else ROOT / "reports" / "figures"
    run = json.loads((base / "run.json").read_text(encoding="utf-8"))
    names = dict(metrics="evaluation_target_metrics.csv", primary="evaluation_primary.csv",
                 secondary="evaluation_secondary.csv", fufa="evaluation_fufa.csv",
                 distribution="evaluation_distribution.csv")
    tables = {k: pd.read_csv(tables_dir / v) for k, v in names.items()}
    for k in ("metrics",):
        tables[k]["direction"] = tables[k].direction.astype(str)
    magnitude = post_hoc_magnitude(base, tables["fufa"])
    magnitude.to_csv(base / "post_hoc_downstream_magnitude.csv", index=False)
    write_summary(tables, run, md, {c: sem_lambda(c) for c in ("DMSO", "Vem")}, magnitude)
    style()
    figures.mkdir(parents=True, exist_ok=True)
    paired_z_figure(tables["metrics"], figures / "evaluation_paired_z.png")
    ceiling_figure(tables["metrics"], figures / "evaluation_r_vs_noise_ceiling.png")
    for p in (md, figures / "evaluation_paired_z.png", figures / "evaluation_r_vs_noise_ceiling.png"):
        print(f"wrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
