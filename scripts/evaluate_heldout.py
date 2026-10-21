"""Held-out evaluation of the frozen CSCN-SEM knockdown predictions.

Usage:
    python scripts/evaluate_heldout.py                                   # the single full run
    python scripts/evaluate_heldout.py --fold 0 --targets A,B            # smoke test (outputs/evaluation_smoke/)
    python scripts/evaluate_heldout.py --deviation "YYYY-MM-DD: reason"  # documented rerun after a bug fix

Rules: config/evaluation_run.json (frozen), implementing config/evaluation_endpoints.json and its
2026-10-21 amendment. Order: (1) verify the predictions and compute every prediction-side check
before any observed response is read; (2) observed target - NTC changes on the evaluation cells
(calib_fold 1 - c); (3) per-target metrics, noise ceilings, primary and supporting tests; (4) secondary
comparisons and the Fufa secondary endpoint; (5) tertiary distribution endpoints.

The full run writes outputs/evaluation/ (never committed) and results/summary/evaluation_*.csv; then
run scripts/summarize_evaluation.py. Printed output: verification counts, shapes, row counts and runtime
only; no metric value is printed.
"""
import argparse
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from aim1 import assert_no_cell_ids  # noqa: E402
from icpdag import dispersion as ds  # noqa: E402
from icpdag import heldout_evaluation as ev  # noqa: E402
from icpdag.sem_inputs import load_panel_cells, panel_nodes, target_nodes  # noqa: E402
from icpdag.sem_prediction import MODELS, fufa_weights, target_delta  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

RUN_CONFIG = ROOT / "config" / "evaluation_run.json"
GUARDED = ["config/evaluation_run.json", "icpdag/heldout_evaluation.py", "scripts/evaluate_heldout.py"]
FUFA = "Fufa2019_AZD_up_SOX10_bound_pigmentation"
PRIMARY, REFERENCE = "cell_specific_cscn_sem", "shared_cscn_sem"
R_MODELS = [m for m in MODELS if m != "zero_downstream_baseline"]
TABLES = ["evaluation_target_metrics.csv", "evaluation_primary.csv", "evaluation_secondary.csv",
          "evaluation_fufa.csv", "evaluation_distribution.csv"]


def say(message):
    print(message, flush=True)


# ---------------------------------------------------------------- configuration and guard
def load_configs():
    run = json.loads(RUN_CONFIG.read_text(encoding="utf-8"))
    if not run.get("frozen"):
        raise ValueError("config/evaluation_run.json must be frozen")
    inputs = run["inputs"]
    endpoints = json.loads((ROOT / inputs["endpoints"]).read_text(encoding="utf-8"))
    if ev.file_sha256(ROOT / inputs["endpoints"]) != inputs["endpoints_sha256"]:
        raise ValueError("config/evaluation_endpoints.json differs from the frozen hash")
    committed = ev.git_file_sha256(ROOT, inputs["prediction_config_commit"], inputs["prediction_config"])
    if committed != inputs["prediction_config_sha256"] or \
            ev.file_sha256(ROOT / inputs["prediction_config"]) != inputs["prediction_config_sha256"]:
        raise ValueError("config/sem_prediction.json differs from the version committed in 81a3bd1")
    sem = json.loads((ROOT / inputs["prediction_config"]).read_text(encoding="utf-8"))
    if run["models"] != MODELS or run["conditions"] != endpoints["conditions"]:
        raise ValueError("Model or condition lists differ from the frozen definitions")
    if run["clip_abs_r"] != endpoints["primary"]["averaging"]["clip_abs_r"] or \
            run["bootstrap_draws"] != endpoints["primary"]["effect_size"]["ci"]["draws"] or \
            run["seed"] != endpoints["primary"]["effect_size"]["ci"]["seed"] or \
            run["noise_ceiling_splits"] != endpoints["primary"]["noise_ceiling"]["splits"]:
        raise ValueError("Run constants differ from config/evaluation_endpoints.json")
    return run, endpoints, sem


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def run_guard(out, deviation):
    tracked = set(git("ls-files", *GUARDED).split())
    if tracked != set(GUARDED):
        raise RuntimeError(f"Not committed: {sorted(set(GUARDED) - tracked)}")
    if git("status", "--porcelain", "--", *GUARDED).strip():
        raise RuntimeError("Evaluation code or run rules have uncommitted changes")
    history = []
    previous, started = out / "run.json", out / "started.json"
    if (previous.exists() or started.exists()) and not deviation:
        raise RuntimeError("outputs/evaluation/ already holds a started or finished run: the evaluation runs once; "
                           "a rerun needs --deviation 'YYYY-MM-DD: reason'")
    if previous.exists():
        old = json.loads(previous.read_text(encoding="utf-8"))
        history = old.pop("history", []) + [old]
    elif started.exists():
        history = [dict(json.loads(started.read_text(encoding="utf-8")), finished=False)]
    head = git("rev-parse", "HEAD").strip()
    out.mkdir(parents=True, exist_ok=True)
    started.write_text(json.dumps(dict(date=dt.datetime.now().isoformat(timespec="seconds"), head=head,
                                       deviation=deviation), indent=2), encoding="utf-8")
    return dict(head=head, history=history, deviation=deviation)


# ---------------------------------------------------------------- prediction-side verification
def verify(run, sem, split, nodes, tnodes):
    inputs = run["inputs"]
    preds = ev.verify_predictions(ROOT / "outputs" / "sem_prediction" / "predictions",
                                  inputs["prediction_config_sha256"], split, nodes, MODELS, run["conditions"])
    zero = ev.included_zero_rows(preds, MODELS)
    if len(zero):
        raise RuntimeError("Prediction rows with included == 0 (stop and report):\n" + zero.to_string(index=False))
    flags = ev.constant_prediction_flags(preds, MODELS, nodes, tnodes)
    structure = pd.read_csv(ROOT / inputs["structure_targets"])
    counts = ev.check_constant_flags(flags, structure)
    ev.check_constant_shift(preds, MODELS, run["rules"]["tertiary"]["mean_shift_models"])
    spec = json.loads((ROOT / inputs["fufa_weights"]).read_text(encoding="utf-8"))
    if spec["program"] != FUFA:
        raise ValueError("Frozen Fufa weights are for another program")
    w, _ = fufa_weights(spec["program_genes"], spec["controls"], spec["scale"], nodes)
    ev.check_fufa_pseudobulk(preds, w)
    return preds, flags, counts


# ---------------------------------------------------------------- observed-side inputs
def load_observed_inputs(run, sem, split, ids, nodes, tnodes, preds):
    inputs = run["inputs"]
    values, meta = load_panel_cells(ROOT / inputs["expression"], ids, run["conditions"], split.target)
    meta["calib_fold"] = meta.calib_fold.astype(int)
    cells = pd.read_csv(ROOT / inputs["cells"], usecols=["cell", "calib_fold"]).set_index("cell")
    if not np.array_equal(cells.calib_fold.reindex(meta.index).to_numpy(), meta.calib_fold.to_numpy()):
        raise ValueError("calib_fold differs between the expression file and data/processed/cells.csv.gz")
    frame = pd.DataFrame(values, index=meta.index, columns=nodes)
    ntc_ids = {}
    for condition in run["conditions"]:
        graph_cells = pd.read_csv(ROOT / inputs["cell_graphs"].replace("<condition>", condition), usecols=["cell_id"])
        ntc = meta.index[(meta.condition.eq(condition) & meta.target.eq("NTC")).to_numpy()]
        if list(graph_cells.cell_id) != list(ntc):
            raise ValueError(f"{condition}: NTC rows differ from the stored graph cells")
        p = preds[condition]
        ntc_ids[condition] = ntc[p["reference_cells"]]
        for t, target in enumerate(p["targets"]):
            delta = target_delta(frame, meta, target, tnodes[target], condition)
            if not np.allclose([delta[0], delta[1]], p["delta"][t], rtol=0, atol=1e-6):
                raise ValueError(f"{condition} {target}: recomputed delta differs from the stored delta")
    return values, meta, ntc_ids


def load_fufa_scores(run, meta):
    inputs = run["inputs"]
    stored = json.loads((ROOT / "outputs" / "aim1" / "config.json").read_text(encoding="utf-8"))
    frozen = json.loads((ROOT / inputs["aim1_config"]).read_text(encoding="utf-8"))
    run_fields = {"argv", "run_name", "targets_subset"}
    a = {k: v for k, v in stored.items() if k not in run_fields}
    b = {k: v for k, v in frozen.items() if k not in run_fields}
    a["paths"] = dict(a["paths"], target_qc=b["paths"]["target_qc"]) \
        if a["paths"].get("target_qc") == "results/summary/phase1_target_qc.csv" else a["paths"]
    if a != b or stored.get("targets_subset"):
        raise ValueError("outputs/aim1 is not the frozen full Aim 1 run")
    scores = pd.read_csv(ROOT / inputs["aim1_scores"], index_col=0)
    scores = scores[["condition", "target", FUFA]]
    groups = pd.read_csv(ROOT / inputs["aim1_group_summaries"])
    groups = groups[groups.readout.eq(FUFA)].set_index(["condition", "target"])
    mine = scores.groupby(["condition", "target"])[FUFA].agg(["size", "mean"]).reindex(groups.index)
    if not (np.array_equal(mine["size"].to_numpy(), groups.n_cells.to_numpy())
            and np.allclose(mine["mean"].to_numpy(), groups["mean"].to_numpy(), rtol=1e-5, atol=1e-6)):
        raise ValueError("Aim 1 per-cell Fufa scores do not reproduce aim1_group_summaries.csv")
    fufa = scores[FUFA].reindex(meta.index)
    if fufa.isna().any():
        raise ValueError("Cells without an Aim 1 Fufa score")
    return fufa.to_numpy(float)


# ---------------------------------------------------------------- analyses
def table_values(table, model, column, targets, condition="pooled"):
    sub = table[table.model.eq(model) & table.condition.eq(condition) & table.direction.eq("combined")]
    return sub.set_index("target")[column].reindex(targets).to_numpy(float)


def test_row(section, analysis, role, version, metric, model, reference, condition, a, b, run, back):
    out = dict(section=section, analysis=analysis, role=role, version=version, metric=metric, model=model,
               reference=reference, condition=condition)
    out.update(ev.paired_test(a, b, run["bootstrap_draws"], run["seed"], back_transform=back))
    return out


def median_rows(table, targets, models, condition, extra):
    rows = []
    for m in models:
        if m in R_MODELS:
            for version, column in (("nan_excluded", "z"), ("r_zero", "z_amended")):
                z = table_values(table, m, column, targets, condition)
                z = z[np.isfinite(z)]
                med = float(np.median(z)) if len(z) else np.nan
                rows.append(dict(extra, model=m, condition=condition, version=version, metric="r",
                                 n_targets=len(z), median_z=med, median_r=float(np.tanh(med))))
        e = table_values(table, m, "rmse", targets, condition)
        rows.append(dict(extra, model=m, condition=condition, version="all", metric="rmse", n_targets=len(e),
                         median_rmse=float(np.median(e))))
    return rows


def primary_tables(table, targets, split, stratum, run, endpoints):
    rows = []
    a = lambda m, col, cond="pooled": table_values(table, m, col, targets, cond)  # noqa: E731
    rows.append(test_row("test", "primary_decision", "primary", "amended", "z", PRIMARY, REFERENCE, "pooled",
                         a(PRIMARY, "z_amended"), a(REFERENCE, "z_amended"), run, True))
    rows.append(test_row("test", "complete_pairs", "supporting", "complete_pairs", "z", PRIMARY, REFERENCE, "pooled",
                         a(PRIMARY, "z"), a(REFERENCE, "z"), run, True))
    for condition in run["conditions"]:
        rows.append(test_row("test", f"per_condition_{condition}", "supporting", "amended", "z", PRIMARY, REFERENCE,
                             condition, a(PRIMARY, "z_amended", condition), a(REFERENCE, "z_amended", condition),
                             run, True))
    rows.append(test_row("test", "rmse", "supporting", "all", "rmse", PRIMARY, REFERENCE, "pooled",
                         a(PRIMARY, "rmse"), a(REFERENCE, "rmse"), run, False))
    alpha = endpoints["primary"]["decision_test"]["alpha"]
    rows[0]["alpha"] = alpha
    rows[0]["reject_h0"] = bool(np.isfinite(rows[0]["p_value"]) and rows[0]["p_value"] < alpha)
    threshold = endpoints["primary"]["nan_rule"]["flag_threshold"]
    for condition in [*run["conditions"], "pooled"]:
        rows += median_rows(table, targets, MODELS, condition, dict(section="median", role="supporting"))
    fold = split.set_index("target").fold
    for f in sorted(fold.unique()):
        members = [t for t in targets if fold[t] == f]
        if members:
            rows += median_rows(table, members, MODELS, "pooled", dict(section="fold_median", role="supporting",
                                                                         fold=int(f)))
    for name in ("connected", "weak"):
        members = [t for t in targets if stratum[t] == name]
        if members:
            rows += median_rows(table, members, MODELS, "pooled", dict(section="stratum_median",
                                                                         role="descriptive", stratum=name))
    flag_any = False
    for m in MODELS:
        for condition in [*run["conditions"], "pooled"]:
            const = table[table.model.eq(m) & table.condition.eq(condition) & table.direction.eq("combined")]
            const = const.set_index("target").constant_prediction.reindex(targets).astype(bool)
            row = dict(section="coverage", role="descriptive", model=m,
                       condition="both" if condition == "pooled" else condition, n_targets=len(targets),
                       n_nonconstant=int((~const).sum()), share_nonconstant=float((~const).mean()))
            if condition == "pooled" and m in R_MODELS:
                z = a(m, "z")
                row.update(n_nan=int(np.isnan(z).sum()), nan_share=float(np.isnan(z).mean()))
                row["nan_flag"] = bool(row["nan_share"] > threshold)
                flag_any |= row["nan_flag"]
            rows.append(row)
    rows.append(dict(section="nan_flag", role="supporting", n_targets=len(targets), nan_flag=bool(flag_any)))
    return pd.DataFrame(rows)


def secondary_table(table, targets, run):
    rows = []
    a = lambda m, col: table_values(table, m, col, targets)  # noqa: E731
    families = {"secondary_z_amended": ("amended", "z_amended", "descriptive"),
                "secondary_z_complete_pairs": ("complete_pairs", "z", "supporting")}
    for family, (version, column, role) in families.items():
        for m in run["rules"]["holm"]["families"]["secondary_z_amended"]:
            rows.append(dict(test_row("test", f"{m}_z_{version}", role, version, "z", PRIMARY, m, "pooled",
                                      a(PRIMARY, column), a(m, column), run, True), holm_family=family))
    for m in run["rules"]["holm"]["families"]["secondary_rmse"]:
        rows.append(dict(test_row("test", f"{m}_rmse", "descriptive", "all", "rmse", PRIMARY, m, "pooled",
                                  a(PRIMARY, "rmse"), a(m, "rmse"), run, False), holm_family="secondary_rmse"))
    frame = pd.DataFrame(rows)
    frame["p_holm"] = np.nan
    for _, idx in frame.groupby("holm_family").groups.items():
        frame.loc[idx, "p_holm"] = ev.holm(frame.loc[idx, "p_value"].to_numpy(float))
    return frame


def fufa_table(preds, meta, fufa, targets, split, run):
    rows = []
    obs = {}
    for condition in run["conditions"]:
        for target in targets:
            obs[condition, target] = np.array([ev.observed_mean_contrast(fufa, meta, target, condition, c)
                                               for c in (0, 1)])
    fold = split.set_index("target").fold
    for i, m in enumerate(MODELS):
        for target in targets:
            predicted = {}
            for condition in run["conditions"]:
                p = preds[condition]
                t = p["targets"].index(target)
                predicted[condition] = np.nanmean(p["fufa_scaled"][i, :, t, :], axis=0)
            for row in ev.fufa_contrast_rows(predicted, {c: obs[c, target] for c in run["conditions"]}):
                rows.append(dict(section="target", model=m, target=target, fold=int(fold[target]), **row))
    frame = pd.DataFrame(rows)
    combined = frame[frame.direction.eq("combined")]
    summary = []
    for contrast in combined.contrast.unique():
        sub = combined[combined.contrast.eq(contrast)]
        for m in MODELS:
            e = sub[sub.model.eq(m)].set_index("target").abs_error.reindex(targets)
            summary.append(dict(section="median", model=m, contrast=contrast, n_targets=len(e),
                                median_abs_error=float(e.median())))
        a = sub[sub.model.eq(PRIMARY)].set_index("target").abs_error.reindex(targets).to_numpy(float)
        b = sub[sub.model.eq(REFERENCE)].set_index("target").abs_error.reindex(targets).to_numpy(float)
        summary.append(dict(section="test", model=PRIMARY, reference=REFERENCE, contrast=contrast, metric="abs_error",
                            **ev.paired_test(a, b, run["bootstrap_draws"], run["seed"])))
    return pd.concat([frame, pd.DataFrame(summary)], ignore_index=True)


def distribution_table(preds, meta, fufa, ntc_ids, targets, split, run):
    seed, draws = run["seed"], run["bootstrap_draws"]
    position = pd.Series(np.arange(len(meta)), index=meta.index)
    fold = split.set_index("target").fold
    aim1 = pd.read_csv(ROOT / run["inputs"]["aim1_dispersion"])
    aim1 = aim1[aim1.statistic.eq("log_var_ratio") & aim1.kind.eq("gene") & aim1.readout.eq(FUFA)]
    aim1 = aim1.set_index(["target", "condition"])[["estimate", "ci_low", "ci_high"]]
    rows, dropped = [], []
    for condition in run["conditions"]:
        p = preds[condition]
        ref_rows = position[ntc_ids[condition]].to_numpy()
        ref_scores, ref_folds = fufa[ref_rows], meta.calib_fold.to_numpy()[ref_rows]
        for c in (0, 1):
            e = 1 - c
            eval_pos = np.flatnonzero(ref_folds == e)
            y = ref_scores[eval_pos]
            ntc_draws = ds.GroupDraws(y, draws, seed, condition, f"NTC|eval{e}")
            for target in targets:
                t = p["targets"].index(target)
                x = fufa[ev.evaluation_rows(meta, target, condition, c)]
                observed = ev.distribution_stats(x, y)
                low, high = ev.observed_lvr_interval(x, y, draws, seed, condition, target, e, ntc_draws=ntc_draws)
                ref = aim1.loc[(target, condition)]
                for i, m in enumerate(MODELS):
                    predicted, n_drop = ev.predicted_distribution_stats(y, p["fufa_scaled"][i, eval_pos, t, c])
                    dropped.append(dict(model=m, condition=condition, direction=c, target=target, n_dropped=n_drop))
                    for statistic, value in predicted.items():
                        row = dict(section="target", model=m, target=target, fold=int(fold[target]),
                                   condition=condition, direction=c, statistic=statistic, predicted=value,
                                   observed=observed[statistic], abs_error=abs(value - observed[statistic]),
                                   n_ntc_eval=len(y), n_included=len(y) - n_drop, n_dropped=n_drop, n_target=len(x))
                        if statistic == "log_var_ratio":
                            row.update(observed_ci_low=low, observed_ci_high=high,
                                       **{k: float(v) for k, v in ev.calibration_flags(value, low, high).items()},
                                       aim1_all_cell_estimate=ref.estimate, aim1_all_cell_ci_low=ref.ci_low,
                                       aim1_all_cell_ci_high=ref.ci_high)
                        rows.append(row)
    frame = pd.DataFrame(rows)
    keys = ["model", "target", "fold", "condition", "statistic"]
    agg = {k: "mean" for k in ("predicted", "observed", "abs_error", "over_dispersion", "covered",
                               "observed_ci_low", "observed_ci_high")}
    agg.update(n_dropped="sum", aim1_all_cell_estimate="first", aim1_all_cell_ci_low="first",
               aim1_all_cell_ci_high="first")
    combined = frame.groupby(keys, sort=False).agg(agg).reset_index()
    combined["section"], combined["direction"] = "target", "combined"
    out = [frame, combined]
    summary = []
    for m in MODELS:
        for condition in [*run["conditions"], "both"]:
            mask = combined.model.eq(m)
            if condition != "both":
                mask &= combined.condition.eq(condition)
            sub = combined[mask]
            for statistic in ("log_var_ratio", "iqr_diff", "centered_wasserstein"):
                s = sub[sub.statistic.eq(statistic)]
                summary.append(dict(section="median", model=m, condition=condition, statistic=statistic,
                                    n_targets=len(s), median_abs_error=float(s.abs_error.median())))
            s = sub[sub.statistic.eq("log_var_ratio")]
            summary.append(dict(section="calibration", model=m, condition=condition, statistic="log_var_ratio",
                                n_rows=len(s), over_dispersion_rate=float(s.over_dispersion.mean()),
                                coverage=float(s.covered.mean()), nominal_over_dispersion=0.025))
    drop = pd.DataFrame(dropped)
    for (m, condition, c), s in drop.groupby(["model", "condition", "direction"], sort=False):
        summary.append(dict(section="dropped_cells", model=m, condition=condition, direction=c,
                            n_targets=len(s), n_dropped=int(s.n_dropped.sum()), max_dropped=int(s.n_dropped.max())))
    out.append(pd.DataFrame(summary))
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------- main
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fold", type=int, help="Smoke test: outer fold of the --targets")
    parser.add_argument("--targets", help="Smoke test: comma-separated held-out targets of --fold")
    parser.add_argument("--deviation", help="Documented rerun after a bug fix: 'YYYY-MM-DD: reason'")
    args = parser.parse_args(argv)
    smoke = bool(args.targets)
    t_start = time.time()
    run, endpoints, sem = load_configs()
    out = ROOT / "outputs" / ("evaluation_smoke" if smoke else "evaluation")
    guard = None if smoke else run_guard(out, args.deviation)
    split = load_split(ROOT / sem["inputs"]["target_split"])
    panel = json.loads((ROOT / sem["inputs"]["panel"]).read_text(encoding="utf-8"))
    ids, nodes = panel_nodes(panel)
    tnodes = target_nodes(panel)

    preds, flags, counts = verify(run, sem, split, nodes, tnodes)
    say("verification passed: hashes, coverage, delta-scaled shapes, Fufa pseudobulk, constant shifts "
        f"({', '.join(run['rules']['tertiary']['mean_shift_models'])})")
    say("targets with a constant prediction (per-unit, 249 genes):\n" + counts.to_string())
    say("rows with included == 0: 0")

    targets = sorted(split.target)
    if smoke:
        chosen = sorted(args.targets.split(","))
        held = set(split.target[split.fold.eq(args.fold)])
        if args.fold is None or set(chosen) - held:
            raise ValueError("Smoke targets must be held out in --fold")
        targets = chosen
    out.mkdir(parents=True, exist_ok=True)
    tables_dir = out / "tables" if smoke else ROOT / "results" / "summary"
    tables_dir.mkdir(parents=True, exist_ok=True)

    values, meta, ntc_ids = load_observed_inputs(run, sem, split, ids, nodes, tnodes, preds)
    say(f"evaluation input: {values.shape[0]} cells x {values.shape[1]} panel genes; delta recomputation matches")
    observed = {}
    for condition in run["conditions"]:
        for target in targets:
            observed[condition, target] = np.stack([ev.observed_change(values, meta, target, condition, c)
                                                    for c in (0, 1)])
    np.savez_compressed(out / "observed_changes.npz", conditions=np.asarray(run["conditions"]),
                        targets=np.asarray(targets), nodes=np.asarray(nodes),
                        observed=np.stack([np.stack([observed[c, t] for t in targets]) for c in run["conditions"]]))
    clip = run["clip_abs_r"]
    rows = []
    for condition in run["conditions"]:
        p = preds[condition]
        for target in targets:
            t = p["targets"].index(target)
            mask = ev.evaluation_mask(nodes, tnodes[target])
            for i, m in enumerate(MODELS):
                for c in (0, 1):
                    rows.append(dict(model=m, target=target, condition=condition, direction=c,
                                     **ev.direction_metrics(p["scaled"][i, t, c, mask], observed[condition, target][c, mask],
                                                            clip)))
    table = ev.combine_metrics(pd.DataFrame(rows))

    splits, ceiling = [], []
    for condition in run["conditions"]:
        for target in targets:
            mask = ev.evaluation_mask(nodes, tnodes[target])
            per_direction = []
            for c in (0, 1):
                rng = ev.noise_ceiling_rng(run["seed"], condition, target, c)
                tv = values[ev.evaluation_rows(meta, target, condition, c)][:, mask]
                nv = values[ev.evaluation_rows(meta, "NTC", condition, c)][:, mask]
                v = ev.split_half_ceiling(tv, nv, run["noise_ceiling_splits"], rng, clip)
                splits += [dict(condition=condition, target=target, direction=c, split=s, value=x)
                           for s, x in enumerate(v)]
                per_direction.append(float(np.median(v)))
                ceiling.append(dict(condition=condition, target=target, direction=c, noise_ceiling_r=per_direction[-1]))
            ceiling.append(dict(condition=condition, target=target, direction="combined",
                                noise_ceiling_r=ev.combine_directions_z(per_direction, clip)))
    pd.DataFrame(splits).to_csv(out / "noise_ceiling_splits.csv.gz", index=False)
    ceiling = pd.DataFrame(ceiling)
    ceiling["direction"] = ceiling.direction.astype(object)
    table["direction"] = table.direction.astype(object)
    table = table.merge(ceiling, on=["condition", "target", "direction"], how="left")
    connected = set(sem["weak_targets"]["connected_targets"])
    stratum = {t: "connected" if t in connected else "weak" for t in targets}
    fold = split.set_index("target").fold
    table["fold"] = table.target.map(fold).astype(int)
    table["stratum"] = table.target.map(stratum)
    table = table[["model", "target", "fold", "stratum", "condition", "direction", "r", "z", "r_amended", "z_amended",
                   "rmse", "constant_prediction", "constant_observed", "noise_ceiling_r"]]

    fufa = load_fufa_scores(run, meta)
    frames = {"evaluation_target_metrics.csv": table,
              "evaluation_primary.csv": primary_tables(table, targets, split, stratum, run, endpoints),
              "evaluation_secondary.csv": secondary_table(table, targets, run),
              "evaluation_fufa.csv": fufa_table(preds, meta, fufa, targets, split, run),
              "evaluation_distribution.csv": distribution_table(preds, meta, fufa, ntc_ids, targets, split, run)}
    for name, frame in frames.items():
        assert_no_cell_ids(frame)
        frame.to_csv(tables_dir / name, index=False)
        say(f"wrote {tables_dir.relative_to(ROOT) / name}: {frame.shape[0]} rows x {frame.shape[1]} columns")
    say(f"observed changes: {len(run['conditions'])} conditions x {len(targets)} targets x 2 directions x "
        f"{len(nodes)} genes; noise-ceiling splits: {len(splits)} rows")

    (out / "config.json").write_text(json.dumps(dict(run_config=run, argv=sys.argv[1:] if argv is None else argv,
                                                     smoke=smoke), indent=2), encoding="utf-8")
    record = dict(date=dt.datetime.now().isoformat(timespec="seconds"), seconds=round(time.time() - t_start, 1),
                  smoke=smoke, targets=len(targets), tables=[str((tables_dir / n).relative_to(ROOT)) for n in TABLES],
                  constant_counts=counts.reset_index().to_dict("records"))
    if guard:
        record.update(head=guard["head"], deviation=guard["deviation"], history=guard["history"])
    (out / "run.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
    say(f"done in {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
