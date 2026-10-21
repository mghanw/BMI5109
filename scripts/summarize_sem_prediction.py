"""Aggregate summaries of the SEM knockdown predictions (no evaluation metrics).

Usage:
    python scripts/summarize_sem_prediction.py [--config config/sem_prediction.json]

Reads the merged predictions (outputs/sem_prediction/predictions/<condition>/fold<f>/) and the
structure tables, and writes results/summary/sem_prediction_model_summary.csv and
sem_prediction_summary.md. Only prediction-side quantities are reported (constant predictions,
genes with a nonzero predicted change, local IDA statistics, exclusions, runtime); no observed
response on evaluation cells, correlation or error is computed.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from aim1 import assert_no_cell_ids  # noqa: E402
from icpdag.sem_prediction import MODELS  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

SUMMARY = ROOT / "results" / "summary"
EXCLUSION_FLAG = 0.05


def load_predictions(cfg, conditions, folds):
    rows, stats, seconds = [], [], []
    base = ROOT / cfg["outputs_dir"] / "predictions"
    for condition in conditions:
        for fold in folds:
            d = base / condition / f"fold{fold}"
            with np.load(d / "predictions.npz") as z:
                models, targets, nodes = list(z["models"]), list(z["targets"]), list(z["nodes"])
                per_unit, included = z["per_unit"], z["included"]
                n_ref = len(z["reference_cells"])
            if models != MODELS:
                raise ValueError(f"{d}: unexpected model list")
            st = pd.read_csv(d / "prediction_stats.csv")
            for i, model in enumerate(models):
                for t, target in enumerate(targets):
                    node = st[(st.model == model) & (st.target == target)].node.iloc[0]
                    vector = per_unit[i, t].copy()
                    vector[nodes.index(node)] = 0.0
                    rows.append(dict(model=model, condition=condition, fold=fold, target=target,
                                     included_cells=int(included[i, t]), reference_cells=n_ref,
                                     genes_nonzero=int((np.abs(vector) > 1e-12).sum()),
                                     constant_prediction=bool(np.all(np.abs(vector) <= 1e-12))))
            stats.append(st)
            merge = json.loads((d / "merge.json").read_text(encoding="utf-8"))
            seconds.append(dict(condition=condition, fold=fold, **merge["seconds"]))
    return pd.DataFrame(rows), pd.concat(stats, ignore_index=True), pd.DataFrame(seconds)


def model_summary(pred, stats):
    st = stats.copy()
    st["valid_choices"] = st.n_candidates - st.n_failed
    st["excluded_share"] = st.excluded / st.units.where(st.units > 0)
    rows = []
    for (model, condition), p in pred.groupby(["model", "condition"], sort=False):
        s = st[(st.model == model) & (st.condition == condition)]
        flagged = sorted(s.target[s.excluded_share > EXCLUSION_FLAG])
        units = s.units.sum()
        rows.append(dict(
            model=model, condition=condition, n_targets=len(p),
            n_constant_prediction=int(p.constant_prediction.sum()),
            constant_targets=";".join(sorted(p.target[p.constant_prediction])),
            genes_nonzero_median=float(p.genes_nonzero.median()), genes_nonzero_min=int(p.genes_nonzero.min()),
            genes_nonzero_max=int(p.genes_nonzero.max()),
            reference_units=int(units), excluded_units=int(s.excluded.sum()),
            excluded_share=float(s.excluded.sum() / units) if units else np.nan,
            max_target_excluded_share=float(s.excluded_share.max()) if len(s) else np.nan,
            n_targets_exclusion_flag=len(flagged), exclusion_flag_targets=";".join(flagged),
            mean_undirected_per_unit=float(s.n_undirected.sum() / units) if units else np.nan,
            mean_valid_choices_per_unit=float(s.valid_choices.sum() / units) if units else np.nan,
            failed_choices=int(s.n_failed.sum()),
            fallback_share_of_valid_choices=(float(s.n_fallback.sum() / s.valid_choices.sum())
                                             if s.valid_choices.sum() else np.nan)))
    return pd.DataFrame(rows)


def either_condition(pred):
    wide = pred.pivot_table(index=["model", "target"], columns="condition", values="constant_prediction",
                            aggfunc="first")
    return wide.any(axis=1).groupby(level="model").sum().reindex(MODELS).astype(int)


def table(frame, columns, fmt=None):
    fmt = fmt or {}
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for _, r in frame.iterrows():
        lines.append("| " + " | ".join(fmt.get(c, "{}").format(r[c]) for c in columns) + " |")
    return lines


def markdown(cfg, summary, either, seconds, lam, shared, structure, orientation):
    md = ["# CSCN-SEM knockdown predictions", "",
          "Predicted pseudobulk expression changes (log1p CP10k) for every pre-registered model, held-out target,",
          "condition and cross-fit direction (config/sem_prediction.json, frozen before the first prediction run).",
          "This summary contains no evaluation metric: no observed response on evaluation cells, correlation or error",
          "was computed. The evaluation follows config/evaluation_endpoints.json including its 2026-10-21 amendment",
          "(constant predictions get r = 0 in the primary decision test; complete pairs are supporting).", "",
          "## Graphs (structure only)"]
    sh = shared[shared.graph.isin(["shared"])]
    for _, r in sh.iterrows():
        md.append(f"- {r.condition} shared graph: {int(r.edges)} edges, {int(r.targets_with_edge)}/56 targets with an "
                  f"incident edge (required {int(r.required_targets_with_edge)}), log10 alpha {r.log10_alpha:.1f}.")
    cl = shared[shared.graph.str.startswith("cluster_")]
    for c, g in cl.groupby("condition"):
        md.append(f"- {c} cluster graphs (K = 6): {int(g.edges.min())}-{int(g.edges.max())} edges, "
                  f"{int(g.targets_with_edge.min())}-{int(g.targets_with_edge.max())} targets with an edge.")
    md += ["- Shared and cluster graphs use a population Spearman partial-correlation test inside the vendored PC-stable",
           "  search (target-coverage matched); the cell-specific graphs use the local CSCN test. This difference in test",
           "  family is a stated limitation.", "",
           "## Orientation and local IDA (structure + training-fold seeds)"]
    agg = structure.groupby(["graph_type", "condition"]).agg(max_U=("max_U", "max"),
                                                              fallback=("share_fold_dag_fallback", "mean"),
                                                              no_choice=("share_units_no_valid_choice", "mean")).reset_index()
    md += table(agg, ["graph_type", "condition", "max_U", "fallback", "no_choice"],
                {"fallback": "{:.3f}", "no_choice": "{:.3f}"})
    md += ["",
           "- fallback: share of fold DAGs completed only under the relaxed rule (colliders may change); no_choice: share",
           "  of (reference cell, held-out target) units without any valid completion.",
           "- The remaining 3.5% of DMSO cell-specific units without a valid choice come from 70 of 1,989 DMSO cells whose",
           "  stored post-collider PDAG already contains a directed cycle among its fixed (collider-oriented) edges, from",
           "  conflicting collider orientations in the CSCN step. No completion can keep every fixed edge, so these cells",
           "  are excluded from the cell-specific and unweighted-ablation pseudobulk means. This is below the 5% exclusion",
           "  flag, so no rule was changed. Vem has no such cells.",
           f"- Seed proposals (eligible training-target responses, q < 0.05) per fold: "
           + "; ".join(f"{c} {', '.join(str(int(x)) for x in g.sort_values('fold').n_eligible_seeds.drop_duplicates())}"
                       for c, g in orientation.drop_duplicates(['condition', 'fold']).groupby('condition')) + ".", "",
           "## Shrinkage lambda (NTC-only cross-validation)"]
    for _, r in lam[lam.selected].iterrows():
        inf = lam[(lam.condition == r.condition) & (lam["lambda"].astype(str) == "inf")].mse.iloc[0]
        md.append(f"- {r.condition}: selected lambda {r['lambda']} (grid boundary: {bool(r.at_grid_boundary)}); "
                  f"MSE {r.mse:.5f} vs {inf:.5f} for lambda = infinity (reference only).")
    md += ["- The NTC-only cross-validation selected the boundary value 10 in both conditions, and pure pooled coefficients",
           "  (lambda = infinity, reference only, not selectable) fit marginally better: there is no detectable gain from",
           "  local coefficients in control cells. The cell-specific model therefore differs from the shared model mainly",
           "  through per-cell structure. Lambda and the grid were not changed.", "",
           "## Predictions per model"]
    s = summary.copy()
    md += table(s, ["model", "condition", "n_constant_prediction", "genes_nonzero_median", "excluded_share",
                    "n_targets_exclusion_flag", "mean_undirected_per_unit", "fallback_share_of_valid_choices"],
                {"genes_nonzero_median": "{:.0f}", "excluded_share": "{:.3f}", "mean_undirected_per_unit": "{:.3f}",
                 "fallback_share_of_valid_choices": "{:.3f}"})
    md += ["", "- n_constant_prediction: held-out targets whose predicted change is 0 for every gene other than the",
           "  perturbed gene (undefined Pearson r; r = 0 under the amended primary rule). Targets constant in at least one",
           "  condition (NaN under the supporting complete-pairs rule): "
           + ", ".join(f"{m} {int(n)}" for m, n in either.items()) + ".",
           "- zero_downstream_baseline is constant by definition (RMSE only).",
           f"- Exclusion flag: (model, target, condition) rows with > {EXCLUSION_FLAG:.0%} of reference cells excluded.", "",
           "## Runtime"]
    total = seconds.drop(columns=["condition", "fold"]).sum()
    md.append("- CPU time over all folds and conditions (14 workers; about 1.6 h wall for the full run): "
              + ", ".join(f"{m} {total.get(m, 0) / 3600:.2f} h" for m in MODELS if m in total) + ".")
    return "\n".join(md) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "sem_prediction.json"))
    args = parser.parse_args(argv)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    split = load_split(ROOT / cfg["inputs"]["target_split"])
    folds = sorted(split.fold.unique())
    pred, stats, seconds = load_predictions(cfg, cfg["conditions"], folds)
    if sorted(pred.target.unique()) != sorted(split.target):
        raise ValueError("Predictions do not cover every eligible target")
    summary = model_summary(pred, stats)
    either = either_condition(pred)
    for frame in (summary,):
        assert_no_cell_ids(frame)
    summary.to_csv(SUMMARY / "sem_prediction_model_summary.csv", index=False, lineterminator="\n")
    text = markdown(cfg, summary, either, seconds, pd.read_csv(SUMMARY / "sem_prediction_lambda_cv.csv"),
                    pd.read_csv(SUMMARY / "sem_prediction_shared_graphs.csv"),
                    pd.read_csv(SUMMARY / "sem_prediction_structure.csv"),
                    pd.read_csv(SUMMARY / "sem_prediction_orientation.csv"))
    (SUMMARY / "sem_prediction_summary.md").write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {SUMMARY / 'sem_prediction_model_summary.csv'} and sem_prediction_summary.md")
    print(summary[["model", "condition", "n_constant_prediction", "genes_nonzero_median", "excluded_share",
                   "n_targets_exclusion_flag"]].to_string(index=False))
    print("constant in at least one condition:", either.to_dict())


if __name__ == "__main__":
    main()
