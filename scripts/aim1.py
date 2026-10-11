"""Aim 1: program scores and descriptive mean / distribution contrasts.

Usage:
    python scripts/aim1.py [--config config/aim1.json] [--targets A,B,...] [--bootstrap B]
                           [--permutations P] [--run-name NAME]

Scores every primary-rule singlet with icpdag.program_scores (frozen DMSO NTC scale), then for the
eligible targets (config/target_split.json) and both conditions computes: mean contrasts
(icpdag.average_effects), dispersion and distribution contrasts, a log-UMI depth sensitivity, a
ligation-1 well cluster bootstrap (primary readout), pseudobulk expression contrasts, highlight
calls and a check against the earlier report. Every per-target table carries the target's fold.
With --targets (smoke test) every output goes to outputs/<run_name>/ only.
"""
import argparse
import json
import logging
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag import dispersion as ds  # noqa: E402
from icpdag import gene_sets as gs  # noqa: E402
from icpdag.average_effects import define_contrasts, estimate_contrasts  # noqa: E402
from icpdag.expression import detected_genes, expression_contrasts, summarize_expression  # noqa: E402
from icpdag.program_scores import program_scores  # noqa: E402
from icpdag.target_split import add_fold, load_split  # noqa: E402

LOG = logging.getLogger("aim1")
NTC = "NTC"


def assert_no_cell_ids(frame):
    """Aggregate outputs must not contain barcode sequences or cell names."""
    pattern = re.compile(r"[ACGT]{8,}")
    for col in frame.columns:
        values = frame[col].dropna().astype(str)
        if values.map(lambda v: bool(pattern.search(v))).any():
            raise ValueError(f"Column {col!r} contains barcode-like strings")


class Writer:
    def __init__(self, directory, split):
        self.directory, self.split = Path(directory), split
        self.directory.mkdir(parents=True, exist_ok=True)

    def __call__(self, frame, name, per_target=True):
        if per_target:
            frame = add_fold(frame, self.split)
            eligible = frame.target.isin(self.split.target)
            if frame.loc[eligible, "fold"].isna().any():
                raise ValueError(f"{name}: eligible target without fold")
        assert_no_cell_ids(frame)
        path = self.directory / name
        frame.to_csv(path, index=False)
        LOG.info("wrote %s (%d rows)", path, len(frame))
        return frame


# ---------------------------------------------------------------- loading and scoring
def load_data(cfg):
    import anndata as ad

    adata = ad.read_h5ad(ROOT / cfg["paths"]["rna"])
    obs = adata.obs.copy()
    for col in ("condition", "target", "guide", "lig1_well", "target_gene_id"):
        obs[col] = obs[col].astype(str)
    counts = sp.csr_matrix(adata.X)
    total = obs.total_counts.to_numpy(float)
    if not np.allclose(np.asarray(counts.sum(axis=1)).ravel(), total):
        raise ValueError("total_counts does not match the matrix row sums")
    return counts, adata.var.copy(), obs


def log_cp10k(counts, total):
    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    X.data = np.log1p(X.data / np.repeat(total, np.diff(X.indptr)).astype(np.float32) * np.float32(1e4))
    return X


def load_programs(path):
    with open(path, encoding="utf-8") as f:
        spec = json.load(f)
    mapping = pd.DataFrame([dict(program=name, readout_role=p["readout_role"], **g)
                            for name, p in spec["programs"].items() for g in p["genes"]])
    return spec, mapping


def library_genes(cfg, var, eligible):
    """Library target symbols (excluded from background controls) and each eligible target's gene."""
    qc = pd.read_csv(ROOT / cfg["paths"]["target_qc"])
    gene_symbol = pd.Series(var.symbol.to_numpy(), index=var.index.str.replace(r"\.\d+$", "", regex=True))
    library_ids = qc.target_gene_id.dropna().map(gs.unversioned).unique()
    library_symbols = set(gene_symbol.reindex(library_ids).dropna())
    target_gene = (qc.dropna(subset=["target_gene_id"]).drop_duplicates("target")
                   .set_index("target").target_gene_id.to_dict())
    target_symbol = {t: gene_symbol.get(gs.unversioned(target_gene[t])) for t in eligible}
    return library_symbols, target_gene, target_symbol


def score_cells(X, var, meta, programs, library_symbols, seed):
    annotation = var[["symbol", "gene_type"]]
    cfg = {"network_genes": [], "gene_readouts": sorted(library_symbols), "seed": seed}
    return program_scores(X, var.symbol.to_numpy(), annotation, meta, programs, cfg)


# ---------------------------------------------------------------- analyses
def mean_contrasts(scores, meta, targets, pairs, cfg):
    contrasts = define_contrasts(targets, pairs)
    effects, groups, _ = estimate_contrasts(scores, meta, contrasts, cfg["bootstrap"], cfg["seed"])
    effects["p_value"] = 2 * stats.norm.sf(np.abs(effects.estimate / effects.iid_cell_se))
    effects["ci_excludes_0"] = (effects.ci_low > 0) | (effects.ci_high < 0)
    effects["q_value"] = np.nan
    for _, idx in effects.groupby(["kind", "condition", "outcome"]).groups.items():
        effects.loc[idx, "q_value"] = ds.bh(effects.loc[idx, "p_value"])
    return effects.rename(columns={"outcome": "readout"}), groups


def group_summaries(scores, meta):
    rows = []
    for (condition, target), idx in meta.groupby(["condition", "target"]).groups.items():
        for readout in scores:
            v = scores.loc[idx, readout].to_numpy()
            q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
            rows.append(dict(condition=condition, target=target, readout=readout, n_cells=len(v), mean=v.mean(),
                             sd=v.std(ddof=1), q05=q[0], q25=q[1], median=q[2], q75=q[3], q95=q[4]))
    return pd.DataFrame(rows)


def depth_sensitivity(scores, meta, targets, primary, cfg, progress):
    log_umi = np.log(meta.total_counts.to_numpy(float))
    fit = (meta.condition.eq("DMSO") & meta.target.eq(NTC)).to_numpy()
    adjusted, fits = {}, []
    for readout in scores:
        adjusted[readout], slope, r2 = ds.residualize_depth(scores[readout], log_umi, fit)
        fits.append(dict(readout=readout, depth_slope_per_log_umi=slope, depth_r2_dmso_ntc=r2))
    adjusted = pd.DataFrame(adjusted, index=scores.index)
    sens = ds.dispersion_contrasts(adjusted, meta, targets, cfg["bootstrap"], cfg["permutations"], cfg["seed"],
                                   tuple(cfg["conditions"]), progress)
    keys = ["statistic", "kind", "target", "condition", "readout"]
    merged = sens.merge(primary[keys + ["estimate", "q_value"]], on=keys, suffixes=("_depth", "_primary"))
    merged = merged.merge(pd.DataFrame(fits), on="readout")
    flips = merged.q_value_primary.lt(0.05) != merged.q_value_depth.lt(0.05)
    merged["significance_flips"] = flips.astype(object).where(merged.q_value_primary.notna(), np.nan)
    depth = pd.DataFrame({"condition": meta.condition, "target": meta.target, "log_umi": log_umi})
    med = depth.groupby(["condition", "target"]).log_umi.median()
    merged["median_log_umi_diff_vs_ntc"] = [
        med.get((c, t), np.nan) - med.get((c, NTC), np.nan) if k == "gene" else np.nan
        for c, t, k in zip(merged.condition, merged.target, merged.kind)]
    return merged, adjusted


def cluster_table(scores, meta, targets, mean_fx, disp, cfg):
    c = cfg["cluster_bootstrap"]
    readout = c["readout"]
    m = meta.rename(columns={c["well_column"]: "well"})[["condition", "target", "well"]]
    draws = ds.cluster_bootstrap(scores[readout].to_numpy(), m, targets, c["repeats"], cfg["seed"],
                                 tuple(cfg["conditions"]))
    for condition, d in draws.items():
        LOG.info("cluster bootstrap %s: %d wells", condition, d["wells"])
    out = ds.cluster_contrasts(draws, targets, readout, tuple(cfg["conditions"]))
    cell_mean = mean_fx[mean_fx.kind.isin(["gene", "interaction"]) & mean_fx.readout.eq(readout)].assign(
        statistic="mean")[["statistic", "kind", "target", "condition", "estimate", "ci_low", "ci_high"]]
    cell_lvr = disp[disp.statistic.eq("log_var_ratio") & disp.kind.eq("gene") & disp.readout.eq(readout)][
        ["statistic", "kind", "target", "condition", "estimate", "ci_low", "ci_high"]]
    cell = pd.concat([cell_mean, cell_lvr]).rename(columns={"ci_low": "cell_ci_low", "ci_high": "cell_ci_high"})
    out = out.merge(cell, on=["statistic", "kind", "target", "condition"], how="left")
    out["cell_ci_excludes_0"] = (out.cell_ci_low > 0) | (out.cell_ci_high < 0)
    out["cluster_ci_excludes_0"] = (out.cluster_ci_low > 0) | (out.cluster_ci_high < 0)
    out["n_wells"] = [draws[c]["wells"] if c in draws else draws[cfg["conditions"][0]]["wells"] for c in out.condition]
    cols = ["statistic", "kind", "target", "condition", "readout", "estimate", "cell_ci_low", "cell_ci_high",
            "cluster_ci_low", "cluster_ci_high", "cell_ci_excludes_0", "cluster_ci_excludes_0",
            "cluster_draws_dropped", "n_wells"]
    return out[cols]


def highlight_table(mean_fx, disp, depth, cluster, targets, cfg):
    h = cfg["highlight"]
    q = h["q"]
    gene = mean_fx[mean_fx.kind.eq("gene")][["target", "condition", "readout", "estimate", "ci_low", "ci_high",
                                             "q_value"]]
    gene = gene.rename(columns={"estimate": "mean_estimate", "ci_low": "mean_ci_low", "ci_high": "mean_ci_high",
                                "q_value": "mean_q"})
    gene["mean_ci_includes_0"] = (gene.mean_ci_low <= 0) & (gene.mean_ci_high >= 0)
    d = disp[disp.kind.eq("gene")]
    wide = d.pivot_table(index=["target", "condition", "readout"], columns="statistic",
                         values=["estimate", "q_value"], aggfunc="first")
    wide.columns = [f"{s}_{'q' if v == 'q_value' else 'estimate'}" for v, s in wide.columns]
    out = gene.merge(wide.reset_index(), on=["target", "condition", "readout"])
    ps = h["primary_statistic"]
    primary_readout = cfg["primary_readout"]
    is_primary = out.readout.eq(primary_readout)
    out["primary_call"] = (out.mean_ci_includes_0 & out[f"{ps}_q"].lt(q)).astype("boolean").where(is_primary)
    any_disp = np.zeros(len(out), bool)
    for s in h["exploratory_or_statistics"]:
        any_disp |= out[f"{s}_q"].lt(q).to_numpy()
    out["exploratory_or_call"] = out.mean_ci_includes_0 & any_disp
    dep = depth[depth.kind.eq("gene") & depth.statistic.eq(ps)][["target", "condition", "readout", "q_value_depth"]]
    out = out.merge(dep.rename(columns={"q_value_depth": f"{ps}_q_depth_adjusted"}),
                    on=["target", "condition", "readout"], how="left")
    out["survives_depth"] = out[f"{ps}_q_depth_adjusted"].lt(q)
    cl = cluster[cluster.kind.eq("gene")].pivot_table(index=["target", "condition"], columns="statistic",
                                                      values=["cluster_ci_low", "cluster_ci_high"], aggfunc="first")
    cl.columns = [f"cluster_{s}_{v.replace('cluster_', '')}" for v, s in cl.columns]
    out = out.merge(cl.reset_index().assign(readout=primary_readout), on=["target", "condition", "readout"],
                    how="left")
    lvr_lo, lvr_hi = out[f"cluster_{ps}_ci_low"], out[f"cluster_{ps}_ci_high"]
    is_primary = out.readout.eq(primary_readout)
    out["survives_cluster"] = (((lvr_lo > 0) | (lvr_hi < 0)) & (out.cluster_mean_ci_low <= 0)
                               & (out.cluster_mean_ci_high >= 0)).astype("boolean").where(is_primary)
    keep = out.readout.eq(primary_readout) | out.exploratory_or_call
    out = out[keep].copy()
    out["readout_role"] = out.readout.map(gs.ROLES)
    return out.sort_values(["readout_role", "readout", "condition", "target"]).reset_index(drop=True)


def background_seed_spread(X, var, meta, programs, library_symbols, targets, pairs, cfg):
    """Spread of report contrasts when only the background-control draw seed changes."""
    r = cfg["reproduction"]
    contrasts = define_contrasts(targets, pairs)
    use = meta.target.isin([NTC, *targets]).to_numpy()
    groups = meta[use].groupby(["condition", "target"]).groups
    rows = []
    for k in range(r["n_background_seeds"]):
        scores, _, _ = score_cells(X, var, meta, programs, library_symbols, cfg["seed"] + 1000 + k)
        sub = scores.loc[use, r["readouts"]]
        means = {key: sub.loc[idx].mean() for key, idx in groups.items()}
        for c in contrasts:
            value = sum(w * means[key] for key, w in c["terms"].items())
            rows += [dict(contrast_id=c["contrast_id"], outcome=o, seed_index=k, estimate=value[o]) for o in sub]
    frame = pd.DataFrame(rows)
    return frame.groupby(["contrast_id", "outcome"]).estimate.agg(
        background_seed_sd="std", background_seed_min="min", background_seed_max="max").reset_index()


def reproduction_table(mean_fx, cfg, spread=None):
    old = pd.read_csv(ROOT / cfg["paths"]["report_effects"])
    r = cfg["reproduction"]
    old = old[old.target.isin(r["targets"] + [NTC]) & old.outcome.isin(r["readouts"])]
    new = mean_fx.rename(columns={"readout": "outcome"})
    merged = old.merge(new[["contrast_id", "outcome", "estimate", "ci_low", "ci_high", "iid_cell_se",
                            "group_cell_counts"]],
                       on=["contrast_id", "outcome"], suffixes=("_report", "_new"), how="left")
    merged["diff"] = merged.estimate_new - merged.estimate_report
    merged["z_vs_report_se"] = merged["diff"] / merged.iid_cell_se_report
    merged["new_in_report_ci"] = merged.estimate_new.between(merged.ci_low_report, merged.ci_high_report)
    merged["same_cells"] = merged.group_cell_counts_report == merged.group_cell_counts_new
    cols = ["contrast_id", "kind", "target", "comparator", "condition", "outcome", "estimate_report", "ci_low_report",
            "ci_high_report", "iid_cell_se_report", "estimate_new", "ci_low_new", "ci_high_new", "diff",
            "z_vs_report_se", "new_in_report_ci", "same_cells"]
    if spread is not None:
        merged = merged.merge(spread, on=["contrast_id", "outcome"], how="left")
        merged["report_within_seed_range"] = merged.estimate_report.between(merged.background_seed_min,
                                                                             merged.background_seed_max)
        cols += ["background_seed_sd", "background_seed_min", "background_seed_max", "report_within_seed_range"]
    return merged[cols].rename(columns={"outcome": "readout"})


def target_in_program(X, var, meta, programs, overlap, library_symbols, targets, cfg):
    rows = []
    for row in overlap.itertuples():
        info = programs[row.program]
        reduced = {row.program: {**info, "genes": [g for g in info["genes"] if g != row.target_symbol]}}
        score, _, _ = score_cells(X, var, meta, reduced, library_symbols, cfg["seed"])
        sub = meta.target.isin([NTC, row.target])
        fx, _ = mean_contrasts(score[sub], meta[sub], [row.target], [], dict(cfg, bootstrap=max(100, cfg["bootstrap"])))
        fx = fx[fx.target.eq(row.target) & fx.kind.isin(["gene", "interaction", "drug"])]
        for f in fx.itertuples():
            rows.append(dict(target=row.target, program=row.program, statistic="mean", kind=f.kind,
                             condition=f.condition, estimate_without_target_gene=f.estimate,
                             ci_low=f.ci_low, ci_high=f.ci_high))
        for condition in cfg["conditions"]:
            x = score.loc[sub & meta.condition.eq(condition) & meta.target.eq(row.target), row.program]
            y = score.loc[sub & meta.condition.eq(condition) & meta.target.eq(NTC), row.program]
            rows.append(dict(target=row.target, program=row.program, statistic="log_var_ratio", kind="gene",
                             condition=condition, estimate_without_target_gene=ds.log_var_ratio(x, y)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- main
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "aim1.json"))
    parser.add_argument("--targets", help="Comma-separated subset of eligible targets (smoke test)")
    parser.add_argument("--bootstrap", type=int)
    parser.add_argument("--permutations", type=int)
    parser.add_argument("--run-name")
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    if args.bootstrap:
        cfg["bootstrap"] = cfg["cluster_bootstrap"]["repeats"] = args.bootstrap
    if args.permutations:
        cfg["permutations"] = args.permutations
    P = cfg["paths"]
    run = args.run_name or ("smoke_aim1" if args.targets else "aim1")
    run_dir = ROOT / P["outputs_dir"] / run
    run_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(run_dir / "log.txt", mode="w", encoding="utf-8")])
    with open(run_dir / "config.json", "w", encoding="utf-8") as f:
        json.dump(dict(cfg, run_name=run, targets_subset=args.targets, argv=sys.argv[1:]), f, indent=2)
    t0 = time.time()
    progress = lambda msg: LOG.info("%s (%.0fs)", msg, time.time() - t0)

    split = load_split(ROOT / P["target_split"])
    eligible = sorted(split.target)
    targets = eligible if not args.targets else args.targets.split(",")
    if set(targets) - set(eligible):
        raise ValueError(f"Not eligible: {sorted(set(targets) - set(eligible))}")
    write = Writer(run_dir if args.targets else ROOT / P["summary_dir"], split)

    counts, var_all, obs = load_data(cfg)
    keep = gs.scoring_genes(var_all)
    var = var_all[keep]
    X = log_cp10k(counts[:, keep], obs.total_counts.to_numpy(float))
    del counts
    meta = obs[["condition", "target", "guide", "lig1_well", "total_counts", "target_gene_id"]].copy()
    progress(f"loaded {X.shape[0]} cells x {X.shape[1]} scoring genes")

    library_symbols, target_gene, target_symbol = library_genes(cfg, var, eligible)
    LOG.info("library target genes excluded from background controls: %d", len(library_symbols))

    spec, mapping = load_programs(ROOT / P["programs"])
    programs = gs.programs_for_scoring(mapping)
    scores, raw_scores, audit = score_cells(X, var, meta, programs, library_symbols, cfg["seed"])
    progress("program scores")

    ntc_dmso = (meta.condition.eq("DMSO") & meta.target.eq(NTC)).to_numpy()
    detect = pd.Series(np.asarray((X[ntc_dmso] > 0).mean(axis=0)).ravel(), index=var.symbol.to_numpy())
    genes_table = mapping.assign(detection_dmso_ntc=mapping.symbol.map(detect),
                                 is_library_target=mapping.symbol.isin(library_symbols),
                                 is_eligible_target=mapping.symbol.isin(set(target_symbol.values())))
    write(genes_table, "aim1_program_genes.csv", per_target=False)
    audit_rows = [dict(program=a["program"], readout_role=gs.ROLES[a["program"]], n_requested=a["n_requested"],
                       n_mapped=a["n_mapped"], missing_genes=";".join(a["missing_genes"]),
                       n_detectable_dmso_ntc=a["detectable_DMSO_NTC"], n_background_controls=len(a["controls"]),
                       center=a["center"], scale=a["scale"], raw_center=a["raw_center"], raw_scale=a["raw_scale"])
                  for a in audit]
    write(pd.DataFrame(audit_rows), "aim1_program_audit.csv", per_target=False)
    pd.concat([meta[["condition", "target", "lig1_well", "total_counts"]], scores], axis=1).to_csv(
        run_dir / "scores.csv.gz")  # cell level: outputs/ only, never committed
    overlap = gs.target_overlap(programs, {t: target_symbol[t] for t in targets})
    LOG.info("perturbed gene inside a program: %s", overlap.to_dict("records"))

    use = meta.target.isin([NTC, *targets]).to_numpy()
    S, M = scores[use], meta[use]
    r = cfg["reproduction"]
    pairs = [p for p in r["pairs"] if set(p) <= set(targets)]
    mean_fx, _ = mean_contrasts(S, M[["condition", "target"]], targets, pairs, cfg)
    progress("mean contrasts")
    flag = {(row.target, row.program) for row in overlap.itertuples()}
    mean_main = mean_fx[~mean_fx.kind.str.startswith("pair")].copy()
    mean_main["target_gene_in_program"] = [(t, o) in flag for t, o in zip(mean_main.target, mean_main.readout)]
    write(mean_main, "aim1_mean_contrasts.csv")
    write(group_summaries(S, M), "aim1_group_summaries.csv")

    disp = ds.dispersion_contrasts(S, M, targets, cfg["bootstrap"], cfg["permutations"], cfg["seed"],
                                   tuple(cfg["conditions"]), progress)
    disp["readout_role"] = disp.readout.map(gs.ROLES)
    disp["target_gene_in_program"] = [(t, o) in flag for t, o in zip(disp.target, disp.readout)]
    write(disp, "aim1_dispersion_contrasts.csv")

    depth, _ = depth_sensitivity(S, M, targets, disp, cfg, progress)
    write(depth, "aim1_dispersion_depth_sensitivity.csv")

    cluster = cluster_table(S, M, targets, mean_fx, disp, cfg)
    write(cluster, "aim1_cluster_bootstrap_fufa.csv")
    progress("cluster bootstrap")

    highlights = highlight_table(mean_fx, disp, depth, cluster, targets, cfg)
    write(highlights, "aim1_highlight_targets.csv")

    rep_targets = [t for t in r["targets"] if t in targets]
    spread = background_seed_spread(X, var, meta[["condition", "target"]], programs, library_symbols,
                                    rep_targets, pairs, cfg)
    progress("background-seed spread")
    rep = reproduction_table(mean_fx, cfg, spread)
    write(rep[rep.estimate_new.notna()], "aim1_report_reproduction.csv")

    if len(overlap):
        tip = target_in_program(X, var, meta[["condition", "target"]], programs, overlap, library_symbols,
                                targets, cfg)
        prim = pd.concat([
            mean_fx[mean_fx.kind.isin(["gene", "interaction", "drug"])].assign(statistic="mean"),
            disp[disp.statistic.eq("log_var_ratio") & disp.kind.eq("gene")]])
        prim = prim.rename(columns={"readout": "program", "estimate": "estimate_primary"})
        tip = tip.merge(prim[["target", "program", "statistic", "kind", "condition", "estimate_primary"]],
                        on=["target", "program", "statistic", "kind", "condition"], how="left")
        write(tip, "aim1_target_in_program_sensitivity.csv")
        progress("target-in-program sensitivity")

    e = cfg["expression"]
    flags = detected_genes(X, meta, e["min_ntc_detection"], tuple(cfg["conditions"]))
    cols = np.flatnonzero(flags.any(axis=1).to_numpy())
    genes = var.iloc[cols].reset_index().rename(columns={var.index.name or "index": "gene_id"})
    LOG.info("expression genes detected in >= %.0f%% of NTC cells (either condition): %d",
             100 * e["min_ntc_detection"], len(cols))
    expr = expression_contrasts(X[use][:, cols], genes, M, targets, target_gene, tuple(cfg["conditions"]))
    expr.to_csv(run_dir / "expression_contrasts.csv.gz", index=False)
    summary, tops = summarize_expression(expr, e["absdiff_cutoff"], e["q"], e["top_genes"])
    summary["n_genes_detected_any_condition"] = len(cols)
    write(summary, "aim1_expression_summary.csv")
    write(tops, "aim1_expression_top_genes.csv")
    progress("expression contrasts; done")


if __name__ == "__main__":
    main()
