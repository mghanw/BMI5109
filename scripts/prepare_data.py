"""Phase 1 data preparation and QC for PerturbFate A375 Dual-omics RNA (GSE291147).

Usage:
    python scripts/prepare_data.py [--config config/prepare_data.json] [--prefix <library prefix>]
                                   [--threshold T] [--rule threshold|paper] [--run-name NAME] [--seed S]

Steps: gene annotation from a pinned GENCODE GTF, construct assignment from sgRNA UMI counts,
first-ligation (sample) well parsing, per drug x target knockdown QC and eligibility, a
well-level expression batch check, and calibration folds made of whole sample wells.
Cell-level outputs go to data/processed/, aggregate tables to results/summary/. With --prefix
(smoke test) every output goes to outputs/<run_name>/ instead.
"""
import argparse
import copy
import gzip
import hashlib
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

NTC = "NTC"
NTC_LABEL = "NO-TARGET"
# Outdated library symbols -> GENCODE v38 gene (unversioned Ensembl ID). The library lists hg19
# coordinates for these constructs; each lies at the promoter of the gene given here.
TARGET_ALIASES = {
    "ATP5G3": ("ATP5MC3", "ENSG00000154518"),
    "CCDC101": ("SGF29", "ENSG00000176476"),
    "SGK223": ("PRAG1", "ENSG00000275342"),
    "DSCR3": ("VPS26C", "ENSG00000157538"),
    "CASC10": ("MIR1915HG", "ENSG00000204682"),
    "KIRREL": ("KIRREL1", "ENSG00000183853"),
}
CLASSES = ["singlet", "multiplet", "unassigned_low", "no_sgrna"]
LOG = logging.getLogger("prepare_data")


# ---------------------------------------------------------------- configuration
def load_config(path, overrides=None):
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    cfg = copy.deepcopy(cfg)
    for key, value in (overrides or {}).items():
        node = cfg
        *parents, leaf = key.split(".")
        for p in parents:
            node = node[p]
        if leaf not in node:
            raise KeyError(f"Unknown config key: {key}")
        node[leaf] = value
    return cfg


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- genes
def load_gene_table(gtf_path):
    """Gene-level records of a GENCODE GTF: gene_id (versioned), symbol, gene_type, chrom."""
    rows = []
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "gene":
                continue
            attr = fields[8]
            rows.append((re.search(r'gene_id "([^"]+)"', attr).group(1),
                         re.search(r'gene_name "([^"]+)"', attr).group(1),
                         re.search(r'gene_type "([^"]+)"', attr).group(1), fields[0]))
    return pd.DataFrame(rows, columns=["gene_id", "symbol", "gene_type", "chrom"])


def annotate_genes(var_names, genes):
    """Annotate matrix features; every versioned ID must exist in the pinned GTF."""
    table = genes.drop_duplicates("gene_id").set_index("gene_id")
    missing = pd.Index(var_names).difference(table.index)
    if len(missing):
        raise ValueError(f"{len(missing)} matrix gene IDs absent from the GTF, e.g. {list(missing[:3])}")
    out = table.loc[list(var_names)].reset_index()
    out["is_par_y"] = out.gene_id.str.endswith("_PAR_Y")
    return out


# ---------------------------------------------------------------- constructs
def parse_construct(name):
    """'ACTN4_1' -> 'ACTN4'; 'NO-TARGET_144' -> 'NTC'. The suffix is a global library index."""
    target, index = name.rsplit("_", 1)
    if not index.isdigit():
        raise ValueError(f"Unexpected construct name: {name}")
    return NTC if target == NTC_LABEL else target


def resolve_targets(library, genes):
    """Map each targeting construct to one non-PAR gene (versioned ID) in the annotation.

    Order: Ensembl ID embedded in sgRNA_pair_id, then exact symbol, then TARGET_ALIASES.
    """
    core = genes[~genes.gene_id.str.endswith("_PAR_Y")].copy()
    core["stable_id"] = core.gene_id.str.split(".").str[0]
    by_id = core.set_index("stable_id")
    rows = []
    for _, r in library.iterrows():
        target = parse_construct(r["id"])
        if target == NTC:
            rows.append(dict(guide=r["id"], target=NTC, target_gene_id=None, target_symbol=None, source="ntc"))
            continue
        stable = re.search(r"(ENSG\d+)", str(r.get("sgRNA_pair_id", "")))
        if stable and stable.group(1) in by_id.index:
            hit, source = by_id.loc[[stable.group(1)]], "library_ensembl_id"
        elif (core.symbol == target).sum() == 1:
            hit, source = core[core.symbol == target], "symbol"
        elif target in TARGET_ALIASES and TARGET_ALIASES[target][1] in by_id.index:
            hit, source = by_id.loc[[TARGET_ALIASES[target][1]]], "alias"
        else:
            raise ValueError(f"Cannot resolve target {target} ({r['id']}) to a unique gene")
        if len(hit) != 1:
            raise ValueError(f"Target {target} resolves to {len(hit)} genes")
        rows.append(dict(guide=r["id"], target=target, target_gene_id=hit.gene_id.iloc[0],
                         target_symbol=hit.symbol.iloc[0], source=source))
    out = pd.DataFrame(rows)
    targeted = out[out.target != NTC]
    if targeted.target.duplicated().any() or targeted.target_gene_id.duplicated().any():
        raise ValueError("Targets and target genes must be one-to-one")
    return out


# ---------------------------------------------------------------- barcodes
def parse_barcode(cell_names, lig1_whitelist):
    """Split '<prefix>.<bc1 10 nt><bc2 10 nt>' into PCR prefix and first-ligation (sample) well."""
    names = pd.Series(cell_names, dtype=str)
    parts = names.str.split(".", n=1, expand=True)
    if parts.shape[1] != 2 or parts[1].str.len().ne(20).any():
        raise ValueError("Cell names must look like <prefix>.<20-nt barcode>")
    index = {bc: i for i, bc in enumerate(lig1_whitelist)}
    lig1 = parts[1].str[:10].map(index)
    if lig1.isna().any():
        raise ValueError(f"{int(lig1.isna().sum())} cells have a first-ligation barcode outside the whitelist")
    return pd.DataFrame({"prefix": parts[0].values, "pcr_well": parts[0].str[-3:].values,
                         "lig1_index": lig1.astype(int).values}, index=names.values)


# ---------------------------------------------------------------- assignment
def top_two(X):
    """Row-wise total, maximum, argmax and second-largest value of a sparse count matrix."""
    X = sp.csr_matrix(X)
    total = np.asarray(X.sum(axis=1)).ravel()
    top = np.asarray(X.max(axis=1).todense()).ravel()
    arg = np.asarray(X.argmax(axis=1)).ravel()
    rest = X - sp.csr_matrix((top, (np.arange(X.shape[0]), arg)), shape=X.shape)
    rest.eliminate_zeros()
    second = np.asarray(rest.max(axis=1).todense()).ravel()
    return total, top, arg, second


def choose_threshold(top_umi, bounds):
    """Presence threshold T at the valley of the top-construct UMI density on log scale.

    For integer UMI t the log-scale density is proportional to count(top == t) * (t + 1);
    the valley is searched within the inclusive bounds after 3-point smoothing.
    """
    lo, hi = bounds
    grid = np.arange(lo - 1, hi + 2)
    counts = np.array([(top_umi == t).sum() for t in grid], dtype=float)
    density = counts * (grid + 1)
    smooth = np.convolve(density, np.ones(3) / 3, mode="same")[1:-1]
    return int(grid[1:-1][np.argmin(smooth)])


def classify_threshold(X, T):
    n_present = np.asarray((sp.csr_matrix(X) >= T).sum(axis=1)).ravel()
    return np.select([n_present == 1, n_present >= 2], ["singlet", "multiplet"], "unassigned_low"), n_present


def classify_paper(total, top, second, min_total, max_to_total, sec_to_max):
    """Authors' rule: total UMI > min_total, top/total >= max_to_total, second/top <= sec_to_max."""
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(total > 0, top / total, 0)
        ratio = np.where(top > 0, second / top, np.inf)
    ok = total > min_total
    return np.select([ok & (frac >= max_to_total) & (ratio <= sec_to_max), ok],
                     ["singlet", "multiplet"], "unassigned_low")


def assign_cells(cells, sg_cells, sg_X, guides, cfg):
    """Cell-level assignment table indexed like `cells` (cells without sgRNA -> no_sgrna)."""
    a = cfg["assignment"]
    total, top, arg, second = top_two(sg_X)
    T = a["threshold"] if a["threshold"] is not None else choose_threshold(top, a["threshold_bounds"])
    cls_t, n_present = classify_threshold(sg_X, T)
    cls_p = classify_paper(total, top, second, a["paper_min_total_umi"], a["paper_max_to_total"], a["paper_sec_to_max"])
    sg = pd.DataFrame({"sg_total": total, "top_guide": np.asarray(guides)[arg], "top_umi": top,
                       "second_umi": second, "n_present": n_present,
                       "class_threshold": cls_t, "class_paper": cls_p}, index=pd.Index(sg_cells))
    out = sg.reindex(pd.Index(cells))
    for col in ["class_threshold", "class_paper"]:
        out[col] = out[col].fillna("no_sgrna")
    for col in ["sg_total", "top_umi", "second_umi", "n_present"]:
        out[col] = out[col].fillna(0).astype(int)
    out["assignment"] = out["class_" + a["rule"]]
    out["guide"] = np.where(out.assignment.eq("singlet"), out.top_guide, None)
    sens = {t: classify_threshold(sg_X, t)[0] for t in a["sensitivity_thresholds"]}
    return out, T, sens


def assignment_summary(assign, condition, T, sens, sg_cells):
    rows = []
    def add(rule, classes):
        for (cond, cls), n in pd.crosstab(condition, classes).stack().items():
            rows.append(dict(rule=rule, condition=cond, assignment_class=cls, n_cells=int(n)))
    add(f"threshold_T{T}_auto", assign.class_threshold)
    add("paper_max0.5_sec0.3", assign.class_paper)
    for t, cls in sens.items():
        add(f"threshold_T{t}", pd.Series(cls, index=pd.Index(sg_cells)).reindex(condition.index).fillna("no_sgrna"))
    out = pd.DataFrame(rows)
    out["assignment_class"] = pd.Categorical(out.assignment_class, CLASSES)
    return out.sort_values(["rule", "condition", "assignment_class"]).reset_index(drop=True)


# ---------------------------------------------------------------- knockdown
def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(len(pv))
        tmp[order] = np.minimum(ranked, 1)
        q[ok] = tmp
    return q


def knockdown_table(cp10k, obs, target_genes, cfg):
    """Per condition x target knockdown QC on singlets.

    cp10k: cells x targets array (CP10k of each target's own gene), rows aligned with obs.
    obs: columns condition, target (singlets only). target_genes: target -> gene_id, ordered as cp10k columns.
    """
    k = cfg["knockdown"]
    targets = list(target_genes)
    rows = []
    for condition in sorted(obs.condition.unique()):
        in_cond = obs.condition.eq(condition).to_numpy()
        ntc = in_cond & obs.target.eq(NTC).to_numpy()
        rows.append(dict(condition=condition, target=NTC, target_gene_id=None, n_cells=int(ntc.sum())))
        for j, target in enumerate(targets):
            kd = in_cond & obs.target.eq(target).to_numpy()
            x_kd, x_ntc = cp10k[kd, j], cp10k[ntc, j]
            det_ntc = float((x_ntc > 0).mean()) if len(x_ntc) else np.nan
            mean_ntc = x_ntc.mean() if len(x_ntc) else np.nan
            ratio = x_kd.mean() / mean_ntc if len(x_kd) and mean_ntc > 0 else np.nan
            p = np.nan
            if len(x_kd) >= 2 and len(x_ntc) >= 2:
                a, b = np.log1p(x_kd), np.log1p(x_ntc)
                if a.var() + b.var() > 0:
                    p = stats.ttest_ind(a, b, equal_var=False, alternative="less").pvalue
            rows.append(dict(condition=condition, target=target, target_gene_id=target_genes[target],
                             n_cells=int(kd.sum()), mean_cp10k_kd=x_kd.mean() if len(x_kd) else np.nan,
                             mean_cp10k_ntc=mean_ntc, detection_kd=float((x_kd > 0).mean()) if len(x_kd) else np.nan,
                             detection_ntc=det_ntc, knockdown_remaining_ratio=ratio, welch_p=p))
    out = pd.DataFrame(rows)
    out["evaluable"] = out.detection_ntc >= k["min_ntc_detection"]
    out["welch_q"] = np.nan
    for condition, idx in out[out.target.ne(NTC)].groupby("condition").groups.items():
        sub = out.loc[idx]
        out.loc[idx, "welch_q"] = bh(np.where(sub.evaluable, sub.welch_p, np.nan))
    out = out.astype({"evaluable": bool})

    def status(max_ratio):
        s = np.select([out.n_cells < 2, ~out.evaluable,
                       (out.knockdown_remaining_ratio <= max_ratio) & (out.welch_q < k["q"])],
                      ["insufficient_cells", "not_evaluable", "adequate"], "inadequate")
        return np.where(out.target.eq(NTC), None, s)
    out["kd_status"] = status(k["max_ratio"])
    out["kd_status_ratio07"] = status(k["sensitivity_max_ratio"])
    out["paper_kd_pass"] = np.where(out.target.eq(NTC), None, out.knockdown_remaining_ratio < k["paper_max_ratio"])
    return add_eligibility(out, k)


def add_eligibility(qc, k):
    """Target is eligible if > min_cells singlets AND adequate knockdown in EVERY condition."""
    t = qc[qc.target.ne(NTC)]
    n_conditions = qc.condition.nunique()

    def per_target(status_col, pass_fn):
        ok = t.assign(ok=(t.n_cells > k["min_cells"]) & pass_fn(t[status_col]))
        g = ok.groupby("target").ok.agg(["all", "size"])
        return (g["all"] & g["size"].eq(n_conditions)).rename(status_col)

    elig = pd.concat([per_target("kd_status", lambda s: s.eq("adequate")),
                      per_target("kd_status_ratio07", lambda s: s.eq("adequate")),
                      per_target("paper_kd_pass", lambda s: s.astype(bool))], axis=1)
    elig.columns = ["eligible", "eligible_ratio07", "eligible_paper_rule"]
    elig["eligibility_differs_from_paper"] = elig.eligible != elig.eligible_paper_rule
    out = qc.merge(elig, left_on="target", right_index=True, how="left")
    return out


# ---------------------------------------------------------------- batch check
def candidate_splits(lig1_index):
    """Neutral two-group layouts of the 48 sample wells of one condition (index within block)."""
    i = np.asarray(lig1_index) % 48
    return {"contiguous_halves": i < 24, "mod12_halves": (i % 12) < 6,
            "parity": i % 2 == 0, "mod8_halves": (i % 8) < 4}


def pseudo_f(scores, labels):
    groups = np.unique(labels)
    grand = scores.mean(axis=0)
    ssb = sum((labels == g).sum() * ((scores[labels == g].mean(axis=0) - grand) ** 2).sum() for g in groups)
    ssw = sum(((scores[labels == g] - scores[labels == g].mean(axis=0)) ** 2).sum() for g in groups)
    return (ssb / (len(groups) - 1)) / (ssw / (len(scores) - len(groups)))


def well_batch_check(pseudobulk_counts, well_info, cfg, rng):
    """PCA + clustering of well-level pseudobulk profiles within each condition.

    pseudobulk_counts: wells x genes summed raw counts; well_info: columns condition, lig1_index,
    median_umi (rows aligned). Returns (per-well table, per-condition split tests).
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    b = cfg["batch_check"]
    counts = np.asarray(pseudobulk_counts, float)
    cp10k = counts / counts.sum(axis=1, keepdims=True) * 1e4
    wells, tests = [], []
    for condition in sorted(well_info.condition.unique()):
        m = well_info.condition.eq(condition).to_numpy()
        x = cp10k[m]
        keep = x.mean(axis=0) >= b["min_cp10k"]
        logx = np.log1p(x[:, keep])
        hvg = np.argsort(logx.var(axis=0))[::-1][:b["n_hvg"]]
        z = logx[:, hvg] - logx[:, hvg].mean(axis=0)
        u, s, _ = np.linalg.svd(z, full_matrices=False)
        n_pcs = min(b["n_pcs"], len(s))
        scores = u[:, :n_pcs] * s[:n_pcs]
        explained = s ** 2 / (s ** 2).sum()
        clusters = fcluster(linkage(scores, "average"), 2, "maxclust")
        info = well_info[m].reset_index(drop=True)
        frame = info.assign(**{f"PC{i + 1}": scores[:, i] for i in range(n_pcs)}, hclust2=clusters)
        wells.append(frame)
        depth = info.median_umi.to_numpy(float)
        for name, labels in {**candidate_splits(info.lig1_index), "hclust2": clusters == 1}.items():
            labels = np.asarray(labels).astype(int)
            observed = pseudo_f(scores, labels)
            null = np.array([pseudo_f(scores, rng.permutation(labels)) for _ in range(b["n_permutations"])])
            tests.append(dict(condition=condition, split=name, n_wells=len(labels),
                              group_sizes=f"{(labels == 0).sum()}/{(labels == 1).sum()}",
                              pseudo_f=observed, perm_p=(1 + (null >= observed).sum()) / (1 + len(null)),
                              depth_diff=depth[labels == 1].mean() - depth[labels == 0].mean(),
                              pc1_var=explained[0], pc2_var=explained[1] if len(explained) > 1 else np.nan,
                              pc1_depth_spearman=stats.spearmanr(scores[:, 0], depth)[0]))
    return pd.concat(wells, ignore_index=True), pd.DataFrame(tests)


def calibration_folds(lig1_index, condition, n_folds, seed):
    """Assign whole sample wells to calibration folds, balanced within each condition."""
    wells = pd.DataFrame({"lig1_index": lig1_index, "condition": condition}).drop_duplicates()
    rng = np.random.default_rng(seed)
    fold = {}
    for _, g in wells.sort_values("lig1_index").groupby("condition"):
        order = rng.permutation(g.lig1_index.to_numpy())
        fold.update({w: i % n_folds for i, w in enumerate(order)})
    return pd.Series(lig1_index).map(fold).to_numpy()


# ---------------------------------------------------------------- safety
def assert_no_cell_ids(frame):
    """Aggregate outputs must not contain barcode sequences or cell names."""
    pattern = re.compile(r"[ACGT]{8,}")
    for col in frame.columns:
        values = frame[col].dropna().astype(str)
        if values.map(lambda v: bool(pattern.search(v))).any():
            raise ValueError(f"Column {col!r} contains barcode-like strings")


# ---------------------------------------------------------------- main
def well_name(i):
    return f"L1_{int(i) + 1:02d}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/prepare_data.json")
    parser.add_argument("--prefix", help="Smoke test on one library prefix; outputs only to outputs/<run>/")
    parser.add_argument("--threshold", type=int)
    parser.add_argument("--rule", choices=["threshold", "paper"])
    parser.add_argument("--run-name")
    parser.add_argument("--seed", type=int)
    args = parser.parse_args(argv)
    overrides = {k: v for k, v in {"assignment.threshold": args.threshold, "assignment.rule": args.rule,
                                   "seed": args.seed}.items() if v is not None}
    cfg = load_config(args.config, overrides)
    P = cfg["paths"]
    run = args.run_name or ("smoke_" + args.prefix.split("_")[-1] if args.prefix else "phase1_full")
    run_dir = Path(P["outputs_dir"]) / run
    run_dir.mkdir(parents=True, exist_ok=True)
    processed = run_dir if args.prefix else Path(P["processed_dir"])
    summary = run_dir if args.prefix else Path(P["summary_dir"])
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(run_dir / "log.txt", mode="w", encoding="utf-8")])
    with open(run_dir / "config.json", "w", encoding="utf-8") as f:
        json.dump(dict(cfg, run_name=run, prefix=args.prefix, argv=sys.argv[1:]), f, indent=2)
    t0 = time.time()
    import anndata as ad

    meta = pd.read_csv(P["metadata"], dtype={"cell_names": str})
    meta["condition"] = meta.conditions.str.strip().map(cfg["conditions"])
    if meta.condition.isna().any():
        raise ValueError("Unmapped condition labels in metadata")
    meta = meta.set_index("cell_names")
    if args.prefix:
        meta = meta[meta.index.str.split(".").str[0] == args.prefix]
        if meta.empty:
            raise ValueError(f"No cells with prefix {args.prefix}")
    LOG.info("cells: %d (%s)", len(meta), meta.condition.value_counts().to_dict())

    whitelist = pd.read_csv(P["lig1_whitelist"], header=None)[0].tolist()
    obs = parse_barcode(meta.index, whitelist).join(meta[["UMI_counts", "gene_number", "condition"]])
    purity = obs.groupby("lig1_index").condition.nunique()
    if purity.gt(1).any():
        raise ValueError("A first-ligation well holds more than one condition")
    LOG.info("sample wells: %s", obs.groupby("condition").lig1_index.nunique().to_dict())

    genes_all = load_gene_table(P["gtf"])
    library = pd.read_csv(P["library"])
    targets = resolve_targets(library, genes_all)
    LOG.info("target resolution: %s", targets[targets.target.ne(NTC)].source.value_counts().to_dict())

    sg = ad.read_h5ad(P["sgrna"])
    sg = sg[sg.obs_names.isin(obs.index)]
    unknown = set(sg.var_names) - set(targets.guide)
    if unknown:
        raise ValueError(f"Guides absent from the library: {sorted(unknown)[:5]}")
    assign, T, sens = assign_cells(obs.index, sg.obs_names, sp.csr_matrix(sg.X), sg.var_names, cfg)
    LOG.info("presence threshold T=%d; %s", T, assign.assignment.value_counts().to_dict())
    obs = obs.join(assign)
    guide_target = targets.set_index("guide").target
    obs["target"] = obs.guide.map(guide_target)
    obs["calib_fold"] = calibration_folds(obs.lig1_index, obs.condition, cfg["calibration_folds"], cfg["seed"])
    obs["lig1_well"] = obs.lig1_index.map(well_name)

    if args.prefix:
        backed = ad.read_h5ad(P["rna"], backed="r")
        rna = backed[np.flatnonzero(backed.obs_names.isin(obs.index))].to_memory()
    else:
        rna = ad.read_h5ad(P["rna"])
    if not obs.index.isin(rna.obs_names).all():
        raise ValueError("Metadata cells missing from the RNA matrix")
    rna = rna[obs.index]
    X = sp.csr_matrix(rna.X)
    LOG.info("RNA matrix %s, nnz=%d (%.0fs)", X.shape, X.nnz, time.time() - t0)

    var = annotate_genes(rna.var_names, genes_all)
    par = var.is_par_y.to_numpy()
    par_nonzero = int((X[:, np.flatnonzero(par)] != 0).sum())
    LOG.info("PAR_Y rows: %d, nonzero entries: %d", par.sum(), par_nonzero)
    total = np.asarray(X.sum(axis=1)).ravel()
    obs["total_counts"] = total
    LOG.info("total counts == metadata UMI_counts: %.4f", np.mean(total == obs.UMI_counts.to_numpy()))

    singlet = obs.assignment.eq("singlet").to_numpy()
    tg = targets[targets.target.ne(NTC)].set_index("target").target_gene_id
    col = pd.Series(np.arange(len(var)), index=var.gene_id)
    Xs = X[singlet][:, col[tg.values].values].toarray().astype(float)
    cp10k = Xs / total[singlet, None] * 1e4
    qc = knockdown_table(cp10k, obs.loc[singlet, ["condition", "target"]], tg.to_dict(), cfg)
    qc = qc.merge(targets[["target", "target_symbol", "source"]].drop_duplicates("target")
                  .rename(columns={"source": "target_gene_source"}), on="target", how="left")
    LOG.info("eligible targets: %d (ratio<=0.7: %d; paper rule: %d)",
             qc.drop_duplicates("target").eligible.eq(True).sum(),
             qc.drop_duplicates("target").eligible_ratio07.eq(True).sum(),
             qc.drop_duplicates("target").eligible_paper_rule.eq(True).sum())

    # Well-level batch check: NTC singlets and all singlets.
    rng = np.random.default_rng(cfg["seed"])
    batch_wells, batch_tests = [], []
    for population, mask in [("NTC", singlet & obs.target.eq(NTC).to_numpy()), ("all_singlets", singlet)]:
        idx = obs.lig1_index.to_numpy()[mask]
        uniq = np.unique(idx)
        ind = sp.csr_matrix((np.ones(mask.sum()), (np.searchsorted(uniq, idx), np.arange(mask.sum()))),
                            shape=(len(uniq), mask.sum()))
        pb = (ind @ X[mask]).toarray()
        sub = obs[mask].groupby("lig1_index")
        info = pd.DataFrame({"lig1_index": uniq, "condition": sub.condition.first().loc[uniq].values,
                             "n_cells": sub.size().loc[uniq].values,
                             "median_umi": sub.total_counts.median().loc[uniq].values,
                             "iqr_umi": (sub.total_counts.quantile(.75) - sub.total_counts.quantile(.25)).loc[uniq].values,
                             "median_genes": sub.gene_number.median().loc[uniq].values,
                             "iqr_genes": (sub.gene_number.quantile(.75) - sub.gene_number.quantile(.25)).loc[uniq].values})
        if len(uniq) < 8:
            LOG.info("batch check skipped for %s (%d wells)", population, len(uniq))
            continue
        w, t = well_batch_check(pb, info, cfg, rng)
        batch_wells.append(w.assign(population=population))
        batch_tests.append(t.assign(population=population))
    if batch_wells:
        wells = pd.concat(batch_wells, ignore_index=True)
        wells.insert(0, "lig1_well", wells.pop("lig1_index").map(well_name))
        tests = pd.concat(batch_tests, ignore_index=True)
        LOG.info("batch split tests:\n%s", tests.to_string())

    # Validation against the earlier report.
    val = cfg["validation"]["report_singlets"]
    comp = []
    for target, (n_dmso, n_vem) in val.items():
        for condition, n_report in [("DMSO", n_dmso), ("Vem", n_vem)]:
            in_c = obs.condition.eq(condition) & obs.top_guide.map(guide_target).eq(target)
            comp.append(dict(target=target, condition=condition, report_n=n_report,
                             threshold_rule_n=int((in_c & obs.class_threshold.eq("singlet")).sum()),
                             paper_rule_n=int((in_c & obs.class_paper.eq("singlet")).sum())))
    comp = pd.DataFrame(comp)
    comp["threshold_minus_report"] = comp.threshold_rule_n - comp.report_n
    LOG.info("report comparison:\n%s", comp.to_string())

    # Outputs.
    summ = assignment_summary(obs[["class_threshold", "class_paper"]], obs.condition, T, sens, sg.obs_names)
    aggregates = {"phase1_target_qc.csv": qc, "phase1_assignment_summary.csv": summ,
                  "phase1_report_comparison.csv": comp}
    if batch_wells:
        aggregates.update({"phase1_well_batch.csv": wells, "phase1_batch_tests.csv": tests})
    summary.mkdir(parents=True, exist_ok=True)
    for name, frame in aggregates.items():
        assert_no_cell_ids(frame)
        frame.to_csv(summary / name, index=False, float_format="%.6g")
    source = dict(annotation=cfg["annotation"], gtf=Path(P["gtf"]).name, gtf_sha256=sha256(P["gtf"]),
                  gtf_url=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{cfg['annotation']['release']}/"
                          f"{Path(P['gtf']).name}",
                  matrix_genes=int(len(var)), matched_versioned=int(len(var)), gtf_genes=int(len(genes_all)),
                  par_y_rows=int(par.sum()), par_y_nonzero_entries=par_nonzero,
                  library=Path(P["library"]).name, library_sha256=sha256(P["library"]),
                  library_repo="https://github.com/JunyueCaoLab/PerturbFate@bafdfff28a6fc36854a5defe32440abc98eaabc0",
                  target_resolution=targets[targets.target.ne(NTC)].source.value_counts().to_dict(),
                  aliases={k: v[0] for k, v in TARGET_ALIASES.items()}, presence_threshold=T)
    with open(summary / "annotation_source.json", "w", encoding="utf-8") as f:
        json.dump(source, f, indent=2)

    var.to_csv(processed / "genes.csv.gz", index=False)
    cells_out = obs.drop(columns=["UMI_counts"]).rename(columns={"gene_number": "n_genes"})
    cells_out["target_gene_id"] = cells_out.target.map(tg)
    cells_out.rename_axis("cell").to_csv(processed / "cells.csv.gz")
    keep_genes = ~par if par_nonzero == 0 else np.ones(len(var), bool)
    out = ad.AnnData(X=X[singlet][:, keep_genes],
                     obs=cells_out[singlet].astype({c: "category" for c in ["prefix", "pcr_well", "condition",
                                                                           "guide", "target", "lig1_well"]}),
                     var=var[keep_genes].set_index("gene_id"))
    for c in out.obs.columns:
        if out.obs[c].dtype == object:
            out.obs[c] = out.obs[c].fillna("").astype(str)
    out.uns["presence_threshold"] = T
    out.uns["assignment_rule"] = cfg["assignment"]["rule"]
    out.write_h5ad(processed / "rna_singlets.h5ad", compression="gzip")
    LOG.info("wrote %s singlets x %d genes; done in %.0fs", out.n_obs, out.n_vars, time.time() - t0)


if __name__ == "__main__":
    main()
