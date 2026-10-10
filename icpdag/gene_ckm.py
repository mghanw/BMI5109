"""Official NMF/CSCN with intervention seeds and gene-space CKM projection."""

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import warnings

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from ._vendor.cscn.core import CSCN
from ._vendor.cscn.ckm import solve_dag_katz
from .intervention_orientation import seed_interventions, valid_extension
from .responses import response_table, evidence_for_condition


def project_to_genes(ckm_nodes, loadings, anchors):
    """Map factor CKM through the same NMF H; preserve measured anchor columns."""
    if ckm_nodes.index.has_duplicates or ckm_nodes.columns.has_duplicates:
        raise ValueError("CKM cell and node identifiers must be unique")
    if loadings.index.has_duplicates or loadings.columns.has_duplicates:
        raise ValueError("NMF factor and gene identifiers must be unique")
    if set(loadings.columns) & set(anchors):
        raise ValueError("Anchor genes must be excluded from the NMF input")
    if set(loadings.index) & set(anchors):
        raise ValueError("Factor and anchor names must be disjoint")
    expected = set(loadings.index) | set(anchors)
    if set(ckm_nodes.columns) != expected:
        raise ValueError("CKM nodes do not match the NMF factors and anchors")
    mapped = ckm_nodes.loc[:, loadings.index].to_numpy(dtype=float) @ loadings.to_numpy(dtype=float)
    genes = pd.DataFrame(mapped, index=ckm_nodes.index, columns=loadings.columns)
    return pd.concat([genes, ckm_nodes.loc[:, anchors]], axis=1)


def _infer_group(job):
    """Each job uses one condition/target reference stratum."""
    with threadpool_limits(limits=1):
        nodes, data, cfg = job["nodes"], job["data"], job["config"]
        model = CSCN(variable_names=nodes, sigmoid_score=cfg["window"],
                     significance_level=cfg["ci_alpha"], max_cond_vars=cfg["max_conditioning"],
                     query_engine="bitmap", pc_var="stable")
        model.run_core(data, use_nmf=False)
        rows, attempts, pdags, dags, values = [], [], [], [], []
        for local in job["focals"]:
            model.clear_task_cache()
            diagnostics = StringIO()
            with redirect_stdout(diagnostics):
                skeleton, separators = model._estimate_pc_skeleton(local)
            if "ICT ERROR:" in diagnostics.getvalue():
                raise RuntimeError(diagnostics.getvalue())
            collider = model._orient_unshielded_colliders(skeleton, separators)
            seeded, records, absent = seed_interventions(
                collider, job["evidence"], model._build_dag_representation
            )
            pdag = model._apply_meek_rules(seeded.copy(), apply_r4=False)
            pa = nx.to_numpy_array(pdag, nodelist=nodes, dtype=np.uint8)
            da = np.zeros_like(pa)
            error = ""
            try:
                dag = model._build_dag_representation(pdag)
                valid = valid_extension(pdag, dag)
                da = nx.to_numpy_array(dag, nodelist=nodes, dtype=np.uint8)
                if not valid:
                    error = "DAG completion inconsistent with PDAG"
            except (ValueError, RuntimeError) as exc:
                valid, error = False, str(exc)
            ckm = np.full(len(nodes), np.nan, dtype=np.float32)
            if valid:
                ptr = np.r_[0, np.cumsum(da.sum(axis=1))].astype(np.int64)
                targets = np.nonzero(da)[1].astype(np.int64)
                order = np.asarray([nodes.index(g) for g in nx.topological_sort(dag)], dtype=np.int64)
                ckm = solve_dag_katz(ptr, targets, order, data[local].astype(np.float32), cfg["katz_alpha"])
            index = int(job["indices"][local])
            rows.append(dict(cell_index=index, condition=job["condition"], target=job["target"],
                             n_reference=len(data), n_edges=skeleton.number_of_edges(),
                             valid_completion=valid, completion_error=error,
                             accepted_seeds=sum(r["status"] == "accepted_seed" for r in records),
                             absent_seed_proposals=absent))
            attempts.extend(dict(cell_index=index, **record) for record in records)
            pdags.append(pa)
            dags.append(da)
            values.append(ckm)
        return dict(rows=rows, attempts=attempts, pdag=np.asarray(pdags),
                    dag=np.asarray(dags), ckm=np.asarray(values))


def fit_gene_ckm(expression, metadata, *, anchors=None, knockdown_qc=None,
                 components=20, window=0.35, max_conditioning=4, ci_alpha=0.05,
                 katz_alpha=0.05, n_states=6, seed=20261011, workers=1):
    """Fit gene-space CKM from nonnegative, normalized cell-by-gene expression.

    Metadata shares the expression index and contains condition/target. With
    knockdown_qc, dictionary_fold marks discovery/validation cells; only discovery
    responses select seeds. NMF and reference graphs use the full input dataset.
    QC columns are condition, target, knockdown_remaining_ratio and
    eligible_both_conditions. Selected genes must include measured anchor genes.
    """
    if expression.index.has_duplicates or expression.columns.has_duplicates or metadata.index.has_duplicates:
        raise ValueError("Cell and gene identifiers must be unique")
    if set(expression.index) != set(metadata.index):
        raise ValueError("Expression and metadata must contain identical cells")
    metadata = metadata.reindex(expression.index).copy()
    if not {"condition", "target"} <= set(metadata) or metadata[["condition", "target"]].isna().any().any():
        raise ValueError("Metadata requires nonmissing condition and target")
    if not np.isfinite(expression.to_numpy(dtype=float)).all() or (expression.to_numpy() < 0).any():
        raise ValueError("NMF expression must be finite and nonnegative")
    if not 0 < window < 1 or not 0 < ci_alpha < 1 or katz_alpha < 0:
        raise ValueError("Invalid window, CI significance level or Katz alpha")
    if workers < 1 or max_conditioning < 0 or n_states < 2:
        raise ValueError("Invalid workers, conditioning order or state count")
    anchors = sorted(set(metadata.target) - {"NTC"}) if anchors is None else list(anchors)
    if len(anchors) != len(set(anchors)) or set(anchors) - set(expression.columns):
        raise ValueError("Anchors must be unique, measured expression columns")
    genes = [g for g in expression if g not in anchors]
    if not 1 <= components <= min(len(expression), len(genes)):
        raise ValueError("NMF components exceed the available cells or nonanchor genes")
    factors = [f"factor_{i:03d}" for i in range(components)]
    if set(factors) & set(anchors):
        raise ValueError("Anchor names collide with factor identifiers")
    config = dict(components=components, window=window, max_conditioning=max_conditioning,
                  ci_alpha=ci_alpha, katz_alpha=katz_alpha, n_states=n_states, seed=seed,
                  nmf_fit_scope="all input cells", reference_scope="condition x target",
                  nmf_random_state=42, state_random_state=seed,
                  orientation="colliders -> discovery intervention seeds -> Meek -> DAG",
                  projection="factor CKM @ native H; anchor CKM columns retained")
    native = CSCN(query_engine="bitmap")
    with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=min(workers, 4)):
        warnings.simplefilter("always")
        native.run_core(expression.loc[:, genes].to_numpy(dtype=np.float32), use_nmf=True,
                        nmf_components=components, nmf_max_iter=10000)
    w = pd.DataFrame(native.data, index=expression.index, columns=factors)
    h = pd.DataFrame(native.loadings.T, index=factors, columns=genes)
    node_values = pd.concat([expression.loc[:, anchors].astype(np.float32), w], axis=1)
    config["nmf_warnings"] = [str(item.message) for item in caught]
    if knockdown_qc is not None:
        if "dictionary_fold" not in metadata or not set(metadata.dictionary_fold) <= {"discovery", "validation"}:
            raise ValueError("Intervention seeds require discovery/validation dictionary_fold labels")
        if knockdown_qc.duplicated(["condition", "target"]).any():
            raise ValueError("Knockdown QC must have one row per condition and target")
        counts = metadata[metadata.dictionary_fold.eq("discovery")].groupby(["condition", "target"]).size()
        for condition in metadata.condition.unique():
            for target in ["NTC", *anchors]:
                if counts.get((condition, target), 0) < 2:
                    raise ValueError(f"Insufficient discovery cells: {condition}/{target}")
        evidence = response_table(node_values, metadata, anchors, knockdown_qc, "discovery")
    else:
        evidence = pd.DataFrame(columns=["condition", "source", "destination", "eligible_seed"])
    jobs = []
    for (condition, target), indices in metadata.groupby(["condition", "target"], sort=True).indices.items():
        if len(indices) < 2:
            raise ValueError("Each reference stratum needs at least two cells")
        proposals = evidence_for_condition(evidence, condition) if len(evidence) else []
        for start in range(0, len(indices), 64):
            jobs.append(dict(nodes=list(node_values), data=node_values.iloc[indices].to_numpy(), config=config,
                             condition=condition, target=target, indices=indices,
                             focals=list(range(start, min(start+64, len(indices)))), evidence=proposals))
    if workers == 1:
        output = [_infer_group(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            output = list(pool.map(_infer_group, jobs))
    rows = pd.DataFrame([r for result in output for r in result["rows"]])
    order = np.argsort(rows.cell_index.to_numpy())
    cells = rows.iloc[order].reset_index(drop=True)
    cells["cell_id"] = expression.index
    arrays = {key: np.concatenate([result[key] for result in output])[order]
              for key in ["pdag", "dag", "ckm"]}
    ckm_nodes = pd.DataFrame(arrays["ckm"], index=expression.index, columns=node_values.columns)
    ckm_genes = project_to_genes(ckm_nodes, h, anchors)
    valid = cells.valid_completion.to_numpy(dtype=bool)
    if not np.array_equal(np.isfinite(ckm_genes.to_numpy()).all(axis=1), valid):
        raise RuntimeError("CKM validity and DAG completion disagree")
    states = metadata.loc[valid, ["condition", "target"]].copy()
    config["state_status"] = "insufficient valid cells"
    if len(states) >= n_states:
        with threadpool_limits(limits=min(workers, 4)):
            scaled = StandardScaler().fit_transform(ckm_genes.loc[valid])
            n_pcs = min(20, *scaled.shape)
            pca = PCA(n_components=n_pcs, random_state=seed, svd_solver="randomized")
            pcs = pca.fit_transform(scaled)
            labels = KMeans(n_clusters=n_states, n_init=20, random_state=seed).fit_predict(pcs)
        states["state"] = [f"S{i+1}" for i in labels]
        states["PC1"] = pcs[:, 0]
        states["PC2"] = pcs[:, 1] if n_pcs > 1 else np.nan
        config.update(state_status="complete", explained_variance_ratio=pca.explained_variance_ratio_.tolist())
    attempts = pd.DataFrame([a for result in output for a in result["attempts"]],
                            columns=["cell_index", "source", "destination", "estimate",
                                     "standardized_estimate", "q_value", "status"])
    return dict(W=w, H=h, node_values=node_values, ckm_nodes=ckm_nodes, ckm_genes=ckm_genes,
                cells=cells, states=states, metadata=metadata, response_evidence=evidence,
                orientation_attempts=attempts, pdag=arrays["pdag"], dag=arrays["dag"], config=config)


def save_result(result, directory):
    """Save one completed analysis to an empty directory."""
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Output directory must be empty")
    directory.mkdir(parents=True, exist_ok=True)
    for name in ["W", "H", "node_values", "ckm_nodes", "ckm_genes", "metadata", "states"]:
        result[name].to_csv(directory / f"{name}.csv.gz")
    for name in ["cells", "response_evidence", "orientation_attempts"]:
        result[name].to_csv(directory / f"{name}.csv.gz", index=False)
    np.savez_compressed(directory / "graphs.npz", pdag=result["pdag"], dag=result["dag"],
                        nodes=np.asarray(result["node_values"].columns, dtype=str))
    (directory / "config.json").write_text(json.dumps(result["config"], indent=2), encoding="utf-8")
