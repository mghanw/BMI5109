"""Inputs for the SEM prediction step: panel expression, stored cell graphs, targets and folds."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .cscn_networks import decode_edges
from .cscn_panel import NTC, log_cp10k, panel_of_size


def panel_nodes(panel):
    """Gene IDs and node names (symbols) of the frozen panel."""
    if not panel.get("frozen"):
        raise ValueError("The network panel must be frozen")
    forced = [g["gene_id"] for g in panel["forced_genes"]]
    ranked = [g for g, _ in panel["fill_genes_ranked"]]
    symbol = {g["gene_id"]: g["symbol"] for g in panel["forced_genes"]}
    symbol.update({g: s for g, s in panel["fill_genes_ranked"]})
    ids = panel_of_size(forced, ranked, panel["panel_size"])
    return ids, [symbol[g] for g in ids]


def target_nodes(panel):
    """Eligible target name -> its panel node (symbol of the perturbed gene)."""
    return {g["target"]: g["symbol"] for g in panel["forced_genes"] if g["role"] == "target"}


def fold_targets(split, fold):
    """(held-out targets, training targets) of an outer fold."""
    held = sorted(split.target[split.fold.eq(fold)])
    return held, sorted(set(split.target) - set(held))


def load_panel_cells(path, gene_ids, conditions, targets):
    """log1p CP10k of the panel genes for NTC cells and the given targets' cells (dense float32).

    Rows of other targets are dropped before normalization. Metadata keeps condition, target,
    calib_fold and lig1_well; total_counts is the whole-transcriptome depth.
    """
    import anndata as ad

    adata = ad.read_h5ad(path)
    obs = adata.obs
    keep = (obs.condition.astype(str).isin(conditions)
            & obs.target.astype(str).isin(set(targets) | {NTC})).to_numpy()
    meta = obs.loc[keep, ["condition", "target", "calib_fold", "lig1_well", "total_counts"]].copy()
    for col in ("condition", "target", "lig1_well"):
        meta[col] = meta[col].astype(str)
    position = pd.Series(np.arange(adata.n_vars), index=adata.var.index).reindex(list(gene_ids))
    if position.isna().any():
        raise ValueError("Panel genes missing from the expression matrix")
    counts = sp.csr_matrix(adata.X[np.flatnonzero(keep)])
    X = log_cp10k(counts[:, position.to_numpy(int)], meta.total_counts.to_numpy(float))
    return np.asarray(X.todense(), dtype=np.float32), meta


class CellGraphs:
    """Stored post-collider PDAGs of one condition (outputs/cscn_panel_networks/<condition>/)."""

    def __init__(self, directory, nodes):
        directory = Path(directory)
        self.cells = pd.read_csv(directory / "cells.csv.gz", keep_default_na=False,
                                 dtype={"error": str, "meek_error": str})
        with np.load(directory / "edges.npz") as z:
            self.edges = {k: z[k] for k in ("cell_index", "source", "destination", "directed")}
            stored = list(z["nodes"])
        if stored != list(nodes):
            raise ValueError("Stored graph nodes do not match the frozen panel")
        if self.cells.error.ne("").any():
            raise ValueError("Cell graphs with CI-test errors are not supported")
        self.nodes = list(nodes)
        self.config = json.loads((directory / "config.json").read_text(encoding="utf-8"))
        bounds = np.searchsorted(self.edges["cell_index"], self.cells.cell_index.to_numpy())
        self.start = dict(zip(self.cells.cell_index, bounds))
        self.stop = dict(zip(self.cells.cell_index, np.r_[bounds[1:], len(self.edges["cell_index"])]))

    def pdag(self, cell_index):
        a, b = self.start[cell_index], self.stop[cell_index]
        return decode_edges(self.edges["source"][a:b], self.edges["destination"][a:b],
                            self.edges["directed"][a:b], self.nodes)

    def n_edges(self, cell_index):
        return int(self.stop[cell_index] - self.start[cell_index])


def structure_graphs(directory, nodes):
    """Shared and cluster starting PDAGs written by scripts/build_sem_graphs.py."""
    with np.load(Path(directory) / "graphs.npz") as z:
        if list(z["nodes"]) != list(nodes):
            raise ValueError("Structure graph nodes do not match the frozen panel")
        names = sorted({key.rsplit("_", 1)[0] for key in z.files if key.endswith("_source")})
        return {name: decode_edges(z[f"{name}_source"], z[f"{name}_destination"], z[f"{name}_directed"], nodes)
                for name in names}


def condition_inputs(root, cfg, condition):
    """Everything the SEM scripts need for one condition (NTC values in stored-graph order)."""
    from .cscn_panel import load_ntc_expression

    panel = json.loads((root / cfg["inputs"]["panel"]).read_text(encoding="utf-8"))
    ids, nodes = panel_nodes(panel)
    values, meta = load_ntc_expression(root / cfg["inputs"]["expression"], ids, [condition])
    graphs = CellGraphs(root / cfg["inputs"]["cell_graphs"] / condition, nodes)
    if list(graphs.cells.cell_id) != list(meta.index):
        raise ValueError(f"{condition}: stored graph cells do not match the NTC expression rows")
    sdir = root / cfg["outputs_dir"] / "structure" / condition
    return dict(panel=panel, ids=ids, nodes=nodes, index={n: i for i, n in enumerate(nodes)},
                values=values.astype(np.float64), meta=meta, graphs=graphs,
                neighbours=np.load(sdir / "neighbours.npy"), labels=np.load(sdir / "clusters.npy"),
                pcs=np.load(sdir / "pcs.npy"), structure=structure_graphs(sdir, nodes))
