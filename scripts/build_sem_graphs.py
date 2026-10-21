"""Structure for the SEM comparators: PCs, neighbourhoods, clusters, shared and cluster graphs (NTC only).

Usage:
    python scripts/build_sem_graphs.py [--config config/sem_prediction.json] [--conditions DMSO,Vem]

Per condition, from the condition's NTC cells only:
  - 20 PCs of panel log1p CP10k, neighbourhoods of round(0.35 n) nearest cells, K = 6 k-means clusters;
  - the shared graph: population PC, target-coverage matched (the sparsest graph in which at least as
    many of the 56 targets have an incident edge as have an edge in >= 10% of the cells' CSCN graphs);
  - one graph per cluster, coverage-matched to the same count over the cluster's cells.
Writes outputs/sem_prediction/structure/<condition>/ and results/summary/sem_prediction_shared_graphs.csv
(target connectivity of every graph next to the cell-specific numbers). No response data is read.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from aim1 import assert_no_cell_ids  # noqa: E402
from icpdag.cscn_networks import encode_edges  # noqa: E402
from icpdag.cscn_panel import load_ntc_expression  # noqa: E402
from icpdag.sem_graphs import clusters, coverage_matched_graph, nearest_neighbours, principal_components  # noqa: E402
from icpdag.sem_inputs import CellGraphs, panel_nodes, target_nodes  # noqa: E402


def target_presence(graphs, target_list):
    """cells x targets: True where the target has an edge in that cell's stored CSCN graph."""
    position = {n: i for i, n in enumerate(graphs.nodes)}
    columns = {position[t]: j for j, t in enumerate(target_list)}
    row = {c: i for i, c in enumerate(graphs.cells.cell_index)}
    present = np.zeros((len(row), len(target_list)), bool)
    for end in ("source", "destination"):
        for cell, node in zip(graphs.edges["cell_index"], graphs.edges[end]):
            if node in columns:
                present[row[cell], columns[node]] = True
    return present


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "sem_prediction.json"))
    parser.add_argument("--conditions")
    args = parser.parse_args(argv)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    panel = json.loads((ROOT / cfg["inputs"]["panel"]).read_text(encoding="utf-8"))
    networks = json.loads((ROOT / cfg["inputs"]["networks"]).read_text(encoding="utf-8"))
    ids, nodes = panel_nodes(panel)
    tnodes = target_nodes(panel)
    conditions = args.conditions.split(",") if args.conditions else cfg["conditions"]
    pcs_cfg, k = cfg["graphs"]["pcs"]["n_components"], cfg["graphs"]["cluster"]["k"]
    connectivity = pd.read_csv(ROOT / cfg["inputs"]["target_connectivity"])
    connected = set(cfg["weak_targets"]["connected_targets"])
    rows = []
    for condition in conditions:
        t0 = time.time()
        values, meta = load_ntc_expression(ROOT / cfg["inputs"]["expression"], ids, [condition])
        graphs = CellGraphs(ROOT / cfg["inputs"]["cell_graphs"] / condition, nodes)
        if list(graphs.cells.cell_id) != list(meta.index):
            raise ValueError(f"{condition}: stored graph cells do not match the NTC expression rows")
        n = len(values)
        scores = principal_components(values, pcs_cfg, cfg["seed"])
        neighbours = nearest_neighbours(scores, int(round(0.35 * n)))
        labels = clusters(scores, k, cfg["seed"])
        edges_per_cell = graphs.cells.n_edges.to_numpy()
        target_list = sorted(tnodes.values())
        presence = target_presence(graphs, target_list)
        out = ROOT / cfg["outputs_dir"] / "structure" / condition
        out.mkdir(parents=True, exist_ok=True)
        np.save(out / "pcs.npy", scores.astype(np.float32))
        np.save(out / "neighbours.npy", neighbours)
        np.save(out / "clusters.npy", labels.astype(np.int16))
        info, arrays = {}, {}
        jobs = [("shared", np.arange(n))] + [(f"cluster_{c + 1}", np.flatnonzero(labels == c)) for c in range(k)]
        for name, rows_used in jobs:
            median_edges = int(np.rint(np.median(edges_per_cell[rows_used])))
            required = int((presence[rows_used].mean(axis=0) >= 0.10).sum())
            g = coverage_matched_graph(values[rows_used], nodes, target_list, required,
                                       networks["max_conditioning"], max_seconds=600)
            src, dst, directed = encode_edges(g["pdag"], nodes)
            arrays.update({f"{name}_source": src, f"{name}_destination": dst, f"{name}_directed": directed})
            info[name] = {key: g[key] for key in ("z_crit", "log10_alpha", "edges", "targets_with_edge",
                                                  "min_targets", "n_obs", "zero_variance", "iterations",
                                                  "level0_edges")}
            info[name]["median_cell_specific_edges"] = median_edges
            info[name]["undirected_share"] = float((~directed).mean()) if len(directed) else float("nan")
            rows.append(dict(condition=condition, graph=name, n_cells=len(rows_used),
                             required_targets_with_edge=required, targets_with_edge=g["targets_with_edge"],
                             n_targets=len(tnodes), edges=g["edges"], median_cell_specific_edges=median_edges,
                             z_crit=g["z_crit"], log10_alpha=g["log10_alpha"], level0_edges=g["level0_edges"],
                             undirected_share=info[name]["undirected_share"],
                             n_zero_variance_genes=len(g["zero_variance"])))
            print(f"{condition} {name}: n={len(rows_used)} coverage {g['targets_with_edge']}/56 (required "
                  f"{required}); edges {g['edges']} (median cell-specific {median_edges}); z_crit {g['z_crit']:.2f} "
                  f"log10 alpha {g['log10_alpha']:.1f} ({time.time() - t0:.0f}s)")
        np.savez_compressed(out / "graphs.npz", nodes=np.asarray(nodes, dtype=str), **arrays)
        (out / "graphs.json").write_text(json.dumps(dict(info, cluster_sizes=np.bincount(labels).tolist(),
                                                         neighbourhood_size=int(neighbours.shape[1])), indent=2),
                                         encoding="utf-8")
        share = connectivity[connectivity.condition.eq(condition)].set_index("target").share_cells_with_edge
        if int((share >= 0.10).sum()) != int((presence.mean(axis=0) >= 0.10).sum()):
            raise RuntimeError("Target coverage does not match cscn_panel_target_connectivity.csv")
        for name, value in [("cell_specific_edge_in_any_cell", int((share > 0).sum())),
                            ("cell_specific_edge_in_10pct_of_cells_this_condition", int((share >= 0.10).sum())),
                            ("cell_specific_connected_stratum_both_conditions", len(connected))]:
            rows.append(dict(condition=condition, graph=name, n_cells=n, targets_with_edge=value,
                             n_targets=len(tnodes), median_cell_specific_edges=int(np.rint(np.median(edges_per_cell)))))
    table = pd.DataFrame(rows)
    assert_no_cell_ids(table)
    path = ROOT / "results" / "summary" / "sem_prediction_shared_graphs.csv"
    table.to_csv(path, index=False, lineterminator="\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
