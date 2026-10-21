"""Fold-specific orientation statistics before any coefficient or prediction (structure + training seeds).

Usage:
    python scripts/sem_structure_report.py [--config config/sem_prediction.json] [--workers W]
                                           [--conditions DMSO,Vem] [--cells N]
                                           [--graph-types cell_specific,shared,cluster,local_association]

For every outer fold and condition: training-target seed evidence (held-out targets excluded), then
for each graph type (cell_specific, shared, cluster, local_association) the seeded + Meek PDAGs and
the number of undirected edges |U| incident to each held-out target, per (reference cell, held-out
target). Needs scripts/build_sem_graphs.py first. Writes results/summary/sem_prediction_structure.csv,
sem_prediction_structure_targets.csv and sem_prediction_orientation.csv (with --cells:
outputs/sem_prediction_structure_smoke/ only). Each (condition, graph type) result is cached in
outputs/sem_prediction/structure/<condition>/ with its rule; graph types not listed in --graph-types
are read from that cache (the rule and cell set must match).
"""
import argparse
from collections import Counter
import json
import multiprocessing as mp
import pickle
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from aim1 import assert_no_cell_ids  # noqa: E402
from icpdag.cscn_networks import decode_edges, reference_model  # noqa: E402
from icpdag.cscn_panel import load_ntc_expression  # noqa: E402
from icpdag.sem_graphs import association_pairs  # noqa: E402
from icpdag.sem_inputs import CellGraphs, fold_targets, load_panel_cells, panel_nodes, target_nodes  # noqa: E402
from icpdag.sem_orientation import (MAX_UNDIRECTED, local_ida_choices, orient, seed_records,  # noqa: E402
                                    training_evidence, undirected_neighbours, undirected_pdag)
from icpdag.target_split import load_split  # noqa: E402

STATUSES = ["accepted_seed", "agrees_collider", "collider_conflict", "bidirectional_response_conflict",
            "cycle_conflict", "new_collider_conflict", "no_consistent_extension"]
RULES = {"cell_specific": "stored CSCN PDAG; unified fallback completion; |U| <= 10",
         "shared": "coverage-matched population PC; unified fallback completion; |U| <= 10",
         "cluster": "coverage-matched population PC per cluster; unified fallback completion; |U| <= 10",
         "local_association": "density-matched association skeleton; unified fallback completion; |U| <= 10"}
# Earlier rules whose results are identical under the current rule when no unit had |U| > 8: the
# association model already used the relaxed completion for seeding, fold DAGs and IDA choices.
EQUIVALENT = {"local_association": "density-matched association skeleton; relaxed completion"}
_W = {}


def _init(state):
    _W.update(state)
    _W["model"] = reference_model(state["values"], state["nodes"], state["cscn"])


def orient_graph(pdag, seeds, held, model):
    """Orientation result summary for one starting PDAG and one fold.

    Per held-out target: |U|, whether it is isolated, whether at least one local IDA choice has a
    valid completion, and whether any valid choice gives it a descendant (else its prediction is
    constant). Structure only; no coefficients.
    """
    result = orient(pdag, seeds, model)
    statuses = Counter(r["status"] for r in result["records"])
    us, isolated, valid, downstream = [], [], [], []
    choices = fallbacks = 0
    for g in held:
        u = len(undirected_neighbours(result["pdag"], g))
        us.append(u)
        isolated.append(pdag.degree(g) == 0)
        if u > MAX_UNDIRECTED:
            valid.append(False)
            downstream.append(False)
            continue
        if u:
            ida = local_ida_choices(result["pdag"], g, model)
            dags = ida["dags"]
            fallbacks += ida["n_fallback"]
        else:
            dags = [result["dag"]] if result["valid"] else []
            fallbacks += int(result["fallback"])
        choices += len(dags)
        valid.append(bool(dags))
        downstream.append(any(nx.descendants(d, g) for d in dags))
    return dict(valid=result["valid"], fallback=result["fallback"], statuses=statuses, absent=result["absent"],
                us=us, isolated=isolated, choice_valid=valid, downstream=downstream, held=list(held),
                choices=choices, choice_fallbacks=fallbacks)


def _run(task):
    nodes, values, neighbours, graphs = _W["nodes"], _W["values"], _W["neighbours"], _W["graphs"]
    tally = Counter()
    for k in task["cells"]:
        if task["kind"] == "cell_specific":
            start = graphs.pdag(k)
        else:
            pairs = association_pairs(values[neighbours[k]], graphs.n_edges(k))
            start = undirected_pdag(nodes, [(nodes[i], nodes[j]) for i, j in pairs])
        for fold, seeds in _W["seeds"].items():
            out = orient_graph(start, seeds, _W["held"][fold], _W["model"])
            add(tally, fold, out, 1)
    return task["kind"], tally


def add(tally, fold, out, weight):
    tally[(fold, "graphs")] += weight
    tally[(fold, "valid")] += weight * int(out["valid"])
    tally[(fold, "fallback")] += weight * int(out["fallback"])
    tally[(fold, "choices")] += weight * out["choices"]
    tally[(fold, "choice_fallbacks")] += weight * out["choice_fallbacks"]
    tally[(fold, "absent")] += weight * out["absent"]
    for status, count in out["statuses"].items():
        tally[(fold, status)] += weight * count
    for g, u, isolated, valid, downstream in zip(out["held"], out["us"], out["isolated"], out["choice_valid"],
                                                 out["downstream"]):
        tally[(fold, "U", min(u, MAX_UNDIRECTED + 1))] += weight
        tally[(fold, "U_sum")] += weight * u
        tally[(fold, "U_max")] = max(tally[(fold, "U_max")], u)
        tally[(fold, "isolated")] += weight * int(isolated)
        tally[(fold, "no_valid_choice")] += weight * int(not valid)
        tally[(fold, "downstream")] += weight * int(downstream)
        tally[("target", g, "units")] += weight
        tally[("target", g, "valid")] += weight * int(valid)
        tally[("target", g, "downstream")] += weight * int(downstream)


def tables(condition, kind, tally, held):
    structure, orientation = [], []
    for fold in sorted(held):
        n_units = sum(tally[(fold, "U", u)] for u in range(MAX_UNDIRECTED + 2))
        row = dict(graph_type=kind, condition=condition, fold=fold, n_held_out_targets=len(held[fold]),
                   n_units=n_units, mean_U=tally[(fold, "U_sum")] / n_units if n_units else np.nan,
                   max_U=tally[(fold, "U_max")], share_target_isolated=tally[(fold, "isolated")] / n_units,
                   share_fold_dag_valid=tally[(fold, "valid")] / tally[(fold, "graphs")],
                   share_fold_dag_fallback=tally[(fold, "fallback")] / tally[(fold, "graphs")],
                   share_choices_fallback=(tally[(fold, "choice_fallbacks")] / tally[(fold, "choices")]
                                           if tally[(fold, "choices")] else np.nan),
                   share_units_no_valid_choice=tally[(fold, "no_valid_choice")] / n_units,
                   share_units_with_downstream=tally[(fold, "downstream")] / n_units)
        row.update({f"U_{u}": tally[(fold, "U", u)] for u in range(MAX_UNDIRECTED + 1)})
        row[f"U_gt_{MAX_UNDIRECTED}"] = tally[(fold, "U", MAX_UNDIRECTED + 1)]
        structure.append(row)
        o = dict(graph_type=kind, condition=condition, fold=fold, n_graphs=tally[(fold, "graphs")],
                 seed_proposals_absent=tally[(fold, "absent")])
        o.update({s: tally[(fold, s)] for s in STATUSES})
        orientation.append(o)
    targets = []
    for fold in sorted(held):
        for g in held[fold]:
            units = tally[("target", g, "units")]
            targets.append(dict(graph_type=kind, condition=condition, fold=fold, node=g, n_cells=units,
                                share_cells_valid_choice=tally[("target", g, "valid")] / units,
                                share_cells_with_downstream=tally[("target", g, "downstream")] / units,
                                constant_prediction_expected=tally[("target", g, "downstream")] == 0,
                                exclusion_flag=1 - tally[("target", g, "valid")] / units > 0.05))
    return structure, orientation, targets


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "sem_prediction.json"))
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--conditions")
    parser.add_argument("--cells", type=int, help="Smoke test: seeded subset of reference cells")
    parser.add_argument("--graph-types", default=",".join(RULES), help="Graph types to recompute")
    args = parser.parse_args(argv)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    panel = json.loads((ROOT / cfg["inputs"]["panel"]).read_text(encoding="utf-8"))
    networks = json.loads((ROOT / cfg["inputs"]["networks"]).read_text(encoding="utf-8"))
    cscn = dict(window=networks["window"], ci_alpha=networks["ci_alpha"], max_conditioning=networks["max_conditioning"])
    ids, nodes = panel_nodes(panel)
    tnodes = target_nodes(panel)
    split = load_split(ROOT / cfg["inputs"]["target_split"])
    folds = sorted(split.fold.unique())
    conditions = args.conditions.split(",") if args.conditions else cfg["conditions"]
    all_values, all_meta = load_panel_cells(ROOT / cfg["inputs"]["expression"], ids, conditions, split.target)
    all_values = pd.DataFrame(all_values, index=all_meta.index, columns=nodes)
    structure, orientation, seed_counts, per_target = [], [], [], []
    for condition in conditions:
        t0 = time.time()
        values, meta = load_ntc_expression(ROOT / cfg["inputs"]["expression"], ids, [condition])
        graphs = CellGraphs(ROOT / cfg["inputs"]["cell_graphs"] / condition, nodes)
        sdir = ROOT / cfg["outputs_dir"] / "structure" / condition
        neighbours = np.load(sdir / "neighbours.npy")
        labels = np.load(sdir / "clusters.npy")
        seeds, held = {}, {}
        for fold in folds:
            held_targets, train = fold_targets(split, fold)
            rows = all_meta.condition.eq(condition) & all_meta.target.isin(set(train) | {"NTC"})
            evidence = training_evidence(all_values[rows], all_meta[rows], train, tnodes, condition)
            evidence.assign(fold=fold).to_csv(sdir / f"evidence_fold{fold}.csv.gz", index=False)
            seeds[fold] = seed_records(evidence)
            held[fold] = [tnodes[t] for t in held_targets]
            seed_counts.append(dict(condition=condition, fold=fold, n_training_targets=len(train),
                                    n_tests=len(evidence), n_eligible_seeds=len(seeds[fold])))
        print(f"{condition}: evidence for {len(folds)} folds ({time.time() - t0:.0f}s); eligible seeds "
              f"{[len(s) for s in seeds.values()]}")
        cells = graphs.cells.cell_index.to_numpy()
        if args.cells:
            cells = np.sort(np.random.default_rng(cfg["seed"]).choice(cells, args.cells, replace=False))
        compute = [k for k in args.graph_types.split(",") if k]
        cache_dir = (ROOT / "outputs/sem_prediction_structure_smoke" / condition if args.cells
                     else ROOT / cfg["outputs_dir"] / "structure" / condition)
        cache_dir.mkdir(parents=True, exist_ok=True)
        tallies = {}
        for kind in RULES:
            if kind in compute:
                continue
            cached = pickle.loads((cache_dir / f"structure_tally_{kind}.pkl").read_bytes())
            max_u = max([v for k, v in cached["tally"].items() if k[-1] == "U_max"], default=0)
            equivalent = cached["rule"] == EQUIVALENT.get(kind) and max_u <= 8
            if (cached["rule"] != RULES[kind] and not equivalent) or cached["cells"] != cells.tolist():
                raise ValueError(f"Cached {kind} result for {condition} was made with another rule or cell set")
            print(f"{condition}: reusing cached {kind} result ({cached['rule']})")
            tallies[kind] = cached["tally"]
        state = dict(nodes=nodes, values=values, neighbours=neighbours, graphs=graphs, seeds=seeds, held=held,
                     cscn=cscn)
        cell_kinds = [k for k in ("cell_specific", "local_association") if k in compute]
        for kind in cell_kinds:
            tallies[kind] = Counter()
        tasks = [dict(kind=kind, cells=cells[i:i + 16].tolist()) for kind in cell_kinds
                 for i in range(0, len(cells), 16)]
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"),
                                 initializer=_init, initargs=(state,)) as pool:
            for kind, tally in pool.map(_run, tasks):
                tallies[kind].update({k: v for k, v in tally.items() if k[-1] != "U_max"})
                for k, v in tally.items():
                    if k[-1] == "U_max":
                        tallies[kind][k] = max(tallies[kind][k], v)
        print(f"{condition}: cell-level orientation done ({time.time() - t0:.0f}s)")
        _init(state)
        with np.load(sdir / "graphs.npz") as z:
            stored = {name: decode_edges(z[f"{name}_source"], z[f"{name}_destination"], z[f"{name}_directed"], nodes)
                      for name in ["shared"] + [f"cluster_{c + 1}" for c in range(labels.max() + 1)]}
        in_subset = np.isin(np.arange(len(values)), cells)
        for kind, members in [("shared", {"shared": in_subset.sum()})] + [
                ("cluster", {f"cluster_{c + 1}": int((labels[in_subset] == c).sum()) for c in range(labels.max() + 1)})]:
            if kind not in compute:
                continue
            tally = Counter()
            for name, weight in members.items():
                for fold in folds:
                    add(tally, fold, orient_graph(stored[name], seeds[fold], held[fold], _W["model"]), weight)
            tallies[kind] = tally
        for kind in compute:
            (cache_dir / f"structure_tally_{kind}.pkl").write_bytes(pickle.dumps(
                dict(rule=RULES[kind], cells=cells.tolist(), tally=tallies[kind])))
        for kind in RULES:
            tally = tallies[kind]
            s, o, tg = tables(condition, kind, tally, held)
            structure += s
            orientation += o
            per_target += tg
        print(f"{condition}: done ({time.time() - t0:.0f}s)")
    structure, orientation = pd.DataFrame(structure), pd.DataFrame(orientation)
    orientation = orientation.merge(pd.DataFrame(seed_counts), on=["condition", "fold"])
    node_target = {v: k for k, v in tnodes.items()}
    per_target = pd.DataFrame(per_target)
    per_target.insert(3, "target", per_target.node.map(node_target))
    out = ROOT / ("results/summary" if not args.cells else "outputs/sem_prediction_structure_smoke")
    out.mkdir(parents=True, exist_ok=True)
    for frame, name in [(structure, "sem_prediction_structure.csv"), (orientation, "sem_prediction_orientation.csv"),
                        (per_target, "sem_prediction_structure_targets.csv")]:
        assert_no_cell_ids(frame)
        frame.to_csv(out / name, index=False, lineterminator="\n")
        print(f"wrote {out / name}")
    over = structure[f"U_gt_{MAX_UNDIRECTED}"].sum()
    print(f"(cell, held-out target) units with |U| > {MAX_UNDIRECTED}: {over}")


if __name__ == "__main__":
    main()
