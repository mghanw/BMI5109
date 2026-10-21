"""Fold-specific intervention orientation of stored PDAGs and local IDA for the held-out target.

Seed evidence comes from training-fold targets only (target - NTC responses of the panel genes);
the held-out targets of a fold never enter it. Seeding uses icpdag.intervention_orientation
unchanged, followed by the official Meek rules and DAG completion. For the held-out target, the
undirected edges incident to it are resolved by enumerating the locally valid parent sets (the
local criterion of IDA, Maathuis et al. 2009) instead of the completion's node sort order.

Completion rule, identical for every graph-based model: the official completion is accepted when it
is a valid extension (same skeleton, fixed edges kept, same unshielded colliders, acyclic); when that
fails, it is accepted under the relaxed rule (same skeleton, every fixed seeded or Meek-oriented edge
kept, acyclic; colliders may change) and counted as a fallback completion. The seed compatibility
check uses the same rule. A choice is dropped only when even the relaxed rule fails.
"""

from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind

from .intervention_orientation import fixed_edges, seed_interventions, valid_extension

NTC = "NTC"
MAX_UNDIRECTED = 10


class LocalIDAError(RuntimeError):
    """The held-out target has more undirected incident edges than the enumeration limit."""


def training_evidence(values, meta, train_targets, target_nodes, condition):
    """Seed evidence for one (fold, condition) from training targets only.

    values: cells x panel nodes (DataFrame); meta: condition and target for the same rows, which
    must be NTC or training-target cells only (raises otherwise). target_nodes: target -> node.
    Same estimates as icpdag.responses.response_table: mean difference and Welch t-test per panel
    gene (destination != source), BH within this (fold, condition), knockdown pass = negative change
    of the target's own node, eligible seed = knockdown pass and q < 0.05.
    """
    allowed = set(train_targets) | {NTC}
    unexpected = set(meta.target) - allowed
    if unexpected:
        raise ValueError(f"Evidence input contains non-training targets: {sorted(unexpected)[:5]}")
    if not values.index.equals(meta.index):
        raise ValueError("Expression and metadata rows must match")
    in_condition = meta.condition.eq(condition).to_numpy()
    columns = list(values.columns)
    position = {c: j for j, c in enumerate(columns)}
    # Contiguous float64 copies make every statistic independent of the frame's memory layout.
    data = np.ascontiguousarray(values.to_numpy(dtype=np.float64))
    control = np.ascontiguousarray(data[in_condition & meta.target.eq(NTC).to_numpy()])
    if len(control) < 2:
        raise ValueError(f"Fewer than two NTC cells in {condition}")
    control_mean = control.mean(axis=0)
    scale = np.clip(control.std(axis=0, ddof=1), 1e-8, None)
    rows = []
    for target in sorted(train_targets):
        treated = np.ascontiguousarray(data[in_condition & meta.target.eq(target).to_numpy()])
        if len(treated) < 2:
            raise ValueError(f"Fewer than two cells for {target} in {condition}")
        source = target_nodes[target]
        with np.errstate(invalid="ignore", divide="ignore"):
            _, pvalues = ttest_ind(treated, control, equal_var=False, axis=0)
        estimate = treated.mean(axis=0) - control_mean
        own = float(estimate[position[source]])
        for j, destination in enumerate(columns):
            if destination == source:
                continue
            rows.append(dict(condition=condition, target=target, source=source, destination=destination,
                             n_target=len(treated), n_control=len(control), estimate=float(estimate[j]),
                             standardized_estimate=float(estimate[j] / scale[j]),
                             p_value=float(pvalues[j]), test_valid=bool(np.isfinite(pvalues[j])),
                             target_log_expression_change=own, knockdown_pass=own < 0))
    frame = pd.DataFrame(rows)
    frame["q_value"] = false_discovery_control(frame.p_value.fillna(1).to_numpy(), method="bh")
    frame.loc[~frame.test_valid, "q_value"] = np.nan
    frame["eligible_seed"] = frame.test_valid & frame.knockdown_pass & frame.q_value.lt(0.05)
    return frame


def fold_seeds(values, meta, train_targets, target_nodes, condition):
    """Seed proposals for one (fold, condition): selects NTC and training-target rows only."""
    rows = meta.condition.eq(condition) & meta.target.isin(set(train_targets) | {NTC})
    evidence = training_evidence(values[rows], meta[rows], train_targets, target_nodes, condition)
    return seed_records(evidence), evidence


def seed_records(evidence):
    """Eligible seed proposals in the format of icpdag.responses.evidence_for_condition."""
    rows = evidence[evidence.eligible_seed]
    return rows[["source", "destination", "estimate", "standardized_estimate", "q_value"]].to_dict("records")


def relaxed_extension(pdag, dag):
    """Same nodes and skeleton, every fixed edge kept, acyclic; colliders may change."""
    pairs = lambda graph: {frozenset((a, b)) for a, b in graph.edges}
    return (set(pdag) == set(dag) and pairs(pdag) == pairs(dag) and fixed_edges(pdag) <= set(dag.edges)
            and nx.is_directed_acyclic_graph(dag))


def complete(pdag, model):
    """Official DAG completion: (dag, valid, fallback); dag is None when not accepted.

    fallback marks a completion accepted only under the relaxed rule (colliders changed).
    """
    try:
        dag = model._build_dag_representation(pdag)
    except (ValueError, RuntimeError):
        return None, False, False
    if valid_extension(pdag, dag):
        return dag, True, False
    if relaxed_extension(pdag, dag):
        return dag, True, True
    return None, False, False


def orient(pdag, seeds, model):
    """Seed with training evidence, propagate with Meek, complete with the official method."""
    seeded, records, absent = seed_interventions(pdag, seeds, model._build_dag_representation, relaxed_extension)
    meek = model._apply_meek_rules(seeded.copy(), apply_r4=False)
    dag, valid, fallback = complete(meek, model)
    return dict(pdag=meek, dag=dag, valid=valid, fallback=fallback, records=records, absent=absent)


def undirected_neighbours(pdag, node):
    return sorted(n for n in pdag.successors(node) if pdag.has_edge(n, node))


def directed_parents(pdag, node):
    return sorted(n for n in pdag.predecessors(node) if not pdag.has_edge(node, n))


def adjacent(pdag, a, b):
    return pdag.has_edge(a, b) or pdag.has_edge(b, a)


def local_parent_sets(pdag, node):
    """Subsets S of the undirected neighbours that create no new unshielded collider at node."""
    undirected = undirected_neighbours(pdag, node)
    parents = directed_parents(pdag, node)
    valid = []
    for size in range(len(undirected) + 1):
        for subset in combinations(undirected, size):
            if all(adjacent(pdag, a, b) for a, b in combinations(subset, 2)) and \
                    all(adjacent(pdag, s, p) for s in subset for p in parents):
                valid.append(subset)
    return undirected, valid


def local_ida_choices(pdag, node, model, max_undirected=MAX_UNDIRECTED):
    """DAGs for every locally valid orientation of node's undirected edges.

    For each parent set S: S -> node, node -> (U minus S), Meek, official completion; choices whose
    completion is not accepted (see complete) are dropped and counted.
    Returns dict(dags, n_undirected, n_candidates, n_failed, n_fallback).
    """
    undirected, subsets = local_parent_sets(pdag, node)
    if len(undirected) > max_undirected:
        raise LocalIDAError(f"{node}: {len(undirected)} undirected incident edges (limit {max_undirected})")
    if not undirected:
        dag, valid, fallback = complete(pdag, model)
        return dict(dags=[dag] if valid else [], n_undirected=0, n_candidates=1, n_failed=int(not valid),
                    n_fallback=int(fallback))
    dags, failed, fallbacks = [], 0, 0
    for subset in subsets:
        graph = pdag.copy()
        for u in undirected:
            if u in subset:
                graph.remove_edge(node, u)
            else:
                graph.remove_edge(u, node)
        graph = model._apply_meek_rules(graph, apply_r4=False)
        dag, valid, fallback = complete(graph, model)
        if valid:
            dags.append(dag)
            fallbacks += int(fallback)
        else:
            failed += 1
    return dict(dags=dags, n_undirected=len(undirected), n_candidates=len(subsets), n_failed=failed,
                n_fallback=fallbacks)


def undirected_pdag(nodes, pairs):
    """PDAG with every listed pair undirected (association skeletons have no collider step)."""
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    for a, b in pairs:
        graph.add_edge(a, b)
        graph.add_edge(b, a)
    return graph
