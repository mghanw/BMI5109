"""Aggregate summaries of the cell-specific CSCN graphs (committed tables; no cell IDs).

Usage:
    python scripts/summarize_cscn_networks.py [--config config/cscn_networks.json] [--run-name cscn_panel_networks]
                                              [--no-figure]

Reads the merged per-condition outputs of scripts/infer_cscn_networks.py
(outputs/<run>/<condition>/cells.csv.gz, edges.npz) and writes to results/summary/:
  cscn_panel_graph_summary.csv        per condition: cells, errors, edges per cell, undirected share,
                                      CI tests, unseeded DAG-completion validity, runtime
  cscn_panel_edge_frequency.csv       per condition and gene pair: share of cells with the adjacency
                                      (and each orientation); pairs in >= 1% of cells
  cscn_panel_target_connectivity.csv  per target and condition: share of cells with >= 1 edge, degree,
                                      panel genes within 2 skeleton steps, isolated-majority flag
  cscn_panel_summary.md
and the figure reports/figures/cscn_panel_graph_overview.{png,svg}. Cells with a CI-test error are
excluded from the graph statistics and counted separately.
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
from icpdag.target_split import add_fold, load_split  # noqa: E402

MIN_EDGE_SHARE = 0.01
WEAK_EDGE_SHARE = 0.01


def load_condition(directory):
    cells = pd.read_csv(directory / "cells.csv.gz", keep_default_na=False, dtype={"error": str, "meek_error": str})
    with np.load(directory / "edges.npz") as z:
        edges = {k: z[k] for k in ("cell_index", "source", "destination", "directed")}
        nodes = list(z["nodes"])
    cells["meek_valid"] = cells.meek_valid.astype(str).eq("True")
    return cells, edges, nodes


def graph_summary(condition, cells, edges, n_targets):
    ok = cells[cells.error.eq("")]
    keep = np.isin(edges["cell_index"], ok.cell_index)
    targets = per_cell_target_presence(edges, ok.cell_index.to_numpy(), n_targets)
    q = ok.n_edges.quantile([0.25, 0.5, 0.75])
    return dict(condition=condition, n_cells=len(cells), n_ci_errors=len(cells) - len(ok),
                edges_per_cell_median=q[0.5], edges_per_cell_q1=q[0.25], edges_per_cell_q3=q[0.75],
                edges_per_cell_mean=ok.n_edges.mean(), cells_without_edges=int(ok.n_edges.eq(0).sum()),
                undirected_share_pooled=float((~edges["directed"][keep]).mean()) if keep.any() else np.nan,
                undirected_share_median=(ok.n_undirected / ok.n_edges.where(ok.n_edges > 0)).median(),
                target_nodes_with_edge_mean=targets.mean(),
                ci_tests_per_cell_median=ok.n_ci_tests.median(),
                unseeded_dag_completion_valid_rate=ok.meek_valid.mean(),
                seconds_per_cell_mean=cells.seconds.mean(), seconds_per_cell_median=cells.seconds.median(),
                cpu_hours=cells.seconds.sum() / 3600)


def per_cell_target_presence(edges, cell_index, n_targets):
    """Fraction of target nodes (indices 0..n_targets-1) with >= 1 edge, per cell."""
    position = {c: i for i, c in enumerate(cell_index)}
    present = np.zeros((len(cell_index), n_targets), bool)
    for end in ("source", "destination"):
        mask = edges[end] < n_targets
        rows = [position.get(c, -1) for c in edges["cell_index"][mask]]
        rows = np.asarray(rows, int)
        good = rows >= 0
        present[rows[good], edges[end][mask][good]] = True
    return present.mean(axis=1)


def edge_frequency(condition, cells, edges, nodes):
    ok = cells.cell_index[cells.error.eq("")].to_numpy()
    keep = np.isin(edges["cell_index"], ok)
    src, dst = (edges[k][keep].astype(np.int64) for k in ("source", "destination"))
    directed = edges["directed"][keep].astype(bool)
    a, b = np.minimum(src, dst), np.maximum(src, dst)
    key = a * len(nodes) + b
    forward = directed & (src == a)
    backward = directed & (src == b)
    frame = pd.DataFrame({"key": key, "fwd": forward, "bwd": backward, "und": ~directed})
    agg = frame.groupby("key").agg(n=("key", "size"), fwd=("fwd", "sum"), bwd=("bwd", "sum"), und=("und", "sum"))
    n = len(ok)
    agg = agg[agg.n / n >= MIN_EDGE_SHARE]
    i, j = np.divmod(agg.index.to_numpy(), len(nodes))
    out = pd.DataFrame({"condition": condition, "gene_a": np.asarray(nodes)[i], "gene_b": np.asarray(nodes)[j],
                        "n_cells": n, "share_cells_adjacent": agg.n.to_numpy() / n,
                        "share_cells_a_to_b": agg.fwd.to_numpy() / n, "share_cells_b_to_a": agg.bwd.to_numpy() / n,
                        "share_cells_undirected": agg.und.to_numpy() / n})
    return out.sort_values(["share_cells_adjacent", "gene_a", "gene_b"], ascending=[False, True, True])


def target_connectivity(condition, cells, edges, nodes, targets):
    """targets: list of (target name, node index)."""
    ok = cells.cell_index[cells.error.eq("")].to_numpy()
    keep = np.isin(edges["cell_index"], ok)
    frame = pd.DataFrame({k: edges[k][keep] for k in ("cell_index", "source", "destination")})
    degree = {t: [] for t, _ in targets}
    two_step = {t: [] for t, _ in targets}
    groups = dict(tuple(frame.groupby("cell_index")))
    for cell in ok:
        g = groups.get(cell)
        neighbours = {}
        if g is not None:
            for s, d in zip(g.source.to_numpy(), g.destination.to_numpy()):
                neighbours.setdefault(int(s), set()).add(int(d))
                neighbours.setdefault(int(d), set()).add(int(s))
        for t, node in targets:
            first = neighbours.get(node, set())
            reach = set(first)
            for n1 in first:
                reach |= neighbours.get(n1, set())
            reach.discard(node)
            degree[t].append(len(first))
            two_step[t].append(len(reach))
    rows = []
    for t, node in targets:
        d = np.asarray(degree[t])
        rows.append(dict(target=t, symbol=nodes[node], condition=condition, n_cells=len(ok),
                         share_cells_with_edge=float((d > 0).mean()), median_degree=float(np.median(d)),
                         mean_degree=float(d.mean()), median_genes_within_2_steps=float(np.median(two_step[t])),
                         isolated_majority=bool((d == 0).mean() > 0.5)))
    return pd.DataFrame(rows)


def figure(edge_counts, connectivity):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from aim1_figures import COLORS, GRID, MARKERS, MUTED, save, style

    style()
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 150 / 25.4), gridspec_kw=dict(width_ratios=[1, 1.15]))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.08, wspace=0.42)
    ax = axes[0]
    top = max(c.max() for c in edge_counts.values())
    bins = np.arange(0, top + 6, 5)
    for condition, counts in edge_counts.items():
        ax.hist(counts, bins=bins, histtype="step", lw=2, color=COLORS[condition],
                label=f"{condition} (median {np.median(counts):.0f})")
    ax.set_xlabel("Edges per cell (skeleton adjacencies)")
    ax.set_ylabel("NTC cells")
    ax.set_title("Edges per cell", loc="left")
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.legend(loc="upper right")
    ax.text(-0.16, 1.035, "a", transform=ax.transAxes, fontsize=9, weight="bold")
    ax = axes[1]
    wide = connectivity.pivot(index="target", columns="condition", values="share_cells_with_edge")
    wide = wide.loc[wide.mean(axis=1).sort_values().index]
    y = np.arange(len(wide))
    ax.hlines(y, wide.min(axis=1), wide.max(axis=1), color=GRID, lw=1, zorder=1)
    for condition in wide.columns:
        ax.scatter(wide[condition], y, s=14, marker=MARKERS[condition], color=COLORS[condition], zorder=3,
                   edgecolors="white", linewidths=0.5)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_yticks(y, wide.index, fontsize=5.5)
    ax.set_ylim(-0.7, len(wide) - 0.3)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Share of NTC cells where the target has >= 1 edge")
    ax.set_title("Target connectivity", loc="left")
    ax.text(0.52, 8, "<- isolated in a majority of cells", ha="left", va="center", fontsize=6, color=MUTED)
    ax.legend(handles=[Line2D([], [], marker=MARKERS[c], ls="", color=COLORS[c], label=c) for c in wide.columns],
              loc="lower right")
    ax.text(-0.25, 1.035, "b", transform=ax.transAxes, fontsize=9, weight="bold")
    fig.suptitle("Cell-specific CSCN graphs on the 250-gene network panel (NTC cells)", x=0.08, y=0.975, ha="left",
                 fontsize=9)
    save(fig, "cscn_panel_graph_overview")


def markdown(summary, connectivity, freq, panel, cfg, timing, sensitivity):
    s = summary.set_index("condition")
    both = connectivity.pivot(index="target", columns="condition", values="isolated_majority").all(axis=1)
    n_targets = connectivity.target.nunique()
    any_cond = connectivity.groupby("condition").isolated_majority.sum()
    share = connectivity.pivot(index="target", columns="condition", values="share_cells_with_edge")
    never = {c: list(share.index[share[c].eq(0)]) for c in share.columns}
    undefined = share.index[share.eq(0).any(axis=1)]
    weak = share[share.lt(WEAK_EDGE_SHARE).any(axis=1)]
    weak_text = ", ".join(f"{t} ({' / '.join(f'{c} {100 * v:.1f}%' for c, v in row.items())})"
                          for t, row in weak.iterrows())
    sel = cfg["selection"]
    lines = [
        "# Cell-specific CSCN graphs on the network panel",
        "",
        "Structure is learned from NTC cells only (both conditions), so the graphs are fold-independent and one",
        "graph set serves all five held-out folds. Stored graph per cell: the PDAG after unshielded-collider",
        "orientation; intervention-based orientation is fold-specific and done when SEMs are fitted.",
        "",
        "## Panel",
        f"- {panel['panel_size']} genes: {panel['n_forced']} forced (56 eligible targets + Fufa program genes; PTEN is",
        f"  both) and {panel['panel_size'] - panel['n_forced']} fill genes from a universe of {panel['n_fill_universe']:,}"
        " genes detected in >= 10% of NTC cells in each condition.",
        "- Fill ranking: Seurat-style mean.var.plot standardization, computed on log1p CP10k with equal-frequency",
        "  mean bins (mean of the two within-condition ranks).",
        f"- {sum(g['low_detection'] for g in panel['forced_genes'])} forced genes (all Fufa) have NTC detection < 5% "
        "in both conditions; they are kept and flagged (cscn_panel_genes.csv).",
        "",
        "## Hyperparameters (graph diagnostics and runtime only, never prediction performance)",
        f"- window {cfg['window']}, CI alpha {cfg['ci_alpha']}, max conditioning {cfg['max_conditioning']}, stable PC, "
        "bitmap queries (official CSCN, pinned).",
        f"- {sel['rationale']}",
        f"- Grids: cscn_panel_timing.csv ({len(timing)} configurations of panel size x max conditioning) and",
        f"  cscn_panel_sensitivity.csv ({len(sensitivity)} configurations of window x CI alpha at P 250).",
        "",
        "## Graphs",
        "| condition | cells | CI errors | edges per cell, median [IQR] | undirected share | target nodes with an edge | "
        "unseeded DAG completion valid | CPU hours |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c, r in s.iterrows():
        lines.append(f"| {c} | {int(r.n_cells):,} | {int(r.n_ci_errors)} | {r.edges_per_cell_median:.0f} "
                     f"[{r.edges_per_cell_q1:.0f}, {r.edges_per_cell_q3:.0f}] | {r.undirected_share_pooled:.2f} | "
                     f"{r.target_nodes_with_edge_mean:.2f} | {r.unseeded_dag_completion_valid_rate:.3f} | "
                     f"{r.cpu_hours:.1f} |")
    lines += [
        "",
        f"- Edge frequency: {len(freq):,} (condition, gene pair) rows present in >= 1% of cells "
        "(cscn_panel_edge_frequency.csv), for the shared and cluster-level comparators.",
        "",
        "## Target connectivity",
        f"- Targets isolated in a majority of cells (no edge in > 50% of cells): DMSO {any_cond.get('DMSO', 0)}, "
        f"Vem {any_cond.get('Vem', 0)} of {n_targets}; in both conditions {int(both.sum())} of {n_targets}.",
        "- Connectivity tracks NTC detection of the target: in a cell where the target is not detected, its",
        "  neighbourhood is the zero tie block and the local CI tests return independence.",
        "- Targets with no edge in any cell: "
        + "; ".join(f"{c}: {', '.join(v) if v else 'none'}" for c, v in never.items()) + ".",
        "- Only targets without an edge in any cell get constant downstream predictions (undefined Pearson r).",
        "  Under the pre-registered rule (NaN if either condition is undefined) the expected NaN share is "
        f"{len(undefined)}/{n_targets} ({100 * len(undefined) / n_targets:.1f}%),",
        "  below the 10% flag.",
        f"- Weakly connected (an edge in < {100 * WEAK_EDGE_SHARE:.0f}% of cells in at least one condition; "
        f"{len(weak)} targets): {weak_text}.",
        "  Weakly connected targets get predictions from very few cells.",
        f"- Isolated in a majority of cells in both conditions: {', '.join(both[both].index)}.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "cscn_networks.json"))
    parser.add_argument("--run-name", default="cscn_panel_networks")
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    with open(ROOT / cfg["paths"]["panel"], encoding="utf-8") as f:
        panel = json.load(f)
    if not panel.get("frozen"):
        raise ValueError("Summaries need the frozen panel")
    run = ROOT / cfg["paths"]["outputs_dir"] / args.run_name
    out = ROOT / cfg["paths"]["summary_dir"]
    split = load_split(ROOT / "config" / "target_split.json")
    target_nodes = [(g["target"], i) for i, g in enumerate(panel["forced_genes"]) if g["role"] == "target"]
    if [i for _, i in target_nodes] != list(range(len(target_nodes))):
        raise ValueError("Target nodes must come first in the panel")
    summaries, freqs, conns, edge_counts = [], [], [], {}
    for condition in cfg["conditions"]:
        run_cfg = json.loads((run / condition / "config.json").read_text(encoding="utf-8"))
        if run_cfg.get("smoke_test") or run_cfg["n_focal"] != run_cfg["n_reference"]:
            raise ValueError(f"{condition}: not a full run")
        cells, edges, nodes = load_condition(run / condition)
        if nodes != run_cfg["nodes"] or len(nodes) != panel["panel_size"]:
            raise ValueError(f"{condition}: nodes do not match the frozen panel")
        summaries.append(graph_summary(condition, cells, edges, len(target_nodes)))
        freqs.append(edge_frequency(condition, cells, edges, nodes))
        conns.append(target_connectivity(condition, cells, edges, nodes, target_nodes))
        edge_counts[condition] = cells.n_edges[cells.error.eq("")].to_numpy()
    summary = pd.DataFrame(summaries)
    freq = pd.concat(freqs, ignore_index=True)
    conn = add_fold(pd.concat(conns, ignore_index=True), split)
    for frame, name in [(summary, "cscn_panel_graph_summary.csv"), (freq, "cscn_panel_edge_frequency.csv"),
                        (conn, "cscn_panel_target_connectivity.csv")]:
        assert_no_cell_ids(frame)
        frame.to_csv(out / name, index=False, lineterminator="\n")
        print(f"wrote {out / name} ({len(frame)} rows)")
    timing = pd.read_csv(out / "cscn_panel_timing.csv")
    sensitivity = pd.read_csv(out / "cscn_panel_sensitivity.csv")
    text = markdown(summary, conn, freq, panel, cfg, timing, sensitivity)
    (out / "cscn_panel_summary.md").write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {out / 'cscn_panel_summary.md'}")
    if not args.no_figure:
        figure(edge_counts, conn)
        print("wrote reports/figures/cscn_panel_graph_overview.{png,svg}")


if __name__ == "__main__":
    main()
