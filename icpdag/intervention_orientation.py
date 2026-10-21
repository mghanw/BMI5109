"""Seed an observed PDAG with compatible intervention-response directions."""

import networkx as nx


def fixed_edges(graph):
    return {(a,b) for a,b in graph.edges if not graph.has_edge(b,a)}


def colliders(graph):
    arrows=fixed_edges(graph)
    result=set()
    for center in graph:
        parents=sorted(a for a,b in arrows if b==center)
        for i,a in enumerate(parents):
            for b in parents[i+1:]:
                if not graph.has_edge(a,b) and not graph.has_edge(b,a):
                    result.add((a,center,b))
    return result


def valid_extension(pdag,dag):
    pairs=lambda graph:{frozenset((a,b)) for a,b in graph.edges}
    return (set(pdag)==set(dag) and pairs(pdag)==pairs(dag)
            and fixed_edges(pdag)<=set(dag.edges)
            and colliders(pdag)==colliders(dag)
            and nx.is_directed_acyclic_graph(dag))


def seed_interventions(collider_graph,evidence,complete,extension_check=valid_extension):
    """Keep PC collider directions; reject conflicting seed sets symmetrically.

    Evidence contains named, experimentally targeted genes, never a gene-to-NMF
    assignment. Nonadjacent responses remain in the evidence table, not the graph.
    extension_check decides whether a completion is consistent with the seeded
    graph (default: valid_extension, which also requires collider preservation).
    """
    original=collider_graph.copy()
    records={}
    absent=0
    proposals={(str(r["source"]),str(r["destination"])):r for r in evidence}
    candidates=set()
    for (source,destination),row in sorted(proposals.items()):
        if not (original.has_edge(source,destination) or original.has_edge(destination,source)):
            absent+=1
            continue
        record=dict(row)
        records[(source,destination)]=record
        if (destination,source) in proposals:
            record["status"]="bidirectional_response_conflict"
        elif not original.has_edge(source,destination):
            record["status"]="collider_conflict"
        elif not original.has_edge(destination,source):
            record["status"]="agrees_collider"
        else:
            record["status"]="accepted_seed"
            candidates.add((source,destination))

    def apply(seeds):
        graph=original.copy()
        graph.remove_edges_from((b,a) for a,b in seeds)
        return graph

    seeded=apply(candidates)
    directed=nx.DiGraph(); directed.add_nodes_from(seeded); directed.add_edges_from(fixed_edges(seeded))
    cyclic=set()
    for component in nx.strongly_connected_components(directed):
        if len(component)>1:
            cyclic.update((a,b) for a,b in candidates if a in component and b in component)
    for edge in cyclic:
        records[edge]["status"]="cycle_conflict"
    candidates-=cyclic
    seeded=apply(candidates)
    incompatible=set()
    for a,center,b in colliders(seeded)-colliders(original):
        incompatible.update({(a,center),(b,center)} & candidates)
    for edge in incompatible:
        records[edge]["status"]="new_collider_conflict"
    candidates-=incompatible
    seeded=apply(candidates)
    try:
        compatible=extension_check(seeded,complete(seeded))
    except (ValueError,RuntimeError):
        compatible=False
    if not compatible and candidates:
        # Do not select winners by gene order or response p-value.
        for edge in candidates:
            records[edge]["status"]="no_consistent_extension"
        candidates.clear()
        seeded=original.copy()
    return seeded,list(records.values()),absent
