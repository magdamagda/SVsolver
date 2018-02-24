from breakpointGraph.edgeType import EdgeType
import networkx as nx
import numpy as np

def prepareConstraintMatrixes(G, breakpoints, config):
    C = []
    A = []
    edges = len(list(G.edges()))
    xi_len = len(config.copy_numbers)
    edge_idx = nx.get_edge_attributes(G, "idx")
    for edge in G.edges():
        idx = edge_idx[edge]
        C.append([0] * (idx * xi_len) + [1] * xi_len + [0] * ((edges - idx - 1) * xi_len))
        A.append(1)

    brpToFlowConstraint = []
    for contig in breakpoints:
        brpToFlowConstraint += [brp for brp in breakpoints[contig] if len(brp.mates) > 0 and brp.pos > -1]

    for brp in brpToFlowConstraint:
        inEdges, outEdges = getInOutEdgesIdx(brp, edge_idx)
        row = [0] * (edges * xi_len)
        for idx in inEdges:
            for j, copy in enumerate(config.copy_numbers):
                row[idx * xi_len + j] = (copy[0] + copy[1]) / 2
        for idx in outEdges:
            for j, copy in enumerate(config.copy_numbers):
                row[idx * xi_len + j] = - (copy[0] + copy[1]) / 2
        C.append(row)
        A.append(0)

    C = np.array(C)
    A = np.array(A)

    return C, A

def getInOutEdgesIdx(G, brp, edge_idx):
    edges_props = nx.get_edge_attributes(G, "props")
    neighbours = [(n, getEdgeProp(edges_props, brp, n)) for n in nx.all_neighbors(G, brp)]
    mates = []
    nexts = []
    same = []
    for n in neighbours:
        if n[1].type == EdgeType.SAME_BREAKPOINT:
            same.append(n[0])
        elif n[1].type == EdgeType.MATE:
            mates.append(n[0])
        elif n[1].type == EdgeType.NEXT_READ:
            nexts.append(n[0])
    if nexts:
        in_edges = [getEdgeProp(edge_idx, brp, nexts[0])]
    else:
        return [], []
        # in_edges = []
    mates_idx = [getEdgeProp(edge_idx, brp, m) for m in mates]
    if same:
        mates_idx.append(getEdgeProp(edge_idx, brp, same[0]))

    return in_edges, mates_idx

def getEdgeProp(edges_props, n1, n2):
    try:
        return edges_props[n1, n2]
    except Exception:
        pass
    try:
        return edges_props[n2, n1]
    except Exception:
        return None