import collections
import numpy as np

from breakpointGraph.breakpointCreator import NOT_MAPPED
from breakpointGraph.edgeProperties import EdgeProperties
from breakpointGraph.edgeType import EdgeType
from objects.breakpoint import Breakpoint
from objects.direction import Direction

import networkx as nx


def create(breakpoints, bam_file_coverage, coverage):
    G = nx.Graph()
    makeEdges(breakpoints, G, bam_file_coverage, coverage)
    subgraphs = getSubgraphs(G, breakpoints)
    removeTwoNodesSubGraphs(G, subgraphs)
    print("Nodes: ", G.nodes())
    print("Edges: ", G.edges())
    return G


def addSameBreakpointWithOppositDirection(brp):
    return Breakpoint(brp.contig, brp.pos, brp.posInRead,
                      Direction.RIGHT if brp.direction == Direction.LEFT else brp.direction == Direction.LEFT)


def getNextLeftBrp(contig, i, breakpoints):
    idx = i
    idx += 1
    while idx < len(breakpoints[contig]) and (breakpoints[contig][idx].direction != Direction.LEFT or breakpoints[contig][i].pos == breakpoints[contig][idx].pos):
        idx += 1
    return breakpoints[contig][idx]


def makeEdges(breakpoints, G, bam_file_coverage, coverage):
    oppositBreakPoints = collections.defaultdict(list)
    for contig in breakpoints:
        # if contig == NOT_MAPPED:
        #    continue
        contig_breakpoints = breakpoints[contig]
        for idx, brp in enumerate(contig_breakpoints):
            # print(contig, brp.pos)
            for m in brp.mates:
                addEdge(G, brp, m, EdgeType.MATE, 0, None, bam_file_coverage, coverage)
            if not brp.terminal and brp.contig != NOT_MAPPED:
                same_brp = addSameBreakpointWithOppositDirection(brp)
                addEdge(G, brp, same_brp, EdgeType.SAME_BREAKPOINT, 0, None, bam_file_coverage, coverage)
                oppositBreakPoints[brp.contig].append(same_brp)
            if brp.contig == NOT_MAPPED and brp.next != None:
                addEdge(G, brp, brp.next, EdgeType.NEXT_READ, abs(brp.posInRead - brp.next.posInRead), brp.next_sup, bam_file_coverage, coverage)

    for contig in breakpoints:
        if contig == NOT_MAPPED:
            continue
        contig_breakpoints = breakpoints[contig]
        contig_breakpoints += oppositBreakPoints[contig]
        contig_breakpoints.sort(key=lambda x: x.pos)
        for idx, brp in enumerate(contig_breakpoints):
            if brp.direction == Direction.RIGHT:
                #print(brp)
                next_read = getNextLeftBrp(contig, idx, breakpoints)
                addEdge(G, brp, next_read, EdgeType.NEXT_READ, abs(next_read.pos - brp.pos), None, bam_file_coverage, coverage)

def addEdge(G, brp1, brp2, edge_type, l, c, bam_file_coverage, coverage):
    if not c is None:
        cov = c
    elif edge_type == EdgeType.MATE:
        cov = brp1.mates[brp2]
    elif edge_type == EdgeType.SAME_BREAKPOINT:
        cov = countCoverageWithoutBreakpoints(bam_file_coverage, brp1.contig, brp1.pos)
    else:
        cov = countCoverage(brp1.contig, brp1.pos, brp2.pos, coverage)
    G.add_edge(brp1, brp2, props = EdgeProperties(cov, edge_type, l))

def check_read(read, pos):
    return read.reference_start < pos -200 and read.reference_end > pos + 200

def countCoverageWithoutBreakpoints(bam_file_coverage, contig, pos):
    return bam_file_coverage.count(contig, min(0, pos - 200), pos + 200, read_callback=lambda read : check_read(read, pos))

def countCoverage(contig, pos1, pos2, coverage):
    #print("countCoverage", contig, pos1, pos2)
    if contig not in coverage:
        return 0
    start = min(pos1, pos2)
    end = max(pos1, pos2)
    l = end - start
    return int(np.rint(np.sum(coverage[contig][start:end])/l))


def getSubgraphs(G, breakpoints):
    result = []
    for contig in breakpoints:
        for brp in breakpoints[contig]:
            brp.flag = False
    # import pdb; pdb.set_trace()
    for brp in G.nodes():
        if not brp.flag and brp.pos > -1:
            brp.flag = True
            q = [brp]
            subNodes = [brp]
            while not len(q) == 0:
                node = q.pop()
                for n in nx.all_neighbors(G, node):
                    if not n.flag:
                        n.flag = True
                        q.append(n)
                        subNodes.append(n)
            result.append(nx.subgraph(G, subNodes))
            #print(str(len(subNodes)))

    return result

def removeTwoNodesSubGraphs(G, subs):
    nodes_to_remove = []
    for sub in subs:
        if len(sub.nodes()) == 2:
            nodes_to_remove += sub.nodes()
    G.remove_nodes_from(nodes_to_remove)