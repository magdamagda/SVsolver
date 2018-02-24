import os.path
import networkx as nx

BREAKPOINT_INFO_FILE_NAME = 'breakpointgraphInfo.txt'
PICKLED_GRAPH_FILE_NAME = 'breakpointgraph.gpickle'

def save(G, breakpoints, output_directory):
    writeGraphInfo(G, breakpoints, output_directory)
    nx.write_gpickle(os.path.join(output_directory, PICKLED_GRAPH_FILE_NAME))


def writeGraphInfo(G, breakpoints, output_directory):
    with open(os.path.join(output_directory, BREAKPOINT_INFO_FILE_NAME), 'w') as f:
        f.write("nodes = {}, edges = {}\n".format(len(G.nodes()), len(G.edges())))
        for contig in breakpoints:
            f.write("{}: {}\n".format(contig, len(breakpoints[contig])))

