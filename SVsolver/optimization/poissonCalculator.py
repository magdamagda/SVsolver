from scipy.stats import poisson
import numpy as np

from breakpointGraph.edgeType import EdgeType
from modules.multiprocmap import parmap
import networkx as nx

ALMOST_ZERO = 1e-11


def get_edge_cov_probs(idx_cov):
    idx,cov, edge_type, copy_numbers, mean_coverage = idx_cov
    probs=[get_poisson_probability_for_range(cov, c[0], c[1], mean_coverage) for c in copy_numbers]
    lprobs=[get_log(p) for p in probs]
    return (idx, cov, probs, lprobs, edge_type)

def calculate(G, config):
    initial_x = np.empty(0)
    alphas = np.empty(0)
    edge_log_probabilities_array = np.empty(0)
    edges_props = nx.get_edge_attributes(G, "props")
    edge_idx = {}
    print("Computing Poisson probs. in parallel...")
    edge_coverages = [(idx, edges_props[e].cov, edges_props[e].type, config.copy_numbers, config.mean_coverage) for idx, e in enumerate(G.edges())]
    edge_cov_probs = parmap(get_edge_cov_probs, edge_coverages,
                            percentage_fun=lambda i, p: print(
                                "Finished {} percent".format(round(p, 2))) if i % 10 == 0 else None)
    #print("Concatenating probabilities...")
    for idx_cov_probs_lp,e in zip(edge_cov_probs,list(G.edges())):
        idx, cov, probs, lprobs, edge_type = idx_cov_probs_lp
        #cov = edges_props[e].cov
        #probabilities = np.array([get_poisson_probability_for_range(cov, c[0], c[1]) for c in copy_numbers]) # poisson probabilities
        probabilities = np.array(probs)
        s=np.sum(probabilities)
        print("Prob.sum={}, cov={}".format(s, cov))
        lprobs -= np.log(s)
        initial_x = np.concatenate((initial_x, probabilities/s))
        edge_log_probabilities_array = np.concatenate((edge_log_probabilities_array,
                                                       np.array(lprobs)))
                                                       #np.array([get_log(p) for p in probabilities])))
                                                # This seems ok since get_log looks similar to get_log_poisson_...
        edge_idx[e] = idx
        #print(idx)
        if edge_type == EdgeType.MATE:
            edge_alphas = config.prior_breakpoint
        else:
            edge_alphas = config.prior_segment
        alphas = np.concatenate((alphas, edge_alphas))

    print("Concatenating probabilities... DONE")
    nx.set_edge_attributes(G, 'idx', edge_idx)

    return initial_x, alphas, edge_log_probabilities_array

def get_poisson_probability(coverage, copy_number, mean_coverage):
    if copy_number == 0:
        lambd = ALMOST_ZERO
    else:
        lambd = copy_number * mean_coverage / 2
    return poisson.pmf(coverage, lambd)

def get_poisson_probability_for_range(coverage, copy_number_start, copy_number_stop, mean_coverage):
    result = 0
    for copy_num in range(copy_number_start, copy_number_stop + 1):
        result += get_poisson_probability(coverage, copy_num, mean_coverage)
    if result == 0:
        return ALMOST_ZERO
    return result

#def get_log_poisson_probability_for_range(coverage, copy_number_start, copy_number_stop):
#    result = np.log(get_poisson_probability_for_range(coverage, copy_number_start, copy_number_stop))
#    if np.isneginf(result):
#        return np.log(ALMOST_ZERO)
#    return result

def get_log(p):
    log = np.log(p)
    if np.isneginf(log):
        return np.log(ALMOST_ZERO)
    return log