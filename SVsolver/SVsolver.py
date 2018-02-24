import sys
import dataPreprocessing.configParser as cp
import dataPreprocessing.cytoBandReader as cbr
import dataPreprocessing.coverageReader as covr
import dataPreprocessing.bamReader as br
import microSVdetector.SVdetector as svd
import microSVdetector.SVgroupper as svg
import resultSaver.microSVresult as microSVresult
import breakpointGraph.breakpointCreator as brpCreator
import breakpointGraph.breakpointMerger as brpMerger
import breakpointGraph.graphCreator as graphCreator
import resultSaver.graphSaver as graphSaver
import optimization.poissonCalculator as poisson
import optimization.constraintsCalculator as constraints
import optimization.augmentedLagrangian as al
import resultSaver.optimResult as optimResult


def readData(configFileName):
    config = cp.parseConfig(configFileName)
    contigs_endpoints = cbr.getContigsEndPoints(config.contigs_names, config.cyto_band)
    coverage = covr.read_coverage(config.contigs_names, config.coverage_directory, config.coverage_file_prefix, contigs_endpoints)
    cbr.correctEndpointsByTelomersLength(contigs_endpoints, coverage)
    segments, bam_file = br.readAndFilterBam(config)
    return config, contigs_endpoints, coverage, segments, bam_file

def microSVdetect(config, segments):
    micro_duplications, micro_deletions, micro_insertions, micro_inversions = svd.detectSV(config, segments)
    micro_duplications, micro_deletions, micro_insertions, micro_inversions = svg.group(micro_duplications, micro_deletions, micro_insertions, micro_inversions, config)
    microSVresult.save(config.output_path, micro_duplications, micro_deletions, micro_insertions, micro_inversions)

def createBrpGraph(config, segments, contigs_endpoints, bam_file, coverage):
    breakpoints = brpCreator.create(config, segments, contigs_endpoints)
    brpMerger.merge(breakpoints, config)
    G = graphCreator.create(breakpoints, bam_file, coverage)
    graphSaver.save(G, breakpoints, config.output_path)
    return G, breakpoints

def optimize(config, G, breakpoints):
    initial_x, alphas, edge_log_probabilities_array = poisson.calculate(G, config)
    C, A = constraints.prepareConstraintMatrixes(G, breakpoints, config)
    x, results = al.optimize(config, C, A, initial_x, alphas, edge_log_probabilities_array)
    optimResult.save(config.output_path, x, results)

def main():
    print("Reading data ...")
    config, contigs_endpoints, coverage, segments, bam_file = readData(sys.argv[1])
    print("Micro SV detecting ...")
    microSVdetect(config, segments)
    print("Breakpoint graph creating ...")
    G, breakpoints = createBrpGraph(config, segments, contigs_endpoints, bam_file, coverage)
    print("Optimization ...")
    optimize(config, G, breakpoints)

if __name__ == "__main__":
    #try:
    main()
    #except Exception as e:
    #    print(e)