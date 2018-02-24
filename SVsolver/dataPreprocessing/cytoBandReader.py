import collections

def getContigsEndPoints(contigsNames, cytoBandFile):
    print("Read cytoBand", cytoBandFile)
    CONTIGS_ENDPOINTS = collections.defaultdict(list)
    with open(cytoBandFile, 'r') as f:
        prev_contig = None
        prev_end = None
        for line in f:
            line = line.split("\t")
            contig = line[0]
            if contig in contigsNames:
                if contig != prev_contig:
                    if prev_contig is not None:
                        CONTIGS_ENDPOINTS[prev_contig].append(prev_end)
                    CONTIGS_ENDPOINTS[contig].append(int(line[1]))

                if line[4] == 'acen\n':
                    if line[3][0] == 'p':
                        CONTIGS_ENDPOINTS[contig].append(int(line[1]))
                    elif line[3][0] == 'q':
                        CONTIGS_ENDPOINTS[contig].append(int(line[2]))

                prev_contig = contig
                prev_end = int(line[2])
        CONTIGS_ENDPOINTS[prev_contig].append(prev_end)

    return CONTIGS_ENDPOINTS

def correctEndpointsByTelomersLength(endpoints, coverage):
    for contig, coverages in coverage.items():
        for idx, cov in enumerate(reversed(coverages)):
            if cov!= 0:
                endpoints[contig][0] += idx
                endpoints[contig][3] -= idx
                break