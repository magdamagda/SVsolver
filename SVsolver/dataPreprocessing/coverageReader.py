import multiprocessing as mp


def get_coverage_for(params):
    contig, coverage_dir, coverage_prefix, contig_len = params
    print(contig)
    with open(coverage_dir + coverage_prefix + contig + ".txt", "r") as f:
        result = [0] * (contig_len + 1)
        line = f.readline()
        while line != '':
            splitted = line.split("\t")
            pos = int(splitted[1])
            if pos <= contig_len:
                result[pos] = int(splitted[2])
            line = f.readline()
        res_tuple = (contig, result)
    return res_tuple


def read_coverage(contigs_names, coverage_dir, coverage_prefix, contigs_endpoints):
    print("Read coverage from", coverage_dir)
    pool = mp.Pool(mp.cpu_count())
    res_iter = pool.imap_unordered(get_coverage_for, [(n, coverage_dir, coverage_prefix, contigs_endpoints[n][3]) for n in contigs_names])
    coverage = {}
    coverage.update(res_iter)
    return coverage
