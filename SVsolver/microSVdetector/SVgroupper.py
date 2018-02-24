import collections

def group(micro_duplications, micro_deletions, micro_insertions, micro_inversions, config):
    return groupSV(micro_duplications, config), groupSV(micro_deletions, config),\
           groupSV(micro_insertions, config), groupSV(micro_inversions, config)


def groupByContigs(sv_list, sv_contig_map, contigs):
    for sv in sv_list:
        sv.support = 1
        if sv.contig in contigs:
            sv_contig_map[sv.contig].append(sv)
    for contig, sv in sv_contig_map.items():
        sv.sort(key=lambda x: x.pos)


def groupSortedSVByPosition(sv_list, maxSvdistanceToMerge, minSvsupport):
    if not sv_list:
        return []
    for i in range(1, len(sv_list)):
        curr = sv_list[i]
        for j in range(i - 1, -1, -1):
            last = sv_list[j]
            if curr.pos - last.pos < maxSvdistanceToMerge:
                if last.support and abs(curr.length - last.length) < maxSvdistanceToMerge:
                    last.pos = (last.pos + curr.pos) / 2
                    last.length = (last.length + curr.length) / 2
                    last.support += 1
                    last.mult = max(last.mult, curr.mult)
                    curr.support = 0
                    break
            else:
                break

    return [sv for sv in sv_list if sv.support >= minSvsupport]


def groupSVmapByPosition(sv_contig_map, maxSvdistanceToMerge, minSvsupport):
    result = []
    for contig, sv in sv_contig_map.items():
        result.extend(groupSortedSVByPosition(sv, maxSvdistanceToMerge, minSvsupport))
    return result


def groupSV(sv_list, config):
    sv_contig_map = collections.defaultdict(list)
    groupByContigs(sv_list, sv_contig_map, config.contigs)
    return groupSVmapByPosition(sv_contig_map, config.max_sv_distance_to_merge, config.min_micro_sv_support)


