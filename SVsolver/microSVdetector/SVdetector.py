from microSVdetector.microSVobjects import *


def detectSV(config, segments):
    micro_inversions = []
    micro_deletions = []
    micro_duplications = []
    micro_insertions = []

    for s in segments:
        chceckForMicroVariations(s, config, micro_duplications, micro_deletions, micro_insertions, micro_inversions)

    return micro_duplications, micro_deletions, micro_insertions, micro_inversions

def chceckForMicroVariations(reads, config, micro_duplications, micro_deletions, micro_insertions, micro_inversions):  # reads have to be sorted by query_alignment_start
    microDup = None
    for i in range(1, len(reads)):
        read = reads[i]
        prev_read = reads[i - 1]
        if not read.reference_name == prev_read.reference_name:
            microDup = None
            continue
        if not read.is_reverse == prev_read.is_reverse:
            if i < len(reads) - 1:
                getInversion(read, prev_read, reads[i + 1], config, micro_inversions)
            microDup = None
            continue
        if not read.is_reverse:
            reference_end = read.reference_end
            reference_start = read.reference_start
            prev_reference_end = prev_read.reference_end
            prev_reference_start = prev_read.reference_start
        else:
            reference_end = read.reference_start
            reference_start = read.reference_end
            prev_reference_end = prev_read.reference_start
            prev_reference_start = prev_read.reference_end
        diff = reference_start - prev_reference_end
        if inNeighbour(read.query_alignment_start, prev_read.query_alignment_end,
                       config.max_alignment_distance) and diff > 0 and diff <= config.max_micro_deletion_length:
            micro_del = microDeletion(read.reference_name, prev_read.reference_end, diff)
            micro_deletions.append(micro_del)
            prev_read.micro_sv.append(micro_del)
            microDup = None
        elif inNeighbour(read.query_alignment_start, prev_read.query_alignment_end,
                         config.max_alignment_distance) or prev_read.query_alignment_end > read.query_alignment_start:
            if diff < 0 and reference_end > prev_reference_start:
                if microDup == None or not inNeighbour(microDup.pos, reference_start, config.max_alignment_distance):
                    microDup = microDuplication(read.reference_name, reference_start,
                                                prev_reference_end - reference_start, 2)
                    micro_duplications.append(microDup)
                else:
                    microDup.mult = microDup.mult + 1
                prev_read.micro_sv.append(microDup)
                if prev_reference_end + config.max_alignment_distance < reference_end:
                    microDup = None
        elif read.query_alignment_start > prev_read.query_alignment_end + config.min_micro_insertion_length:
            if diff >= - config.max_alignment_distance and diff < config.max_alignment_distance:
                micro_insertion = microInsertion(read.reference_name, prev_read.reference_end,
                                                 read.query_alignment_start - prev_read.query_alignment_end)
                micro_insertions.append(micro_insertion)
                prev_read.micro_sv.append(micro_insertion)
            microDup = None
        else:
            microDup = None


def getInversion(read, prev_read, next_read, config, micro_inversions):
    if inNeighbour(read.query_alignment_start, prev_read.query_alignment_end, config.max_alignment_distance) and inNeighbour(
            read.query_alignment_end, next_read.query_alignment_start, config.max_alignment_distance):
        if not read.is_reverse == next_read.is_reverse and prev_read.is_reverse == next_read.is_reverse:
            inversion = microInversion(read.reference_name, min(read.reference_start, read.reference_end),
                                       read.query_alignment_length)
            micro_inversions.append(inversion)
            prev_read.micro_sv.append(inversion)
            read.micro_sv.append(inversion)

def inNeighbour(pos1, pos2, eps):
    return abs(pos1-pos2) < eps