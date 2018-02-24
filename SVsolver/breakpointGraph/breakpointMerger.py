from breakpointGraph.breakpointCreator import NOT_MAPPED, inNeighbour


def merge(breakpoints, config):
    joinBreakPoints(breakpoints, config.max_breakpoint_distance_to_merge)
    removeUnsupportedBreakpoints(breakpoints, config.min_breakpoint_support)
    mergeNotMapped(breakpoints)


def mergeBreakPoints(brp1, brp2):
    brp1.addReads(brp2.reads)
    brp1.addMates(brp2.mates)
    for mate in brp2.mates:
        # import pdb; pdb.set_trace()
        mate.removeMate(brp2)
        mate.addMate(brp1)
    brp1.pos = int((brp1.pos + brp2.pos) / 2)
    return brp1


def joinBreakPoints(breakpoints, max_breakpoint_distance_to_merge):
    for contig in breakpoints:
        if contig == NOT_MAPPED:
            continue
        breakpoints = breakpoints[contig]
        breakpoints.sort(key=lambda x: x.pos)
        prev_breakpoint = None
        for i in range(0, len(breakpoints)):
            brp = breakpoints[i]
            j = i - 1
            while j >= 0 and (inNeighbour(breakpoints[j].pos, brp.pos, max_breakpoint_distance_to_merge) or breakpoints[j].pos == -1):
                prev_breakpoint = breakpoints[j]
                if prev_breakpoint.pos != -1 and prev_breakpoint.direction == brp.direction and not brp.terminal:
                    mergeBreakPoints(prev_breakpoint, brp)
                    brp.pos = -1
                    break
                j = j - 1
    for contig in breakpoints:
        breakpoints[contig] = [brp for brp in breakpoints[contig] if brp.pos > -1]


def removeUnsupportedBreakpoints(breakpoints, min_breakpoint_support):
    for contig in breakpoints:
        if contig == NOT_MAPPED:
            continue
        breakpoints = breakpoints[contig]
        for brp in breakpoints:
            if brp.pos > -1 and not brp.terminal:
                for mate in list(brp.mates):
                    common_reads = brp.reads.intersection(mate.reads)
                    if (brp.mates[mate] < min_breakpoint_support) and mate.contig != NOT_MAPPED:
                        mate.removeMate(brp)
                        brp.removeMate(mate)
                        mate.reads = mate.reads.difference(common_reads)
                        brp.reads = brp.reads.difference(common_reads)
                        if mate.next:
                            mate.next.next = None
                            # if len(mate.reads) < MIN_BREAKPOINT_COV:
                            #    mate.pos = -1
                            # if len(brp.reads) < MIN_BREAKPOINT_COV:
                            #    brp.pos = -1
    for contig in breakpoints:
        if contig == NOT_MAPPED:
            continue
        breakpoints = breakpoints[contig]
        for brp in breakpoints:
            if brp.pos > -1:
                if (len(brp.mates) == 0 or len(brp.reads) < min_breakpoint_support) and not brp.terminal:
                    brp.pos = -1
                    for mate in list(brp.mates):
                        mate.removeMate(brp)
                        common_reads = brp.reads.intersection(mate.reads)
                        mate.reads = mate.reads.difference(common_reads)
                        if len(mate.mates) == 0:
                            mate.pos = -1
                        if mate.next:
                            mate.next.next = None
    for contig in breakpoints:
        breakpoints[contig] = [brp for brp in breakpoints[contig] if brp.pos > -1]


def mergeNotMapped(breakpoints):
    for contig in breakpoints:
        if contig == NOT_MAPPED:
            continue
        breakpoints = breakpoints[contig]
        for brp in breakpoints:
            not_mapped = [m for m in brp.mates if m.contig == NOT_MAPPED]
            without_next = None
            without_next_sup = 0
            next_mates = {}
            for m in not_mapped:
                if m.next is None:
                    without_next_sup += 1
                    if without_next is None:
                        without_next = m
                    else:
                        m.pos = -1
                        del brp.mates[m]
                else:
                    next_mate = list(m.next.mates)[0]
                    if not next_mate in next_mates:
                        next_mates[next_mate] = [m, 1]
                    else:
                        next_mates[next_mate][1] += 1
                        m.pos = -1
                        del brp.mates[m]

            if not without_next is None:
                brp.mates[without_next] = without_next_sup

            for key, value in next_mates.items():
                brp.mates[value[0]] = value[1]
                value[0].next_sup = value[1]
                value[0].next.next_sup = value[1]
    breakpoints[NOT_MAPPED] = [brp for brp in breakpoints[NOT_MAPPED] if brp.pos > -1]