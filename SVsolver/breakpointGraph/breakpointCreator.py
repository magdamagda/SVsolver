from objects.breakpoint import Breakpoint
from objects.direction import Direction

NOT_MAPPED = "notmapped"

def create(config, segments, CONTIGS_ENDPOINTS):
    breakpoints = {NOT_MAPPED : []}
    for s in segments:
        createBreakPoints(s, config.include_not_mapped, config.min_aligment_length, breakpoints, CONTIGS_ENDPOINTS)
    return breakpoints

def createBreakPoints(reads, include_not_mapped, min_aligment_length, breakpoints, CONTIGS_ENDPOINTS):

    if reads[0].query_alignment_start > min_aligment_length:
        brp1 = Breakpoint(NOT_MAPPED, 0, 0, Direction.NONE)
        brp2 = Breakpoint(reads[0].reference_name, reads[0].reference_start, reads[0].query_alignment_start,
                          getDirection(reads[0]))
        addMateBreakpoints(brp1, brp2, reads[0].query_name, breakpoints, CONTIGS_ENDPOINTS)
    if reads[-1].read_length - reads[-1].query_alignment_end > min_aligment_length:
        brp1 = Breakpoint(NOT_MAPPED, 0, reads[-1].query_alignment_end, Direction.NONE)
        brp2 = Breakpoint(reads[-1].reference_name, reads[-1].reference_end, reads[-1].query_alignment_end,
                          getDirection(reads[-1]))
        addMateBreakpoints(brp1, brp2, reads[0].query_name, breakpoints, CONTIGS_ENDPOINTS)
    for i in range(len(reads) - 1):
        read = reads[i]
        next_read = reads[i + 1]
        if len(read.micro_sv) > 0:
            continue
        # Nice breakpoints without large gaps between aligned segments or overlapped segments
        if inNeighbour(read.query_alignment_end, next_read.query_alignment_start,
                       min_aligment_length) or read.query_alignment_end > next_read.query_alignment_start:
            brp1 = Breakpoint(read.reference_name, read.reference_end, read.query_alignment_end,
                              getReversedDirection(read))
            brp2 = Breakpoint(next_read.reference_name, next_read.reference_start, next_read.query_alignment_start,
                              getDirection(next_read))
            addMateBreakpoints(brp1, brp2, read.query_name, breakpoints, CONTIGS_ENDPOINTS)
        elif include_not_mapped and read.query_alignment_end < next_read.query_alignment_start:
            # print("not mapped length:", next_read.query_alignment_start - read.query_alignment_end)
            # TODO count such cases
            brp1 = Breakpoint(read.reference_name, read.reference_end, read.query_alignment_end,
                              getReversedDirection(read))
            brp2 = Breakpoint(NOT_MAPPED, 0, read.query_alignment_end,
                              Direction.NONE)
            addMateBreakpoints(brp1, brp2, read.query_name, breakpoints, CONTIGS_ENDPOINTS)

            brp3 = Breakpoint(NOT_MAPPED, 0, next_read.query_alignment_end,
                              Direction.NONE)
            brp4 = Breakpoint(next_read.reference_name, next_read.reference_start, next_read.query_alignment_start,
                              getDirection(next_read))

            addMateBreakpoints(brp3, brp4, read.query_name, breakpoints, CONTIGS_ENDPOINTS)

            brp2.setNext(brp3)
            brp3.setNext(brp2)

        else:  # not inNeighbour(...) and read.query_alignment_end >= next_read.query_alignment_start:
            # TODO log such cases
            print(read.reference_name, next_read.reference_name, str(read.query_alignment_end),
                  str(next_read.query_alignment_start), str(read.reference_end), str(next_read.reference_start),
                  str(read.query_alignment_length), str(next_read.query_alignment_length))


def getDirection(read):
    return getDirectionForPos(read.reference_start, read.reference_end)


def getReversedDirection(read):
    return getDirectionForPos(read.reference_end, read.reference_start)


def getDirectionForPos(start, end):
    if start < end:
        return Direction.RIGHT
    else:
        return Direction.LEFT


def addMateBreakpoints(brp1, brp2, read_name, breakpoints, CONTIGS_ENDPOINTS):
    brp1.addRead(read_name)
    brp2.addRead(read_name)
    brp1.addMate(brp2)
    brp2.addMate(brp1)
    addBreakPointToContigs(brp1, breakpoints, CONTIGS_ENDPOINTS)
    addBreakPointToContigs(brp2, breakpoints, CONTIGS_ENDPOINTS)


def addBreakPointToContigs(breakpoint, breakpoints, CONTIGS_ENDPOINTS):
    if breakpoint.contig not in breakpoints:
        breakpoints[breakpoint.contig] = []
        if not breakpoint.contig == NOT_MAPPED:
            endpoints = CONTIGS_ENDPOINTS[breakpoint.contig]
            breakpoints[breakpoint.contig].append(Breakpoint(breakpoint.contig, endpoints[0], None, Direction.RIGHT, True))
            breakpoints[breakpoint.contig].append(Breakpoint(breakpoint.contig, endpoints[1], None, Direction.LEFT, True))
            breakpoints[breakpoint.contig].append(Breakpoint(breakpoint.contig, endpoints[2], None, Direction.RIGHT, True))
            breakpoints[breakpoint.contig].append(Breakpoint(breakpoint.contig, endpoints[3], None, Direction.LEFT, True))
            if breakpoint.pos >= endpoints[1] and breakpoint.pos <= endpoints[2]:
                return
    breakpoints[breakpoint.contig].append(breakpoint)


def inNeighbour(pos1, pos2, eps):
    return abs(pos1 - pos2) < eps