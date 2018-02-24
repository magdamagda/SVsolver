from enum import Enum


class EdgeType(Enum):
    MATE = 1
    NEXT_READ = 2 # next by pos in chrom
    SAME_BREAKPOINT = 3