import collections

class Breakpoint:
    breakpoint_id = 0

    def __init__(self, contig, pos, posInRead, direct, terminal=False):
        self.contig = contig
        self.pos = pos
        self.posInRead = posInRead
        self.direction = direct
        self.reads = set()
        self.mates = collections.defaultdict(int)
        self.id = Breakpoint.breakpoint_id
        Breakpoint.breakpoint_id += 1
        self.flag = False
        self.terminal = terminal
        self.used_support = 0
        self.next = None  # only for breakpoints in NOT_MAPPED
        self.next_sup = 0

    def __hash__(self):
        return self.id

    def addRead(self, readName):
        self.reads.add(readName)

    def addMate(self, halfBreakPoint):
        self.mates[halfBreakPoint] += 1

    def removeMate(self, mate):
        if mate in self.mates:
            del self.mates[mate]

    def addReads(self, readName):
        self.reads = self.reads | readName

    def addMates(self, mates):
        for m in mates:
            self.mates[m] += 1

    def setNext(self, n):
        self.next = n
        self.next_sup = 1

    def __str__(self):
        return str(self.id) + " " + self.contig + " " + str(self.pos) + " " + str(self.direction) + " " + str(
            self.posInRead)