class microSV:
    def __init__(self, contig, pos, length):
        self.contig = contig
        self.pos = pos
        self.length = length
        self.support = 1
        self.mult = 1

    def __str__(self):
        return "{}, {}, {}".format(self.contig, self.pos, self.length)

class microDeletion(microSV):
    def getType(self):
        return "deletion"

class microDuplication(microSV):
    def __init__(self, contig, pos, length, mult):
        super().__init__(contig, pos, length)
        self.mult = mult

    def getType(self):
        return "duplication"

    def __str__(self):
        return "{}, {}, {}, mult: {}".format(self.contig, self.pos, self.length, self.mult)

class microInsertion(microSV):
    def getType(self):
        return "insertion"


class microInversion(microSV):
    def getType(self):
        return "inversion"