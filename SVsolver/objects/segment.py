class Segment():
    def __init__(self, read):
        self.read_length = read.infer_query_length()
        if read.is_reverse:
            self.query_alignment_start = read.infer_query_length() - read.query_alignment_end
            self.query_alignment_end = self.query_alignment_start + read.query_alignment_length
            self.reference_start = read.reference_end
            self.reference_end = read.reference_start
        else:
            self.query_alignment_start = read.query_alignment_start
            self.query_alignment_end = read.query_alignment_end
            self.reference_start = read.reference_start
            self.reference_end = read.reference_end
        self.query_alignment_length = read.query_alignment_length
        self.reference_name = read.reference_name
        self.query_name = read.query_name
        self.micro_sv = []
        self.is_reverse = read.is_reverse