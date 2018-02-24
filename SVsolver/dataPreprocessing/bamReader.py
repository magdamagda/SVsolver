import pysam

from objects.segment import Segment


def readAndFilterBam(config):
    print("Read bam", config.bam_file)

    bam_file_all = pysam.AlignmentFile(filepath_or_object=config.bam_file, mode='rb')
    index = pysam.IndexedReads(bam_file_all)
    index.build()

    reads_names = set()
    for ch in config.contigs:
        for read in bam_file_all.fetch(ch, 0, ):
            reads_names.add(read.query_name)

    result = []
    for read_name in reads_names:
        reads = [r for r in index.find(read_name) if validRead(r, config) and not r.is_secondary]
        readsInChromosoms = [r for r in reads if r.reference_name in config.contigs]
        if len(reads) < 2 or not readsInChromosoms:
            continue
        segments = [Segment(read) for read in reads]
        segments.sort(key=lambda x: x.query_alignment_start)
        read_aligment = calculateAlignment(segments)
        # print("alignment", read_aligment * 100, "%")
        if read_aligment < config.min_read_aligment:
            continue
        result.append(segments)

    return result, bam_file_all

def validRead(read, config):
    return read.reference_name in config.contigs_names and read.query_alignment_length > config.min_aligment_length

def calculateAlignment(reads): # sorted
    s = 0
    last_pos = None
    for r in reads:
        if last_pos is None or last_pos <= r.query_alignment_start:
            s += r.query_alignment_length
            last_pos = r.query_alignment_end
        else:
            s += max(0, r.query_alignment_end - last_pos)
            last_pos = max(last_pos, r.query_alignment_end)
    return s/reads[0].read_length