import csv
import os.path

MICRO_SV_FILE_NAME = "microSV.csv"

def save(fileDirectory, micro_duplications, micro_deletions, micro_insertions, micro_inversions):
    with open(os.path.join(fileDirectory, MICRO_SV_FILE_NAME), 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=';')
        csvwriter.writerow(["type", "support", "contig", "position", "length", "multiplicity"])
        writeMicroSVToCSV(micro_duplications, csvwriter)
        writeMicroSVToCSV(micro_deletions, csvwriter)
        writeMicroSVToCSV(micro_insertions, csvwriter)
        writeMicroSVToCSV(micro_inversions, csvwriter)

def writeMicroSVToCSV(microSv, csvwriter):
    for sv in microSv:
        csvwriter.writerow([sv.getType(), sv.support, sv.contig, sv.pos, sv.length, sv.mult])