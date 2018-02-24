import csv
import numpy as np
import os.path

OPTYM_LOG_FILE_NAME = 'optymalisationLog.csv'
PICKLED_RESULT_FILE_NAME = "result.npy"
RESULT_FILE_NAME = 'result.txt'

def save(outputDirectory, x, results):
    saveOptimLog(outputDirectory, results)
    saveOptimResult(outputDirectory, x)

def saveOptimLog(outputDirectory, results):
    with open(os.path.join(outputDirectory, OPTYM_LOG_FILE_NAME), 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=';')
        csvwriter.writerow(["goal function", "error", "goal function change", "error change", "result difference"])
        for r in results:
            csvwriter.writerow([str(x) for x in r])

def saveOptimResult(outputDirectory, x):
    np.save(os.path.join(outputDirectory, PICKLED_RESULT_FILE_NAME), x)
    with open(os.path.join(outputDirectory, RESULT_FILE_NAME), 'w') as f:
        for i in x:
            f.write(str(i) + "\n")