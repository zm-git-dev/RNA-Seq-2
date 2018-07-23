#!/bin/env python2
import sys
import math
from numpy import mean
from numpy import median

countsFile = sys.argv[1]  

if countsFile == "help" or countsFile == "-h" or countsFile == "--help" or countsFile == "-help" or countsFile.split(".")[-1] != "txt" :
    print ("\nDescription: NormalizeCounts.py normalizes \"counts.txt\" files created with the use of FeatureCounts from the Subread \
package\nnormalization is based on DEseq2 as explaned in the video https://www.youtube.com/watch?v=UFB993xufUU\n\
tijs bliek, 7/20/2018, Amsterdam\
\n usage: ./NormalizeCounts.py [FILENAME] \n\n[FILENAME] \tShould be FeatureCounts counts table (tab delimited .txt),\
\n\t\tcontaining 6 columns of feature information and next multiple columns containing raw counts.")
    sys.exit()

counts = open(countsFile, "r")
samples = []
numbers = {}

for line in counts:
    line = line.split("\t")
    if len(line) > 2:
        if line[0] == "Geneid":
            for i in range(6, len(line)):
                samples.append(line[i])
                numbers[line[i]] = []
        else:
            values = map(int, line[6:])
            if 0 not in values:
                coef = mean(list(map(lambda x: math.log(x), values)))
                norm = list(map(lambda x: math.log(x) - coef, values))
                for i in range(len(norm)):
                    numbers[samples[i]].append(norm[i])
counts.close()

coefs = []
for i in samples:
    coefs.append(math.pow(math.e ,median(numbers[i])))

counts = open(countsFile, "r")
normCounts = "_norm.".join(countsFile.split("."))
uitFile = open(normCounts, "w")


for lin in counts:
    line = lin.split("\t")
    if len(line) < 2 or line[0] == "Geneid":
        uitFile.write(lin)
    else:
        regel = line[:6]
        for x,y in zip(line[6:], coefs):
            regel.append(str(float(x)/y))
        uitFile.write("\t".join(regel) + "\n")
uitFile.close()
counts.close()
