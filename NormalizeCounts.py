import math
from numpy import mean
from numpy import median

countsFile = "counts.txt"
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
