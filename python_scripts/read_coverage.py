__author__ = 'Margarita'

import sys


filename = sys.argv[1]
#filename = "ecoli.fasta"
f = open(filename, 'r')
data = f.readlines()
f.close()

f = open("coverage.txt", "w")

for line in data:
    if line[0] == '>':
        f.write(line.strip()[1:] + '\t')
        line = line.split('_')
        f.write(line[5] + '\n')

f.close()