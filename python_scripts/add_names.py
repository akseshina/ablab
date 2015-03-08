__author__ = 'Margarita'

import sys


filename = sys.argv[1]
name = sys.argv[2]
#filename = "ecoli.fasta"
f = open(filename, 'r')
data = f.readlines()
f.close()

f = open(name + ".fasta", "w")

for line in data:
    if line[0] == '>':
        f.write(line.strip() + '_' + name + '\n')
    else:
        f.write(line)
f.close()