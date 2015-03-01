__author__ = 'Margarita'
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

input_file = sys.argv[1]
id_file = sys.argv[2]
output_file = sys.argv[3]

wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
print "Found %i unique identifiers in %s" % (len(wanted), id_file)

index = SeqIO.index(input_file, "fasta")
records = (index[r] for r in wanted)
count = SeqIO.write(records, output_file, "fasta")
assert count == len(wanted)

print "Saved %i records from %s to %s" % (count, input_file, output_file)