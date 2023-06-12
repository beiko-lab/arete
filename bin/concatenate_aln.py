#!/usr/bin/env python
import os
import sys
from collections import defaultdict

sequences = defaultdict()
alnlengths = defaultdict()
allgenomes = defaultdict()

alignments = defaultdict()

print("Reading alignments...", file=sys.stderr)

for alnName in sys.stdin:
    alnName = alnName.strip()
    filepath = alnName + ".aln"
    if os.stat(filepath).st_size != 0:
        curSeq = ""
        with open(filepath, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    curSeq = line[1:]
                    allgenomes[curSeq] = 1
                else:
                    alignments.setdefault(alnName, {}).setdefault(curSeq, "")
                    alignments[alnName][curSeq] += line

        alnlengths[alnName] = len(alignments[alnName][curSeq])

genomelist = sorted(allgenomes.keys())

concatenated = {}

print("Building concatenated alignment...", file=sys.stderr)

for outgenome in genomelist:
    print(">" + outgenome)
    for alignment in sorted(alignments.keys()):
        if alignment in alnlengths:
            if outgenome in alignments[alignment]:
                print(alignments[alignment][outgenome], end="")
            else:
                print("-" * alnlengths[alignment], end="")
    print("\n")
