# Example for alnsequence objects

from modeller import *

env = environ()

aln = alignment(env, file='../commands/toxin.ali')
print "Alignment contains %d sequences:" % len(aln)
for seq in aln:
    print "   Sequence %s from %s contains %d residues" \
          % (seq.code, seq.source, len(seq))
