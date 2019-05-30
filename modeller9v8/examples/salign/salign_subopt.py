from modeller import *

log.verbose()
env = environ()
aln = alignment(env, file='fm07254_test.ali', alignment_format='PIR')
aln.salign(feature_weights=(1., 0, 0, 0, 0, 0), gap_penalties_1d=(-450, -50),
           n_subopt = 5, subopt_offset = 15)

# Convert suboptimal alignment output file into actual alignments
f = open('suboptimal_alignments.out')
for (n, aln) in enumerate(aln.get_suboptimals(f)):
    aln.write(file='fm07254_out%d.ali' % n)
