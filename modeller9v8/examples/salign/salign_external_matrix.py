# Reads an external matrix
from modeller import *

log.verbose()
env = environ()

aln = alignment(env, file='1dubA-1nzyA.ali', align_codes='all')

aln.salign(alignment_type='pairwise', output='',
           rr_file='$(LIB)/blosum62.sim.mat',
           #rr_file='$(LIB)/as1.sim.mat',
           #max_gap_length=20,
           gap_function=False,
           input_weights_file='external.mtx',   # External weight matrix
           #weights_type='DISTANCE',            # type of ext. wgt. mtx
           # ensure appropriate gap penalites for the ext. matrix
           #feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_1d=(30, 26),
           #write_weights=True, output_weights_file='score.mtx',
           feature_weights=(1., 0., 0., 0., 0., 1.),
           gap_penalties_1d=(-500, -300))
aln.write(file='output.ali', alignment_format='PIR')
aln.write(file='output.pap', alignment_format='PAP')
