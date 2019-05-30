# align2d/align using salign

# parameters to be input by the user
# 1.  gap_penalties_1d
# 2.  gap_penalties_2d
# 3.  input alignment file

from modeller import *
log.verbose()
env = environ()
env.io.atom_files_directory = ['../atom_files']

aln = alignment(env, file='align2d_in.ali', align_codes='all')
aln.salign(rr_file='$(LIB)/as1.sim.mat',  # Substitution matrix used
           output='',
           max_gap_length=20,
           gap_function=True,              # If False then align2d not done
           feature_weights=(1., 0., 0., 0., 0., 0.),
           gap_penalties_1d=(-100, 0),
           gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0.0, 0.0),
           # d.p. score matrix
           #write_weights=True, output_weights_file='salign.mtx'
           similarity_flag=True)   # Ensuring that the dynamic programming
                                   # matrix is not scaled to a difference matrix
aln.write(file='align2d.ali', alignment_format='PIR')
aln.write(file='align2d.pap', alignment_format='PAP')
