# Example for: alignment.segment_matching()

from modeller import *

log.level(1, 1, 1, 1, 0)
env = environ()

aln = alignment(env, file='ednf2.pap', align_codes=('7rsa', 'edn', 'templ'),
                alignment_format='PAP', remove_gaps=True)

aln.segment_matching(file='segmatch.dat',
                     align_block=1, rr_file='$(LIB)/as1.sim.mat',
                     segment_shifts=(-8, 8, 0, 0),
                     segment_growth_n=(0, 0, 0, 0),
                     segment_growth_c=(0, 0, 0, 0),
                     min_loop_length=(0,2,0),
                     segment_report=1000000, segment_cutoff=0,
                     root_name='segmatch', file_ext='.ali', file_id='default')
