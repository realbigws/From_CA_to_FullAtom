# Illustrates the SALIGN iterative multiple structure alignment
from modeller import *
import modeller.salign

log.none()
env = environ()
env.io.atom_files_directory = ['.', '../atom_files']

aln = alignment(env)
for (code, chain) in (('1is4', 'A'), ('1uld', 'D'), ('1ulf', 'B'),
                      ('1ulg', 'B'), ('1is5', 'A')):
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

modeller.salign.iterative_structural_align(aln)

aln.write(file='1is3A-it.pap', alignment_format='PAP')
aln.write(file='1is3A-it.ali', alignment_format='PIR')
