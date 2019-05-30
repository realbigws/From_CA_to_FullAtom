# This will add models to the alignment.

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']

aln = alignment(env, file='toxin.ali', align_codes='2ctx', remove_gaps=False)
code = "2abx"
mdl = model(env, file=code, model_segment=('1:A', '74:A'))
for n in range(1, 4):
    aln.append_model(mdl, align_codes=code+"_9999%04d" % n,
                          atom_files=code+".B9999%04d" % n)
aln.write(file='toxin-expand.ali')
