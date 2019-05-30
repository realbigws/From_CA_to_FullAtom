# Example for: alignment.describe()

# Describe the sequences and structures in the alignment.

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
aln = alignment(env, file='toxin.ali', align_codes=('2ctx', '2abx'))
aln.describe()
