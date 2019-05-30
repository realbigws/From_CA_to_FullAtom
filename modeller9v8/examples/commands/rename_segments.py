# Example for: model.rename_segments()

# This will assign new PDB single-character chain id's to all the chains
# in the input PDB file (in this example there are two chains).

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
mdl = model(env, file='2abx')

# Assign new segment names and write out the new model:
mdl.rename_segments(segment_ids=('X', 'Y'))
mdl.write(file='2abx-renamed.pdb')
