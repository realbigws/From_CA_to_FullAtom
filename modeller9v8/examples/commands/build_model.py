# Example for: model.build()
# This will build a model for a given sequence in an extended conformation.

from modeller import *
env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read the sequence from a file (does not have to be part of an alignment):
aln = alignment(env, file='toxin.ali', align_codes='1fas')
# Calculate its molecular topology:
mdl = model(env)
mdl.generate_topology(aln['1fas'])
# Calculate its Cartesian coordinates using internal coordinates and
# parameters if necessary:
mdl.build(initialize_xyz=True, build_method='INTERNAL_COORDINATES')

# Add PDB remarks for readability
mdl.remark = """REMARK   4 Extended-chain model of 1fas
REMARK   4 Built from internal coordinates only"""

# Write the coordinates to a PDB file:
mdl.write(file='1fas.ini')
