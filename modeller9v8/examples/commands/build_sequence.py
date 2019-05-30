# This demonstrates the use of alignment.append_sequence() and
# model.build_sequence() to build residue sequences from one-letter codes

from modeller import *
env = environ()

# Read parameters (needed to build models from internal coordinates)
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')

# Create a new empty alignment and model:
aln = alignment(env)
mdl = model(env)

# Build a model from one-letter codes, and write to a PDB file:
mdl.build_sequence("AFVVTDNCIK/CKYTDCVEVC")
mdl.write(file='sequence.pdb')

# Build an alignment from one-letter codes
aln.append_sequence("AF---VVTDN---CIKCK------")
aln.append_sequence("-------AFVVTDN--CI--K-CK")
# Set alignment information, and write to file:
aln[0].code = 'seq1'
aln[1].code = 'seq2'
aln.write(file='sequence.ali')
