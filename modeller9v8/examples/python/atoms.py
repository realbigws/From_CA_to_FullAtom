# Example for 'atom' objects

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, "1fas")

# 'mdl.atoms' is a list of all atoms in the model
print "Name of C-alpha atom in residue 4: ", mdl.atoms['CA:4'].name
a = mdl.atoms[0]
print "Coordinates of first atom: %.3f, %.3f, %.3f" % (a.x, a.y, a.z)

# Each 'residue' object lists its own atoms, as does each chain
a = mdl.residues['10'].atoms[0]
print "Biso for first atom in residue 10: %.3f" % a.biso

a = mdl.chains[0].residues[-1].atoms[-1]
print "Biso for last atom in last residue in first chain: %.3f" % a.biso
