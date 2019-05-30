# This demonstrates relating PDB residue numbers with residue indices.

from modeller import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['../atom_files']

mdl = model(env, file='2abx')

print "Mapping from residue indices to PDB residue and chain names:"
for r in mdl.residues:
    print "%6d   %3s:%s   %s" % (r.index, r.num, r.chain.name, r.pdb_name)
