# Example for 'chain' objects

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, "1b3q")

# Print existing chain IDs and lengths:
print "Chain IDs and lengths: ", [(c.name, len(c.residues)) for c in mdl.chains]

# Set new chain IDs:
mdl.chains['A'].name = 'X'
mdl.chains['B'].name = 'Y'

# Write out chain sequences:
for c in mdl.chains:
    c.write(file='1b3q%s.chn' % c.name, atom_file='1b3q',
            align_code='1b3q%s' % c.name)
