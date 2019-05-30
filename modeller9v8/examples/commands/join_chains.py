# Example for: Chain.join()

# This will take a model containing two chains and join them into one.

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']

mdl = model(env)
mdl.read(file='2abx')

# Join the B chain onto the end of the A chain
mdl.chains['A'].join(mdl.chains['B'])

# Renumber all residues in the new chain starting from 1
for num, residue in enumerate(mdl.chains['A'].residues):
    residue.num = '%d' % (num + 1)

mdl.write(file='2abx-join.pdb')
