# This script illustrates the use of the swap_atoms_in_res
# argument to the selection.superpose() command:

# Need to make sure that the topologies of the two molecules
# superposed are exactly the same:

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

atfil = '../atom_files/pdb1fdn.ent'
mdl = complete_pdb(env, atfil)
aln = alignment(env)
aln.append_model(mdl, align_codes='orig')

mdl2 = model(env, file='1fdn.swap.atm')
aln.append_model(mdl2, align_codes='swap')
atmsel = selection(mdl)
atmsel.superpose(mdl2, aln, swap_atoms_in_res='')
atmsel.superpose(mdl2, aln, swap_atoms_in_res='DEFHLNQRVY', fit=False)
atmsel.superpose(mdl2, aln, swap_atoms_in_res='', fit=True)
