# Example for: group_restraints()
from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Allow calculation of statistical (dynamic_modeller) potential
env.edat.dynamic_modeller = True

mdl = complete_pdb(env, "1fas")

# Read Fiser/Melo loop modeling potential
gprsr = group_restraints(env, classes='$(LIB)/atmcls-melo.lib',
                         parameters='$(LIB)/melo1-dist.lib')
# Read DOPE loop modeling potential
#gprsr = group_restraints(env, classes='$(LIB)/atmcls-mf.lib',
#                         parameters='$(LIB)/dist-mf.lib')
# Read DOPE-HR loop modeling potential
#gprsr = group_restraints(env, classes='$(LIB)/atmcls-mf.lib',
#                         parameters='$(LIB)/dist-mfhr.lib')


# Use this potential for the 1fas model
mdl.group_restraints = gprsr

# Evaluate the loop score of PDB residues 1 through 10
atmsel = selection(mdl.residue_range('1', '10'))
atmsel.energy()
