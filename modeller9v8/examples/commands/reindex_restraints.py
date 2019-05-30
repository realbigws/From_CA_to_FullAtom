# Example for: restraints.reindex()

# This will reindex restraints obtained previously for a simpler topology so
# that they will now apply to a more complicated topology.

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.io.atom_files_directory = ['../atom_files']
tpl = env.libs.topology
par = env.libs.parameters

# Generate the model for the simpler topology (CA only in this case):
tpl.read(file='$(LIB)/top_ca.lib')
par.read(file='$(LIB)/par_ca.lib')

code = '1fas'
mdl = complete_pdb(env, code)
mdl.write(file=code+'.ca')

# Generate the restraints for the simpler topology:
sel = selection(mdl)
mdl.restraints.make(sel, restraint_type='stereo', spline_on_site=False)
mdl.restraints.write(file='1fas-ca.rsr')
sel.energy()

# Generate the model for the more complicated topology:
tpl.read(file='$(LIB)/top_heav.lib')
par.read(file='$(LIB)/par.lib')

mdl.read(file=code)
aln = alignment(env)
aln.append_model(mdl, atom_files=code, align_codes=code)
aln.append_model(mdl, atom_files=code+'.ini', align_codes=code+'-ini')
mdl.clear_topology()
mdl.generate_topology(aln[code+'-ini'])
mdl.transfer_xyz(aln)
mdl.write(file='1fas.ini')

mdl2 = model(env, file='1fas.ca')
mdl.restraints.reindex(mdl2)
mdl.restraints.write(file='1fas.rsr')
sel = selection(mdl)
sel.energy()
