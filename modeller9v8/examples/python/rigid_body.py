from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
mdl = model(env, file='1fas')

# Keep residues 1-10 in chain A rigid:
r = rigid_body(mdl.residue_range('1:A', '10:A'))
mdl.restraints.rigid_bodies.append(r)

# Randomize the coordinates of the whole model; the rigid body remains rigid
sel = selection(mdl)
sel.randomize_xyz(deviation=4.0)
mdl.write(file='1fas.ini')
