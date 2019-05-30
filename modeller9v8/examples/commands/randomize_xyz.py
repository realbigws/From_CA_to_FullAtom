# Example for: selection.randomize_xyz()

# This will randomize the X,Y,Z of the model:

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']

mdl = model(env, file='1fas')

# Act on all atoms in the model
sel = selection(mdl)

# Change all existing X,Y,Z for +- 4 angstroms:
sel.randomize_xyz(deviation=4.0)
mdl.write(file='1fas.ini1')

# Assign X,Y,Z in the range from -100 to 100 angstroms:
sel.randomize_xyz(deviation=-100.0)
mdl.write(file='1fas.ini2')
