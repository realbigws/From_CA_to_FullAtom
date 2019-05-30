# Example for: model.orient()

# This will orient the model along the principal axes of the inertia ellipsoid:

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
mdl = model(env)
mdl.read(file='1fas')
r = mdl.orient()
mdl.write(file='1fas.ini')

print "Translation: ", r.translation
