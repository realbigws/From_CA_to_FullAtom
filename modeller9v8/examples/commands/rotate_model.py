# Example for: selection.transform(), selection.translate(),
# selection.rotate_origin()

# This will orient a model as specified:

from modeller import *

# Read the structure:
env = environ()
env.io.atom_files_directory = ['../atom_files']
mdl = model(env, file='1fas')
# Select all atoms
s = selection(mdl)

# Translate 1 angstrom along the x axis:
s.translate([1, 0, 0])

# Transform with a rotation matrix (no change in this example):
s.transform([[1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]])

# Rotate 90 degrees about the axis, through the origin:
s.rotate_origin([1, 1, 1], 90)

mdl.write(file='1fas.ini')
