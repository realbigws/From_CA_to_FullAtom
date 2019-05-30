from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients

env = environ()

env.io.atom_files_directory = ['../atom_files']
log.verbose()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read in the model
mdl = complete_pdb(env, "1fdn")
rsr = mdl.restraints

# Select all C-alpha atoms
allat = selection(mdl)
allca = allat.only_atom_types('CA')

# Create a pseudo atom that is the center of all C-alphas, and activate it
center = pseudo_atom.gravity_center(allca)
rsr.pseudo_atoms.append(center)

# Constrain every C-alpha to be no more than 10 angstroms from the center
for at in allca:
    r = forms.upper_bound(group=physical.xy_distance,
                          feature=features.distance(at, center),
                          mean=10.0, stdev=0.1)
    rsr.add(r)

# Constrain the gravity center to the x=0 plane
r = forms.gaussian(group=physical.xy_distance,
                   feature=features.x_coordinate(center),
                   mean=0.0, stdev=0.1)
rsr.add(r)

# Keep sensible stereochemistry
rsr.make(allat, restraint_type='stereo', spline_on_site=False)

# Optimize with CG
cg = conjugate_gradients()
cg.optimize(allat, max_iterations=100, output='REPORT')
mdl.write(file='1fas.ini')
