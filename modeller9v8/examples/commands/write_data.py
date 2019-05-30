# Example for: model.write_data()

# This will calculate solvent accessibility, dihedral angles,
# residue-residue neighbors, secondary structure, and local mainchain
# curvature for a structure in the PDB file.

from modeller import *

log.verbose()

# Get topology library for radii and the model without waters and HETATMs:
env = environ()
env.io.atom_files_directory = ['../atom_files']
env.io.hetatm = False
env.io.water = False

env.libs.topology.read(file='$(LIB)/top_heav.lib')
mdl = model(env, file='1fas')

# Calculate residue solvent accessibilities, dihedral angles,
# residue neighbors, and secondary structure, and write to files:
myedat = energy_data()
myedat.radii_factor = 1.0 # The default is 0.82 (for soft-sphere restraints)
mdl.write_data(file='1fas', edat=myedat, output='PSA DIH NGH SSM')

# Calculate local mainchain curvature
mdl.write_data(output='CRV')

# Use the calculated curvature data (in Biso)
print "The following residues have local mainchain curvature"
print "greater than 90 degrees:"
print [r for r in mdl.residues if r.atoms[0].biso * 10 > 90.0]
