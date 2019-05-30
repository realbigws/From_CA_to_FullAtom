from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = model(env, file='1fdn')

# New empty selection
s = selection()

# Add all atoms from residues 4 through 10 (chain A) inclusive (PDB numbering)
s.add(mdl.residue_range('4:A', '10:A'))

# Selection of all atoms currently within 5A of atom CA in residue 1 in chain A
# (this destroys the previous selection):
s = mdl.atoms['CA:1:A'].select_sphere(5)

# Is the CB:1:A atom in the selection?
print mdl.atoms['CB:1:A'] in s

# Alternative ways of selecting the same atom:
print mdl.chains['A'].residues['1'].atoms['CB'] in s
print mdl.residues['1:A'].atoms['CB'] in s

# All atoms currently within 5A of atom CA:1:A, OR currently within 3A of the
# point (1,10,1):
s = mdl.atoms['CA:1:A'].select_sphere(5) | mdl.point(1,10,1).select_sphere(3)

# All atoms currently within 5A of atom CA:1:A, AND also currently within 3A
# of the point (1,10,1):
s = mdl.atoms['CA:1:A'].select_sphere(5) & mdl.point(1,10,1).select_sphere(3)

# All atoms currently within 5A of atom CA:1:A, OR currently within 3A of the
# point (1,10,1), but not BOTH:
s = mdl.atoms['CA:1:A'].select_sphere(5) ^ mdl.point(1,10,1).select_sphere(3)

# Create a selection containing the CA atom from residue 1, chain A,
# and all of residue 2 (PDB numbering)
s = selection(mdl.atoms['CA:1:A'], mdl.residues['2:A'])

# All residues EXCEPT 5-10 in chain A (i.e. all atom selection minus the
# selection of residues 5-10, otherwise known as an inverted selection):
s = selection(mdl) - selection(mdl.residue_range('5:A', '10:A'))

# Selection of residues 1, 4, 8 and 10-15 (PDB numbering) from chain A:
s = selection(mdl.residues['1:A'], mdl.residues['4:A'], mdl.residues['8:A'],
              mdl.residue_range('10:A', '15:A'))

# Print the center of mass (note: not mass weighted)
print s.mass_center

# Rotate by 90 degrees about the z axis through the origin (0,0,0)
# (right handed rotation)
s.rotate_origin([0,0,1], 90)

# The same thing, except that the axis passes through the center of mass:
s.rotate_mass_center([0,0,1], 90)

# Translate by 5 angstroms along the x axis
s.translate([5.0, 0, 0])

# Equivalent (but less efficient, as it involves calculating the COM)
s.x += 5.0
