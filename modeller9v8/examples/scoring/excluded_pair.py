# Demonstrate the use of excluded pairs.

# In this example we approximate a disulfide linkage by creating a distance
# restraint between two SG atoms in CYS residues. Since these atoms are in
# different residues, ordinarily Modeller will calculate a van der Waals
# (soft sphere) interaction between them. We use an excluded pair to prevent
# this interaction from being calculated, as otherwise it will conflict
# with the new distance restraint.

# Note that this is an example only; ordinarily a DISU patch would be used
# to create a disulfide linkage. The DISU patch has the advantage that it
# restrains the angles and dihedrals involved with the SG-SG bond, and also
# excludes these atom pairs from van der Waals interaction.

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

code = '1fas'
mdl = complete_pdb(env, code)

atom1 = mdl.atoms['SG:3']
atom2 = mdl.atoms['SG:22']
mdl.restraints.add(forms.gaussian(group=physical.xy_distance,
                                  mean=2.0, stdev=0.1,
                                  feature=features.distance(atom1, atom2)))
mdl.restraints.excluded_pairs.append(excluded_pair(atom1, atom2))

# Retain stereochemistry
atmsel = selection(mdl)
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

# Optimize the model with CG
cg = conjugate_gradients(output='REPORT')
cg.optimize(atmsel, max_iterations=100)

mdl.write(file=code+'.expair.pdb')
