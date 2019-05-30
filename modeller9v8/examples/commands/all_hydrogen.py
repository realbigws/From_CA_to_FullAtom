# This will read a specified atom file, generate all hydrogen atoms,
# add atomic radii and charges, and write the model to a PDB file in
# the GRASP format. This can be used with GRASP to display electrostatic
# properties without assigning charges and radii in GRASP.

from modeller import *
from modeller.scripts import complete_pdb

log.verbose()
env = environ()
env.io.atom_files_directory = ['../atom_files']

env.libs.topology.read(file='$(LIB)/top_allh.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

def patch_disulfides(mdl):
    """Patch topology to remove sulfhydril hydrogens"""
    for ids in [ ('17', '39'),
                 ( '3', '22'),
                 ('53', '59'),
                 ('41', '52') ]:
        mdl.patch(residue_type='DISU', residues=[mdl.residues[r] for r in ids])

mdl = complete_pdb(env, "1fas", patch_disulfides)

mdl.write(file='1fas.ini1', model_format='GRASP')
mdl.write(file='1fas.ini2', model_format='PDB')
