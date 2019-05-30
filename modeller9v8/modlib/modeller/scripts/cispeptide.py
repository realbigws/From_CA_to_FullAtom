from modeller import forms, features, physical

__docformat__ = "epytext en"

def cispeptide(rsr, atom_ids1, atom_ids2):
    """To create cis-peptide bond stereochemical restraints.

       Note: you must have done the STEREOCHEMICAL restraints already.

       Example (cis-peptide between residues 4 and 5):

       >>> a = mdl.atoms
       >>> cispeptide(mdl.rsr,
       ...            atom_ids1=(a['O:4'], a['C:4'], a['N:5'], a['CA:5']),
       ...            atom_ids2=(a['CA:4'], a['C:4'], a['N:5'], a['CA:5']))
    """

    rsr.unpick(atom_ids1)
    rsr.add(forms.gaussian(group=physical.dihedral,
                           feature=features.dihedral(atom_ids1),
                           mean=3.141593, stdev=0.087))

    rsr.unpick(atom_ids2)
    rsr.add(forms.gaussian(group=physical.dihedral,
                           feature=features.dihedral(atom_ids2),
                           mean=0.0, stdev=0.087))
