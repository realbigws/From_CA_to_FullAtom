"""Methods to perturb the initial model. Usually called by setting
   automodel.rand_method."""

from modeller.selection import selection

def xyz(atmsel):
    """Randomize all coordinates of residues with defined topology.
       (BLK residues are not randomized since they are generally kept rigid.)"""
    mdl = atmsel.get_model()
    atmsel = atmsel - atmsel.only_no_topology()
    if len(atmsel) > 0:
        atmsel.randomize_xyz(deviation=mdl.deviation)

def dihedrals(atmsel):
    """Randomize dihedrals"""
    mdl = atmsel.get_model()
    aln = alignment(mdl.env, file=mdl.alnfile,
                    align_codes=mdl.knowns+[mdl.sequence])
    # Just in case, generate topology again (ROTATE needs bonds)
    # (could replace with GENERATE_TOPOLOGY if no special patches)
    mdl.create_topology(aln)

    # Optimize all dihedral angles:
    atmsel.rotate_dihedrals(change='RANDOMIZE', deviation=mdl.deviation,
                            dihedrals=('phi', 'psi', 'omega', 'chi1', 'chi2',
                                       'chi3', 'chi4'))
