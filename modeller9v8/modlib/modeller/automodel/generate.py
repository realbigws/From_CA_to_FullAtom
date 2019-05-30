"""Functions to build a structure from a sequence in the alignment file.
   Usually called by setting automodel.generate_method."""

from modeller.alignment import alignment

def generate_xyz(mdl, aln):
    """Build coordinates from scratch (ignore templates)"""

    mdl.read_top_par()
    mdl.create_topology(aln)
    mdl.build(initialize_xyz=True, build_method='3D_INTERPOLATION')


def transfer_xyz(mdl, aln):
    """Build structure by copying equivalent coordinates from the templates"""

    # If initial malign3d was requested, orient the template structures but
    # then restore the original alignment
    if mdl.initial_malign3d:
        aln.clear()
        aln.append(file=mdl.alnfile, align_codes=mdl.knowns)
        aln.malign3d(fit=False, gap_penalties_3d=(0, 4))
        mdl.read_alignment(aln)
    mdl.read_top_par()
    mdl.create_topology(aln)
    mdl.transfer_xyz(aln, cluster_cut=-1.0)
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')


def read_xyz(mdl, aln):
    """Read in the initial structure from an existing file"""

    # Create two copies of the template sequence:
    a = alignment(mdl.env)
    a.append(file=mdl.alnfile, align_codes=[mdl.sequence]*2)

    # Use the initial model as the first structure in the alignment:
    a[0].prottyp = 'structureX'
    a[0].atom_file = a[0].code = mdl.my_inifile

    mdl.read_top_par()
    mdl.create_topology(aln)

    # Get model coordinates from the initial model, making sure that the
    # sequence is correct and any missing atoms are filled in:
    mdl.transfer_xyz(a, cluster_cut=-1.0)
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
