import _modeller
from modeller.restraints import Restraints
import model_topology
import alignment
from modeller.util.modobject import modobject
from modeller import modfile, coordinates, alnsequence
from modeller.atom_type import AtomType

__docformat__ = "epytext en"

class Atom(coordinates.Atom):
    """A single atom in a protein structure"""

    def __get_dvx(self):
        dvx = _modeller.mod_model_dvx_get(self.mdl.modpt)
        return _modeller.mod_float1_get(dvx, self._num)
    def __get_dvy(self):
        dvy = _modeller.mod_model_dvy_get(self.mdl.modpt)
        return _modeller.mod_float1_get(dvy, self._num)
    def __get_dvz(self):
        dvz = _modeller.mod_model_dvz_get(self.mdl.modpt)
        return _modeller.mod_float1_get(dvz, self._num)
    def __get_vx(self):
        vx = _modeller.mod_model_vx_get(self.mdl.modpt)
        return _modeller.mod_float1_get(vx, self._num)
    def __get_vy(self):
        vy = _modeller.mod_model_vy_get(self.mdl.modpt)
        return _modeller.mod_float1_get(vy, self._num)
    def __get_vz(self):
        vz = _modeller.mod_model_vz_get(self.mdl.modpt)
        return _modeller.mod_float1_get(vz, self._num)
    def __get_mass(self):
        return self.type.mass
    def __get_charge(self):
        charge = _modeller.mod_model_charge_get(self.mdl.modpt)
        return _modeller.mod_float1_get(charge, self._num)
    def __set_charge(self, val):
        charge = _modeller.mod_model_charge_get(self.mdl.modpt)
        _modeller.mod_float1_set(charge, self._num, val)
    def __get_type(self):
        iattyp = _modeller.mod_model_iattyp_get(self.mdl.modpt)
        return AtomType(self.mdl, _modeller.mod_int1_get(iattyp, self._num)-1)
    def __get_gprsr_class(self):
        iatta = _modeller.mod_model_iatta_get(self.mdl.modpt)
        return _modeller.mod_int1_get(iatta, self._num)

    dvx = property(__get_dvx, doc="Objective function derivative, dF/dx")
    dvy = property(__get_dvy, doc="Objective function derivative, dF/dy")
    dvz = property(__get_dvz, doc="Objective function derivative, dF/dz")
    vx = property(__get_vx, doc="x component of velocity")
    vy = property(__get_vy, doc="y component of velocity")
    vz = property(__get_vz, doc="z component of velocity")
    mass = property(__get_mass, doc="Atomic mass")
    charge = property(__get_charge, __set_charge, doc="Electrostatic charge")
    type = property(__get_type, doc="Atom type")
    gprsr_class = property(__get_gprsr_class, doc="group_restraints class")


class model(coordinates.Coordinates):
    """Holds a model of a protein"""
    __modpt = None
    env = None
    __gprsr = None
    _atom_class = Atom

    def __init__(self, env, **vars):
        coordinates.Coordinates.__init__(self)
        self.__modpt = _modeller.mod_model_new(self)
        self.env = env.copy()
        self.group_restraints = self.env.group_restraints
        if len(vars) > 0:
            self.read(**vars)

    def __getstate__(self):
        d = coordinates.Coordinates.__getstate__(self)
        # Don't pickle cached information
        for key in ('dope_restraints', 'dopehr_restraints'):
            if key in d:
                del d[key]
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_model_new(self)
        # Ensure that group_restraints are loaded at the Fortran level
        self.group_restraints = self.__gprsr

    def __del__(self):
        if self.modpt:
            _modeller.mod_model_free(self.modpt)

    def __repr__(self):
        return "Model containing %s, %s, and %s" \
               % (repr(self.chains), repr(self.residues), repr(self.atoms))

    def __str__(self):
        return "<%s>" % repr(self)

    def __get_modpt(self):
        return self.__modpt
    def __get_seqpt(self):
        return _modeller.mod_model_seq_get(self.__modpt)
    def __get_cdpt(self):
        return _modeller.mod_model_cd_get(self.__modpt)

    def get_atom_indices(self):
        """Get the indices of all atoms in this model"""
        return (range(1, self.natm+1), self)

    def read(self, file, model_format='PDB',
             model_segment=('FIRST:@', 'LAST:'), io=None):
        """Read coordinates from a file"""
        if io is None:
            io = self.env.io
        return _modeller.mod_model_read(self.modpt, io.modpt,
                                        self.env.libs.modpt, file,
                                        model_format, model_segment)

    def write(self, file, model_format='PDB', no_ter=False):
        """Write coordinates to a file"""
        if isinstance(file, str):
            file = modfile.File(file, 'w')
        return _modeller.mod_model_write(self.modpt, self.env.libs.modpt, (),
                                         file.file_pointer, model_format,
                                         no_ter, True)

    def build_sequence(self, sequence, special_patches=None,
                       patch_default=True):
        """Build an extended chain from a string of one-letter residue codes"""
        a = alignment.alignment(self.env)
        a.append_sequence(sequence)
        self.clear_topology()
        self.generate_topology(a[0], patch_default=patch_default)
        if special_patches:
            special_patches(self)
        self.build(initialize_xyz=True, build_method='INTERNAL_COORDINATES')
        self.prottyp = 'structure'

    def use_lennard_jones(self):
        """Set up to use Lennard Jones rather than soft sphere (the default) for
           van der Waals interactions"""
        edat = self.env.edat
        edat.contact_shell = 8.00
        edat.dynamic_sphere = False
        edat.dynamic_lennard = True

    def assess_ga341(self):
        """Assess with the GA341 method"""
        return _modeller.mod_assess_ga341(self.modpt, self.env.libs.modpt)

    def fast_rmsd(self, mdl):
        """Calculate the RMSD between this model and the input."""
        return _modeller.mod_model_fast_rmsd(self.modpt, mdl.modpt)

    def orient(self):
        """Center and orient the model"""
        import orient
        retval = _modeller.mod_model_orient(self.modpt)
        return orient.OrientData(*retval)

    def build(self, build_method, initialize_xyz):
        """Build coordinates from topology"""
        return _modeller.mod_model_build(self.modpt, self.env.libs.modpt,
                                         build_method, initialize_xyz)

    def build_ic(self, initialize_xyz):
        """Build coordinates from internal coordinates"""
        return _modeller.mod_model_build_ic(self.modpt, self.env.libs.modpt,
                                            initialize_xyz)

    def saxs_intens(self, saxsd, filename, fitflag=False):
        """Calculate SAXS intensity from model."""
        return _modeller.mod_saxs_intens(self.modpt, saxsd.modpt, filename,
                                         fitflag)

    def saxs_chifun(self, transfer_is, edat=None):
        """Calculate SAXS score from model.
           @param transfer_is: if False: I(s) is not computed, i.e. the
                  intensity which is in memory is used
        """
        if edat is None:
            edat = self.env.edat
        return _modeller.mod_saxs_chifun(edat.modpt, self.modpt, transfer_is)

    def saxs_pr(self, saxsd, filename):
        """Calculate P(r) from model
           mdl.saxs_pr(saxsd, filename, sample_dr)

           mdl          model
           sample_dr    sampling of P(r) in Angstrom"""
        return _modeller.mod_saxs_pr(self.modpt, saxsd.modpt, filename)

    def transfer_xyz(self, aln, cluster_cut=-1.0, cluster_method='RMSD',
                     io=None):
        """Copy coordinates from template structures"""
        if io is None:
            io = self.env.io
        return _modeller.mod_transfer_xyz(self.modpt, aln.modpt, io.modpt,
                                          self.env.libs.modpt, cluster_cut,
                                          cluster_method)


    def res_num_from(self, mdl, aln):
        """Copy residue numbers from the given model"""
        return _modeller.mod_model_res_num_from(self.modpt, mdl.modpt,
                                                aln.modpt, self.env.libs.modpt)

    def reorder_atoms(self):
        """Standardize atom order to match the current topology library"""
        return _modeller.mod_model_reorder_atoms(self.modpt,
                                                 self.env.libs.modpt)

    def rename_segments(self, segment_ids, renumber_residues=[]):
        """Relabel residue numbers in each chain/segment"""
        return _modeller.mod_model_rename_segments(self.modpt, segment_ids,
                                                   renumber_residues)

    def to_iupac(self):
        """Make dihedral angles satisfy the IUPAC convention"""
        return _modeller.mod_model_to_iupac(self.modpt, self.env.libs.modpt)

    def clear_topology(self):
        """Clear all covalent topology (atomic connectivitiy) and sequence"""
        return _modeller.mod_model_topology_clear(self.modpt)

    def generate_topology(self, alnseq, patch_default=None, io=None):
        """Generate covalent topology (atomic connectivity) and sequence"""
        if not isinstance(alnseq, alnsequence.Sequence):
            raise TypeError("""Must use an alignment 'Sequence' object here.
For example, replace 'generate_topology(aln, sequence="foo")' with
'generate_topology(aln["foo"])'""")
        if io is None:
            io = self.env.io
        if patch_default is None:
            patch_default = self.env.patch_default
        func = _modeller.mod_model_topology_generate
        return func(self.modpt, alnseq.aln.modpt, io.modpt,
                    self.env.libs.modpt, alnseq._num, patch_default)

    def patch(self, residue_type, residues):
        """Patch the model topology"""
        if not hasattr(residues, '__iter__'):
            residues = [residues]
        for r in residues:
            if not isinstance(r, coordinates.Residue) or r.mdl is not self:
                raise TypeError("expecting one or more 'Residue' objects")
        return _modeller.mod_model_patch(mdl=self.modpt,
                                         libs=self.env.libs.modpt,
                                         residue_type=residue_type,
                                         residue_ids=[r.index \
                                                      for r in residues])

    def patch_ss(self):
        """Guess disulfides from the current structure"""
        return _modeller.mod_model_patch_ss(self.modpt, self.env.libs.modpt)

    def patch_ss_templates(self, aln, io=None):
        """Guess disulfides from templates"""
        if io is None:
            io = self.env.io
        return _modeller.mod_model_patch_ss_templates(self.modpt, aln.modpt,
                                                      io.modpt,
                                                      self.env.libs.modpt)

    def write_data(self, output, file=None, surftyp=1, neighbor_cutoff=6.0,
                   accessibility_type=8, probe_radius=1.4,
                   psa_integration_step=0.1,
                   dnr_accpt_lib='${LIB}/donor_acceptor.lib', edat=None):
        """Write derivative model data"""
        if edat is None:
            edat = self.env.edat
        return _modeller.mod_model_write_data(self.modpt, edat.modpt,
                                              self.env.libs.modpt, surftyp,
                                              neighbor_cutoff,
                                              accessibility_type, output, file,
                                              probe_radius,
                                              psa_integration_step,
                                              dnr_accpt_lib)

    def make_region(self, atom_accessibility=1.0, region_size=20):
        """Define a random surface patch of atoms"""
        return _modeller.mod_model_make_region(self.modpt, self.env.libs.modpt,
                                               atom_accessibility, region_size)

    def color(self, aln):
        """Color according to the alignment"""
        return _modeller.mod_model_color(self.modpt, aln.modpt)

    def loops(self, aln, minlength, maxlength, insertion_ext, deletion_ext):
        """Returns a list of all loops (insertions or deletions) in the model,
           as defined by the alignment"""
        return self.get_insertions(aln, minlength, maxlength, insertion_ext) \
               + self.get_deletions(aln, deletion_ext)

    def get_insertions(self, aln, minlength, maxlength, extension):
        """Returns a list of all insertions in the model, as defined by the
           alignment."""
        return self.__get_insdel(aln, _modeller.mod_alignment_next_insert,
                                 minlength, maxlength, extension)

    def get_deletions(self, aln, extension):
        """Returns a list of all deletions in the model, as defined by the
           alignment."""
        return self.__get_insdel(aln, _modeller.mod_alignment_next_delete,
                                 extension)

    def assess_normalized_dope(self):
        """Assess the model, and return a normalized DOPE score (z score)"""
        from modeller.selection import selection
        import normalized_dope
        sel = selection(self)
        dope_score = sel.assess_dope()
        scorer = normalized_dope.DOPEScorer(self)
        z_score = scorer.get_z_score(dope_score)
        print ">> Normalized DOPE z score: %.3f" % z_score
        return z_score

    def assess_normalized_dopehr(self):
        """Assess the model, and return a normalized DOPE-HR score (z score)"""
        from modeller.selection import selection
        import normalized_dope
        sel = selection(self)
        dope_score = sel.assess_dopehr()
        scorer = normalized_dope.DOPEHRScorer(self)
        z_score = scorer.get_z_score(dope_score)
        print ">> Normalized DOPE-HR z score: %.3f" % z_score
        return z_score

    def get_normalized_dope_profile(self):
        """Return a normalized DOPE per-residue profile"""
        from modeller.selection import selection
        import normalized_dope
        import physical
        sel = selection(self)
        edat = sel.get_dope_energy_data()
        oldgprsr = self.group_restraints
        self.group_restraints = sel.get_dope_potential()
        try:
            profile = sel.get_energy_profile(edat, physical.nonbond_spline)
        finally:
            self.group_restraints = oldgprsr
        scorer = normalized_dope.DOPEScorer(self)
        return scorer.get_profile(profile)

    def write_psf(self, file, xplor=True):
        """Write the molecular topology to a PSF file, in either X-PLOR format
           (the default) or CHARMM format."""
        model_topology.write_psf(file, self, xplor)

    def find_atoms(self, residue_type, atom_names):
        return model_topology.FindAtoms(self, residue_type, atom_names)

    def find_chi1_dihedrals(self, residue_type):
        return model_topology.FindDihedrals(self, residue_type, 5)
    def find_chi2_dihedrals(self, residue_type):
        return model_topology.FindDihedrals(self, residue_type, 6)
    def find_chi3_dihedrals(self, residue_type):
        return model_topology.FindDihedrals(self, residue_type, 7)
    def find_chi4_dihedrals(self, residue_type):
        return model_topology.FindDihedrals(self, residue_type, 8)

    def __get_insdel(self, aln, func, *args):
        _modeller.mod_alnsequence_check_model(aln.modpt, len(aln) - 1,
                                              self.modpt, self.env.libs.modpt)
        l = []
        pos = 0
        while pos >= 0:
            (pos, start, end) = func(aln.modpt, pos, *args)
            if start > 0 and end > 0:
                l.append(coordinates.ResidueList(self, start - 1,
                                                 end - start + 1))
        return l

    def __get_seq_id(self):
        return _modeller.mod_model_seq_id_get(self.modpt)
    def __set_seq_id(self, val):
        return _modeller.mod_model_seq_id_set(self.modpt, val)
    def __get_remark(self):
        return _modeller.mod_model_remark_get(self.modpt)
    def __set_remark(self, val):
        return _modeller.mod_model_remark_set(self.modpt, val)
    def __get_header(self):
        return _modeller.mod_model_header_get(self.modpt)
    def __set_header(self, val):
        return _modeller.mod_model_header_set(self.modpt, val)
    def __get_last_energy(self):
        return _modeller.mod_model_last_energy_get(self.modpt)
    def __set_last_energy(self, val):
        return _modeller.mod_model_last_energy_set(self.modpt, val)
    def __get_restraints(self):
        return Restraints(self)
    def __set_group_restraints(self, val):
        self.__gprsr = val
        if val:
            _modeller.mod_model_gprsr_set(self.modpt, self.env.libs.modpt,
                                          val._modpt)
        else:
            _modeller.mod_model_gprsr_unset(self.modpt, self.env.libs.modpt)
    def __get_group_restraints(self):
        return self.__gprsr
    def __get_bonds(self):
        return model_topology.BondList(self,
                                       _modeller.mod_model_topology_nbnd_get,
                                       _modeller.mod_model_topology_iatb_get, 2)
    def __get_angles(self):
        return model_topology.BondList(self,
                                       _modeller.mod_model_topology_nang_get,
                                       _modeller.mod_model_topology_iata_get, 3)
    def __get_dihedrals(self):
        return model_topology.BondList(self,
                                       _modeller.mod_model_topology_ndih_get,
                                       _modeller.mod_model_topology_iatd_get, 4)
    def __get_impropers(self):
        return model_topology.BondList(self,
                                       _modeller.mod_model_topology_nimp_get,
                                       _modeller.mod_model_topology_iati_get, 4)


    modpt = property(__get_modpt, doc="Internal Modeller object")
    cdpt = property(__get_cdpt, doc="Internal Modeller coordinates object")
    seqpt = property(__get_seqpt, doc="Internal Modeller sequence object")
    seq_id = property(__get_seq_id, __set_seq_id,
                      doc="Sequence identity between model and best template")
    header = property(__get_header, __set_header, doc="PDB header")
    last_energy = property(__get_last_energy, __set_last_energy,
                           doc="Energy from last energy or optimize")
    remark = property(__get_remark, __set_remark, doc="PDB REMARK line(s)")
    restraints = property(__get_restraints,
                          doc="All restraints acting on this model")
    group_restraints = property(__get_group_restraints, __set_group_restraints,
                                doc="Group restraints active for this model")
    bonds = property(__get_bonds, doc="All defined bonds")
    angles = property(__get_angles, doc="All defined angles")
    dihedrals = property(__get_dihedrals, doc="All defined dihedrals")
    impropers = property(__get_impropers, doc="All defined impropers")
