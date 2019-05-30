import _modeller
import util.modutil as modutil
import util.modlist as modlist
from modeller.util.modobject import modobject
import alignment
import pssmdb
import residue

__docformat__ = "epytext en"

class profile(modobject):
    """Holds a profile of multiple sequences"""
    __modpt = None
    env = None

    def __init__(self, env, aln=None, **vars):
        self.__modpt = _modeller.mod_profile_new(self)
        self.env = env.copy()
        if aln is not None:
            _modeller.mod_profile_from_aln(aln=aln.modpt, prf=self.modpt,
                                           libs=self.env.libs.modpt)
        elif len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_profile_new(self)

    def __del__(self):
        if self.__modpt:
            _modeller.mod_profile_free(self.__modpt)

    def __get_modpt(self):
        return self.__modpt

    def __len__(self):
        return _modeller.mod_profile_nseq_get(self.modpt)

    def __opt_file(self, fname):
        """Processing for optional files for some routines"""
        if fname is None or fname == "":
            return ('', False)
        else:
            return (fname, True)

    def read(self, file, profile_format):
        """Reads a profile from a specified file"""
        return _modeller.mod_profile_read(self.modpt, self.env.libs.modpt, file,
                                          profile_format)

    def write(self, file, profile_format):
        """Write profile to a file"""
        return _modeller.mod_profile_write(self.modpt, self.env.libs.modpt,
                                           file, profile_format)

    def build(self, sdb, gap_penalties_1d=(-900., -50.), matrix_offset=0.,
              rr_file='$(LIB)/as1.sim.mat', n_prof_iterations=3,
              max_aln_evalue=0.1, matrix_scaling_factor=0.0069,
              check_profile=True, output_score_file=None, gaps_in_target=False,
              score_statistics=True, pssm_weights_type='HH1', pssm_file=None):
        """Refines the profile with sequences from the given database"""
        (output_score_file, output_scores) = self.__opt_file(output_score_file)
        (pssm_file, write_pssm) = self.__opt_file(pssm_file)
        func = _modeller.mod_profile_build
        return func(self.modpt, sdb.modpt, self.env.libs.modpt,
                    gap_penalties_1d, matrix_offset, rr_file,
                    n_prof_iterations, max_aln_evalue, matrix_scaling_factor,
                    check_profile, output_scores, output_score_file,
                    gaps_in_target, score_statistics, pssm_weights_type,
                    write_pssm, pssm_file)

    def to_alignment(self):
        """Converts the profile to alignment format"""
        aln = alignment.alignment(self.env)
        _modeller.mod_profile_to_aln(prf=self.modpt, aln=aln.modpt,
                                     libs=self.env.libs.modpt)
        return aln

    def scan(self, profile_list_file, matrix_offset=0., profile_format='TEXT',
             rr_file='$(LIB)/as1.sim.mat', gap_penalties_1d=(-900., -50.),
             matrix_scaling_factor=0.0069, max_aln_evalue=0.1,
             aln_base_filename='alignment', score_statistics=True,
             output_alignments=True, output_score_file=None,
             pssm_weights_type='HH1', summary_file='ppscan.sum',
             ccmatrix_offset=-200, score_type='CCMAT', psm=None):
        """Compare the profile against a database of profiles"""
        (output_score_file, output_scores) = self.__opt_file(output_score_file)
        (summary_file, write_summary) = self.__opt_file(summary_file)
        if psm is None:
            psm = pssmdb.pssmdb(self.env)
        func = _modeller.mod_profile_scan
        return func(self.modpt, psm.modpt, self.env.libs.modpt,
                    profile_list_file, matrix_offset, profile_format, rr_file,
                    gap_penalties_1d, matrix_scaling_factor, max_aln_evalue,
                    aln_base_filename, score_statistics, output_alignments,
                    output_scores, output_score_file, pssm_weights_type,
                    write_summary, summary_file, ccmatrix_offset, score_type)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return Sequence(self, indx)
        else:
            return [self[ind] for ind in ret]

    def __get_positions(self):
        return PositionList(self)
    def __get_filename(self):
        return _modeller.mod_profile_filename_get(self.modpt)

    modpt = property(__get_modpt)
    filename = property(__get_filename,
                        doc="Name of the file this profile was read from")
    positions = property(__get_positions, doc="Profile alignment positions")


class Sequence(object):
    """A single sequence in a profile"""

    def __init__(self, prf, num):
        self._prf = prf
        self._num = num

    def __len__(self):
        nres = _modeller.mod_profile_nres_get(self._prf.modpt)
        return _modeller.mod_int1_get(nres, self._num)

    def __get_code(self):
        return _modeller.mod_profile_code_get(self._prf.modpt, self._num)
    def __set_code(self, val):
        _modeller.mod_profile_code_set(self._prf.modpt, self._num, val)
    def __get_prottyp(self):
        return _modeller.mod_profile_prottyp_get(self._prf.modpt, self._num)
    def __set_prottyp(self, val):
        _modeller.mod_profile_prottyp_set(self._prf.modpt, self._num, val)
    def __get_iter(self):
        iter = _modeller.mod_profile_iter_get(self._prf.modpt)
        return _modeller.mod_int1_get(iter, self._num)
    def __get_neqv(self):
        neqv = _modeller.mod_profile_neqv_get(self._prf.modpt)
        return _modeller.mod_int1_get(neqv, self._num)
    def __get_fid(self):
        fid = _modeller.mod_profile_fid_get(self._prf.modpt)
        return _modeller.mod_float1_get(fid, self._num)
    def __get_evalue(self):
        evalue = _modeller.mod_profile_evalue_get(self._prf.modpt)
        return _modeller.mod_double1_get(evalue, self._num)

    code = property(__get_code, __set_code, doc="Code")
    prottyp = property(__get_prottyp, __set_prottyp,
                       doc="Protein type (sequence/structure")
    iter = property(__get_iter, doc="Number of iterations")
    neqv = property(__get_neqv, doc="Number of equivalences")
    fid = property(__get_fid)
    evalue = property(__get_evalue, doc="E value")


class PositionList(modlist.FixList):
    """A list of L{Position} objects."""

    def __init__(self, prf):
        self.__prf = prf
        modlist.FixList.__init__(self)

    def __len__(self):
        return _modeller.mod_profile_npos_get(self.__prf.modpt)

    def _getfunc(self, indx):
        return Position(self.__prf, indx)

class Position(object):
    """A profile position"""

    def __init__(self, prf, indx):
        self.__prf = prf
        self.__indx = indx

    def get_residue(self, seq):
        """Get the residue in C{seq} that is at this profile position, or None
           if a gap is present."""
        prf = self.__prf
        if not isinstance(seq, Sequence):
            raise TypeError("Expected a profile 'Sequence' object for seq")
        if seq._prf != prf:
            raise ValueError("seq must be a sequence in the same profile")
        sprofile = _modeller.mod_profile_sprofile_get(prf.modpt)
        irestyp = _modeller.mod_int2_get(sprofile, self.__indx, seq._num)
        if irestyp == 21:   # 21 = gap
            return None
        else:
            return Residue(prf, 0, irestyp)

class Residue(residue.Residue):
    """A single residue in a profile sequence"""

    def __init__(self, prf, num, type):
        residue.Residue.__init__(self, prf, num)
        self.type = type
