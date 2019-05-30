"""Classes for handling databases of protein sequences"""

import _modeller
import util.modutil as modutil
import residue
from modeller.util.modobject import modobject

__docformat__ = "epytext en"

class sequence_db(modobject):
    """Holds a database of protein sequences"""
    __modpt = None
    env = None

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_sequence_db_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_sequence_db_new(self)

    def __del__(self):
        if self.__modpt:
            _modeller.mod_sequence_db_free(self.__modpt)

    def __repr__(self):
        if len(self) == 1:
            return "Database of 1 sequence"
        else:
            return "Database of %d sequences" % len(self)

    def __str__(self):
        return "<%s>" % repr(self)

    def __get_modpt(self):
        return self.__modpt

    def __len__(self):
        return _modeller.mod_sequence_db_nchn_get(self.modpt)

    def read(self, chains_list, seq_database_file, seq_database_format,
             clean_sequences=True, minmax_db_seq_len=(0, 999999)):
        """Reads in a database from a file"""
        return _modeller.mod_sequence_db_read(self.modpt, self.env.libs.modpt,
                                              chains_list, seq_database_file,
                                              seq_database_format,
                                              clean_sequences,
                                              minmax_db_seq_len)

    def convert(self, chains_list, seq_database_file, seq_database_format,
                outfile, clean_sequences=True, minmax_db_seq_len=(0, 999999)):
        """Converts a database to binary format"""
        return _modeller.mod_sequence_db_convert(self.modpt,
                                                 self.env.libs.modpt,
                                                 chains_list, seq_database_file,
                                                 seq_database_format, outfile,
                                                 clean_sequences,
                                                 minmax_db_seq_len)

    def write(self, chains_list, seq_database_file, seq_database_format):
        """Writes out a database to a file"""
        return _modeller.mod_sequence_db_write(self.modpt, self.env.libs.modpt,
                                               chains_list, seq_database_file,
                                               seq_database_format)

    def search(self, aln, seq_database_file, search_group_list,
               search_randomizations=0, search_top_list=20, off_diagonal=100,
               overhang=0, gap_penalties_1d=(-900., -50.),
               signif_cutoff=(4.0, 5.0), rr_file='$(LIB)/as1.sim.mat',
               matrix_offset=0., fast_search_cutoff=1.0, data_file=False,
               search_sort='LONGER', output='LONG',
               alignment_features='INDICES CONSERVATION',
               local_alignment=False, fast_search=False, io=None, **vars):
        """Search for similar sequences"""
        if io is None:
            io = self.env.io
        func = _modeller.mod_sequence_db_search
        return func(self.modpt, aln.modpt, io.modpt, self.env.libs.modpt,
                    search_randomizations, search_top_list, off_diagonal,
                    overhang, gap_penalties_1d, signif_cutoff, rr_file,
                    matrix_offset, fast_search_cutoff, data_file,
                    search_group_list, search_sort, output, alignment_features,
                    seq_database_file, local_alignment, fast_search)

    def filter(self, seqid_cut, output_grp_file, output_cod_file,
               gap_penalties_1d=(-900., -50.), matrix_offset=0.,
               rr_file='$(LIB)/as1.sim.mat', max_diff_res=30):
        """Cluster sequences by sequence-identity"""
        return _modeller.mod_sequence_db_filter(self.modpt, self.env.libs.modpt,
                                                gap_penalties_1d, matrix_offset,
                                                rr_file, seqid_cut,
                                                max_diff_res, output_grp_file,
                                                output_cod_file)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return Sequence(self, indx)
        else:
            return [self[ind] for ind in ret]

    modpt = property(__get_modpt)


class Sequence(object):
    """A single sequence in the database"""

    def __init__(self, sdb, num):
        self.sdb = sdb
        self.env = sdb.env
        self.num = num

    def __len__(self):
        nseq = _modeller.mod_sequence_db_nseqdb_get(self.sdb.modpt)
        return _modeller.mod_int1_get(nseq, self.num)

    def __repr__(self):
        if len(self) == 1:
            return "Sequence of 1 residue"
        else:
            return "Sequence of %d residues" % len(self)
    def __str__(self):
        return "<%s>" % repr(self)

    def __get_code(self):
        return _modeller.mod_sequence_db_code_get(self.sdb.modpt, self.num)
    def __set_code(self, val):
        _modeller.mod_sequence_db_code_set(self.sdb.modpt, self.num, val)
    def __get_prottyp(self):
        return _modeller.mod_sequence_db_prottyp_get(self.sdb.modpt, self.num)
    def __set_prottyp(self, val):
        _modeller.mod_sequence_db_prottyp_set(self.sdb.modpt, self.num, val)
    def __get_resol(self):
        resoldb = _modeller.mod_sequence_db_resoldb_get(self.sdb.modpt)
        return _modeller.mod_float1_get(resoldb, self.num)
    def __set_resol(self, val):
        resoldb = _modeller.mod_sequence_db_resoldb_get(self.sdb.modpt)
        _modeller.mod_float1_set(resoldb, self.num, val)
    def __get_residues(self):
        return ResidueList(self)

    code = property(__get_code, __set_code, doc="Alignment code")
    prottyp = property(__get_prottyp, __set_prottyp, doc="Protein type")
    resol = property(__get_resol, __set_resol, doc="Resolution")
    residues = property(__get_residues, doc="All residues in this sequence""")


class ResidueList(object):
    """A list of residues (L{Residue} objects) in a single sequence in
       a sequence database."""

    def __init__(self, seq, offset=0, length=None):
        self.seq = seq
        self.offset = offset
        self.length = length

    def __len__(self):
        if self.length is not None:
            return self.length
        else:
            return len(self.seq)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return Residue(self.seq, ret + self.offset)
        else:
            return [self[ind] for ind in ret]


class Residue(residue.Residue):
    """A single residue in a sequence database."""

    def __get_type(self):
        return _modeller.mod_sequence_db_restype_get(self.mdl.sdb.modpt,
                                                     self.mdl.num, self._num)
    def __set_type(self, val):
        _modeller.mod_sequence_db_restype_set(self.mdl.sdb.modpt,
                                              self.mdl.num, self._num, val)
    type = property(__get_type, __set_type, doc="Integer residue type")
