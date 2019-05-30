"""Classes for handling alignment sequences (without structures)"""

import _modeller
from modeller import sequence

__docformat__ = "epytext en"

class Residue(sequence.Residue):
    """A single residue in an aligned sequence"""

    def get_position(self):
        """Get the position in the alignment of this residue"""
        invaln = _modeller.mod_alignment_invaln_get(self.mdl.aln.modpt)
        num = _modeller.mod_int2_get(invaln, self._num, self.mdl._num)
        return self.mdl.aln.positions[num - 1]

    def get_aligned_residue(self, seq):
        """Get the residue in C{seq} that is aligned with this one, or None"""
        alnpos = self.get_position()
        return alnpos.get_residue(seq)

    def add_leading_gaps(self, ngap=1):
        """Add gaps immediately before this residue."""
        return _modeller.mod_alnresidue_add_gaps(self.mdl.aln.modpt, self._num,
                                                 self.mdl._num + 1, ngap)

    def add_trailing_gaps(self, ngap=1):
        """Add gaps immediately after this residue."""
        return _modeller.mod_alnresidue_add_gaps(self.mdl.aln.modpt,
                                                 self._num + 1,
                                                 self.mdl._num + 1, ngap)

    def remove_leading_gaps(self, ngap=1):
        """Remove gaps immediately before this residue."""
        return _modeller.mod_alnresidue_remove_gaps(self.mdl.aln.modpt,
                                                    self._num,
                                                    self.mdl._num + 1, ngap)

    def remove_trailing_gaps(self, ngap=1):
        """Remove gaps immediately after this residue."""
        return _modeller.mod_alnresidue_remove_gaps(self.mdl.aln.modpt,
                                                    self._num + 1,
                                                    self.mdl._num + 1, ngap)

    def get_leading_gaps(self):
        """Get the number of gaps in the alignment immediately preceding this
           residue."""
        mypos = self.get_position().num
        if self._num > 0:
            prepos = self.mdl.residues[self._num - 1].get_position()
            prepos = prepos.num
            return mypos - prepos - 1
        else:
            return mypos

    def get_trailing_gaps(self):
        """Get the number of gaps in the alignment immediately following this
           residue."""
        mypos = self.get_position().num
        try:
            postpos = self.mdl.residues[self._num + 1].get_position()
            postpos = postpos.num
        except IndexError:
            postpos = self.mdl.naln
        return postpos - mypos - 1

    def __get_type(self):
        irestyp = _modeller.mod_sequence_irestyp_get(self.mdl.seqpt)
        return _modeller.mod_int1_get(irestyp, self._num)
    def __set_type(self, val):
        return _modeller.mod_alnresidue_type_set(self.mdl.aln.modpt, self._num,
                                                 self.mdl._num + 1, val,
                                                 self.mdl.aln.env.libs.modpt)
    type = property(__get_type, __set_type, doc="Integer residue type")


class Sequence(sequence.Sequence):
    """A single sequence (without structure) within an alignment"""
    aln = None
    _num = None
    env = None
    _residue_class = Residue

    def __init__(self, aln, num):
        self.aln = aln
        self.env = self.aln.env
        self._num = num
        sequence.Sequence.__init__(self)

    def __len__(self):
        return self.nres

    def __repr__(self):
        return "Sequence %s" % repr(self.code)
    def __str__(self):
        return "<%s>" % repr(self)

    def transfer_res_prop(self):
        """Transfer residue properties of predicted secondary structure"""
        return _modeller.mod_transfer_res_prop(self.aln.modpt, self._num)

    def get_num_equiv(self, seq):
        """Get the number of identical aligned residues between this sequence
           and C{seq}."""
        neqv = 0
        for res in self.residues:
            other_res = res.get_aligned_residue(seq)
            if other_res is not None and res.type == other_res.type:
                neqv += 1
        return neqv

    def get_sequence_identity(self, seq):
        """Get the % sequence identity between this sequence and C{seq}, defined
           as the number of identical aligned residues divided by the length
           of the shorter sequence."""
        return 100.0 * self.get_num_equiv(seq) / min(len(self), len(seq))

    def __get_naln(self):
        return _modeller.mod_alignment_naln_get(self.aln.modpt)
    def __get_code(self):
        return _modeller.mod_alignment_codes_get(self.aln.modpt, self._num)
    def __set_code(self, val):
        _modeller.mod_alignment_codes_set(self.aln.modpt, self._num, val)
    def __get_atom_file(self):
        alnseq = self.__get_alnseqpt()
        return _modeller.mod_alnsequence_atom_files_get(alnseq)
    def __set_atom_file(self, val):
        alnseq = self.__get_alnseqpt()
        _modeller.mod_alnsequence_atom_files_set(alnseq, val)
    def __get_seqpt(self):
        return _modeller.mod_alignment_sequence_get(self.aln.modpt, self._num)
    def __get_alnseqpt(self):
        return _modeller.mod_alignment_alnsequence_get(self.aln.modpt,
                                                       self._num)

    code = property(__get_code, __set_code, doc="Alignment code")
    atom_file = property(__get_atom_file, __set_atom_file, doc="PDB file name")
    naln = property(__get_naln, doc="Length of alignment (including gaps)")
    seqpt = property(__get_seqpt)
