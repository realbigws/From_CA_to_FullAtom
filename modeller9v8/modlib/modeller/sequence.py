"""Classes for handling residue sequences, either in alignments or in
   protein structures."""

import _modeller
import residue
from modeller.util.modobject import modobject
import modeller.util.modutil as modutil
from modeller.util import modlist

__docformat__ = "epytext en"

class ChainList(object):
    """A list of chain objects."""

    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return self.seq.nseg

    def __repr__(self):
        ln = len(self)
        s = "%d chain" % ln
        if ln != 1:
            s += 's'
        return s

    def __str__(self):
        return "<List of " + repr(self) + ">"

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx, self.seq._indxseg, (0, None))
        if isinstance(ret, int):
            return self.seq._chain_class(self.seq, ret)
        else:
            return [self[ind] for ind in ret]


class Chain(object):
    """A single chain in a sequence"""

    def __init__(self, seq, num):
        self.seq = seq
        self._num = num

    def __repr__(self):
        return "<Chain %s>" % repr(self.name)

    def filter(self, structure_types='structure', minimal_resolution=99.0,
               minimal_chain_length=30, max_nonstdres=10,
               chop_nonstd_termini=True, minimal_stdres=30):
        """Does this chain pass all filter criteria?"""
        f = _modeller.mod_chain_filter
        return f(self.seq.seqpt, self._num, structure_types, minimal_resolution,
                 minimal_chain_length, max_nonstdres, chop_nonstd_termini,
                 minimal_stdres)

    def atom_file_and_code(self, filename):
        """Return suitable atom_file and align_codes for this chain, given
           a model filename."""
        return _modeller.mod_chain_atom_file_and_code(self.seq.seqpt, self._num,
                                                      filename)

    def join(self, chain):
        if chain.seq is not self.seq:
            raise ValueError("Chains must be from the same sequence")
        if chain._num <= self._num:
            raise ValueError("Chain to be joined must come after this chain")
        _modeller.mod_chain_join(self.seq.seqpt, self._num, chain._num)

    def _get_resind(self):
        iress1 = _modeller.mod_sequence_iress1_get(self.seq.seqpt)
        iress2 = _modeller.mod_sequence_iress2_get(self.seq.seqpt)
        return (_modeller.mod_int1_get(iress1, self._num) - 1,
                _modeller.mod_int1_get(iress2, self._num))

    def __get_residues(self):
        (startres, endres) = self._get_resind()
        suffix = ":%s" % self.name
        return self.seq._residue_list_class(self.seq, startres,
                                            endres - startres, suffix)
    def __get_name(self):
        return _modeller.mod_sequence_segid_get(self.seq.seqpt, self._num)
    def __set_name(self, val):
        return _modeller.mod_sequence_segid_set(self.seq.seqpt, self._num, val)
    def __get_index(self):
        return self._num + 1

    residues = property(__get_residues,
                        doc="List of all residues in this chain")
    name = property(__get_name, __set_name, doc="Chain ID")
    index = property(__get_index, doc="Internal residue index")


class SeqRange(modlist.FixList):
    """The starting and ending residue numbers and chain IDs for a sequence"""

    def __init__(self, seq):
        self.__seq = seq
        modlist.FixList.__init__(self)

    def __len__(self):
        return 2

    def _getfunc(self, indx):
        return _modeller.mod_sequence_rng_get(self.__seq.seqpt, indx)

    def _setfunc(self, indx, val):
        _modeller.mod_sequence_rng_set(self.__seq.seqpt, indx, val)


class Residue(residue.Residue):
    """A single residue in a sequence"""

    def __get_type(self):
        irestyp = _modeller.mod_sequence_irestyp_get(self.mdl.seqpt)
        return _modeller.mod_int1_get(irestyp, self._num)
    def __set_type(self, val):
        irestyp = _modeller.mod_sequence_irestyp_get(self.mdl.seqpt)
        _modeller.mod_int1_set(irestyp, self._num, val)
    def __get_chain(self):
        ichain = _modeller.mod_sequence_chain_for_res(self.mdl.seqpt, self._num)
        return self.mdl.chains[ichain]
    def __get_index(self):
        return self._num + 1

    type = property(__get_type, __set_type, doc="Integer residue type")
    chain = property(__get_chain, doc="Chain containing this residue")
    index = property(__get_index, doc="Internal residue index")


class ResidueList(object):
    """A list of residue objects."""

    def __init__(self, seq, offset=0, length=None, suffix=""):
        self.seq = seq
        self.offset = offset
        self.length = length
        self.suffix = suffix

    def __len__(self):
        if self.length is not None:
            return self.length
        else:
            return self.seq.nres

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return self.seq._residue_class(self.seq, ret + self.offset)
        else:
            return [self[ind] for ind in ret]


class Sequence(modobject):
    """A residue sequence"""

    #: class for residues in this sequence
    _residue_class = Residue

    #: class for list of residues in this sequence
    _residue_list_class = ResidueList

    #: class for chains in this sequence
    _chain_class = Chain

    def _indxseg(self, offset, length, indx):
        if isinstance(indx, str):
            newindx = _modeller.mod_sequence_find_chain(indx,
                                                        self.seqpt) - 1 - offset
            if newindx < 0 or (length is not None and newindx >= length):
                raise KeyError("No such chain: %s" % indx)
            return newindx
        raise TypeError("Chain IDs must be numbers or strings")

    def __get_nseg(self):
        return _modeller.mod_sequence_nseg_get(self.seqpt)
    def __get_nres(self):
        return _modeller.mod_sequence_nres_get(self.seqpt)
    def __get_range(self):
        return SeqRange(self)
    def __set_range(self, val):
        modlist.set_fixlist(self.range, val)
    def __get_source(self):
        return _modeller.mod_sequence_source_get(self.seqpt)
    def __set_source(self, val):
        _modeller.mod_sequence_source_set(self.seqpt, val)
    def __get_name(self):
        return _modeller.mod_sequence_name_get(self.seqpt)
    def __set_name(self, val):
        _modeller.mod_sequence_name_set(self.seqpt, val)
    def __get_prottyp(self):
        return _modeller.mod_sequence_prottyp_get(self.seqpt)
    def __set_prottyp(self, val):
        return _modeller.mod_sequence_prottyp_set(self.seqpt, val)
    def __get_resol(self):
        return _modeller.mod_sequence_resol_get(self.seqpt)
    def __set_resol(self, val):
        return _modeller.mod_sequence_resol_set(self.seqpt, val)
    def __get_rfactr(self):
        return _modeller.mod_sequence_rfactr_get(self.seqpt)
    def __set_rfactr(self, val):
        return _modeller.mod_sequence_rfactr_set(self.seqpt, val)
    def __get_residues(self):
        return self._residue_list_class(self)
    def __get_chains(self):
        return ChainList(self)

    nseg = property(__get_nseg, doc="Number of chains/segments")
    nres = property(__get_nres, doc="Number of residues")
    range = property(__get_range, __set_range,
                     doc="Residue number and chain ID range")
    source = property(__get_source, __set_source, doc="Source organism")
    name = property(__get_name, __set_name, doc="Protein name")
    prottyp = property(__get_prottyp, __set_prottyp,
                       doc="Protein sequence type")
    resolution = property(__get_resol, __set_resol, doc="Resolution")
    rfactor = property(__get_rfactr, __set_rfactr, doc="R factor")
    residues = property(__get_residues, doc="List of residues")
    chains = property(__get_chains, doc="List of all chains/segments")
    segments = chains
