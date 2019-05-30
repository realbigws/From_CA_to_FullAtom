"""Classes for handling template structures in alignments"""

import _modeller
from modeller import sequence, alnsequence, coordinates, modfile

__docformat__ = "epytext en"

class Residue(coordinates.Residue, alnsequence.Residue):
    """A single residue in a template structure."""

    def __get_curv(self):
        x = _modeller.mod_structure_curvn_get(self.mdl._strucpt)
        return _modeller.mod_float1_get(x, self._num)
    def __set_curv(self, val):
        x = _modeller.mod_structure_curvn_get(self.mdl._strucpt)
        _modeller.mod_float1_set(x, self._num, val)
    curvature = property(__get_curv, __set_curv, doc="Mainchain curvature")


class Structure(coordinates.Coordinates, alnsequence.Sequence):
    """A single template structure within an alignment"""

    __num = None
    __read_coord = False
    _residue_class = Residue

    def __init__(self, aln, num):
        self.__num = num
        alnsequence.Sequence.__init__(self, aln, num)
        coordinates.Coordinates.__init__(self)

    def __repr__(self):
        return "Structure %s" % repr(self.code)
    def __str__(self):
        return "<%s>" % repr(self)

    def write(self, file):
        """Write current coordinates to a PDB file"""
        if isinstance(file, str):
            file = modfile.File(file, 'w')
        return _modeller.mod_alnstructure_write(self.cdpt, self.seqpt,
                                                file.file_pointer,
                                                self.aln.env.libs.modpt)

    def reread(self):
        """Reread coordinates from the atom file. Useful for restoring the
           original template orientation."""
        self.__read_coord = False
        _modeller.mod_alnstructure_invalidate(self.aln.modpt, self.__num)
        self.__get_cdpt()   # Force reread of structural info

    def __get_strucpt(self):
        if not self.__read_coord:
            aln = self.aln
            _modeller.mod_alnstructure_read(aln.modpt, self.__num,
                                            aln.env.io.modpt,
                                            aln.env.libs.modpt)
            self.__read_coord = True
        return _modeller.mod_alignment_structure_get(self.aln.modpt, self.__num)

    def __get_cdpt(self):
        return _modeller.mod_structure_cd_get(self.__get_strucpt())

    _strucpt = property(__get_strucpt)
    cdpt = property(__get_cdpt)
