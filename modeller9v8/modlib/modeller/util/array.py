"""Basic classes for accessing Modeller pointer arrays"""

import _modeller
from modeller.util.modlist import FixList

class BaseArray(object):
    """Virtual base class for all pointer arrays"""

    def __init__(self, modpt, readonly=False):
        self.__modpt = modpt
        self.readonly = readonly

    def __check_indx(self, indx):
        if self._rank == 1:
            if not isinstance(indx, int):
                raise TypeError("Index should be an integer")
        else:
            if not isinstance(indx, (list, tuple)) or len(indx) != self._rank:
                raise TypeError("Index should be a %d-element tuple"
                                % self._rank)

    def __getitem__(self, indx):
        self.__check_indx(indx)
        return self._getfunc(self.__modpt, indx)

    def __setitem__(self, indx, val):
        if self.readonly:
            raise ValueError("array is read-only")
        self.__check_indx(indx)
        return self._setfunc(self.__modpt, indx, val)

class Float3DArray(BaseArray):
    """3D array of floating-point values"""
    _rank = 3
    def _getfunc(self, modpt, indx):
        return _modeller.mod_float3_get(modpt, *indx)
    def _setfunc(self, modpt, indx, val):
        _modeller.mod_float3_set(modpt, indx[0], indx[1], indx[2], val)

class Double1DArray(FixList):
    """1D array of double-precision floating-point values"""
    def __init__(self, modpt, dimfunc):
        self.modpt = modpt
        self.dimfunc = dimfunc
    def _getfunc(self, indx):
        return _modeller.mod_double1_get(self.modpt, indx)
    def _setfunc(self, indx, val):
        _modeller.mod_double1_set(self.modpt, indx, val)
    def __len__(self):
        return self.dimfunc()

class Double2DArray(BaseArray):
    """2D array of double-precision floating-point values"""
    _rank = 2
    def _getfunc(self, modpt, indx):
        return _modeller.mod_double2_get(modpt, *indx)
    def _setfunc(self, modpt, indx, val):
        _modeller.mod_double2_set(modpt, indx[0], indx[1], val)
