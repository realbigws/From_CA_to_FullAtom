import modutil

__docformat__ = "epytext en"

class FixList(object):
    """A fixed-length list of Modeller data"""
    def __len__(self):
        raise NotImplementedError

    def _getfunc(self, indx):
        """Get the data at C{indx}"""
        raise NotImplementedError

    def _setfunc(self, indx, val):
        """Set the data at C{indx} to C{val}"""
        raise TypeError("this list is read-only")

    def __str__(self):
        return str([val for val in self])
    def __repr__(self):
        return str(self)

    def _lookup_index(self, indx, require_inrange=True):
        return modutil.handle_seq_indx(self, indx,
                                       require_inrange=require_inrange)

    def __getitem__(self, indx):
        ret = self._lookup_index(indx)
        if isinstance(ret, int):
            return self._getfunc(ret)
        else:
            return [self[i] for i in ret]

    def __setitem__(self, indx, val):
        ret = self._lookup_index(indx)
        if isinstance(ret, int):
            self._setfunc(ret, val)
        else:
            if not isinstance(val, (list, tuple)) or len(val) != len(ret):
                raise ValueError("need a list of size %d for assignment"
                                 % len(ret))
            for (i, v) in zip(ret, val):
                self[i] = v


class AppendList(FixList):
    """A list of Modeller data which can be added to but not reordered"""

    def _setfunc(self, indx, val):
        raise TypeError("this list can only be appended to")

    def append(self, obj):
        """Append obj to the end of the list"""
        raise NotImplementedError

    def extend(self, objs):
        """Extend list by adding elements from objs"""
        for obj in objs:
            self.append(obj)


class VarList(FixList):
    """A variable-length list of Modeller data"""

    def _setdimfunc(self, newdim):
        raise NotImplementedError

    def clear(self):
        """Clear all items from the list"""
        self._setdimfunc(0)

    def append(self, obj):
        """Append obj to the end of the list"""
        dim = len(self)
        self._setdimfunc(dim + 1)
        self[dim] = obj

    def extend(self, objs):
        """Extend list by adding elements from objs"""
        dim = len(self)
        self._setdimfunc(dim + len(objs))
        for (n, obj) in enumerate(objs):
            self[dim + n] = obj

    def pop(self, index=-1):
        """Remove and return item at index (default last)"""
        if len(self) == 0:
            raise IndexError("pop from empty list")
        x = self[index]
        del self[index]
        return x

    def insert(self, indx, obj):
        """Insert object before index"""
        indx = self._lookup_index(indx, require_inrange=False)
        if not isinstance(indx, int):
            raise TypeError("an integer is required")
        dim = len(self)
        self._setdimfunc(dim+1)
        for n in range(dim, indx, -1):
            self[n] = self[n-1]
        self[indx] = obj

    def __delitem__(self, indx):
        ret = self._lookup_index(indx)
        if isinstance(ret, int):
            dim = len(self)
            for n in range(ret+1, dim):
                self[n-1] = self[n]
            self._setdimfunc(dim-1)
        else:
            ret.sort(reverse=True)
            # Not very efficient, particularly for long lists!
            for i in ret:
                del self[i]


def set_fixlist(lst, obj):
    """Helper function to set members that are FixList objects"""
    if len(obj) != len(lst):
        raise ValueError("expecting a tuple or list of length %d" % len(lst))
    for (n, val) in enumerate(obj):
        lst[n] = val

def set_varlist(lst, obj):
    """Helper function to set members that are VarList objects"""
    lst.clear()
    lst.extend(obj)

def del_varlist(lst):
    """Helper function to delete members that are VarList objects"""
    lst.clear()


class LinkList(VarList):
    """A variable-length linked list of Modeller data"""

    def _delfunc(self, indx):
        raise NotImplementedError

    def _insfunc(self, indx, val):
        raise NotImplementedError

    def _setfunc(self, indx, val):
        self._delfunc(indx)
        self._insfunc(indx, val)

    def append(self, obj):
        """Append obj to the end of the list"""
        self.insert(len(self), obj)

    def extend(self, objs):
        """Extend list by adding elements from objs"""
        for obj in objs:
            self.append(obj)

    def insert(self, indx, obj):
        """Insert object before index"""
        indx = self._lookup_index(indx, require_inrange=False)
        if not isinstance(indx, int):
            raise TypeError("an integer is required")
        self._insfunc(indx, obj)

    def __delitem__(self, indx):
        ret = self._lookup_index(indx)
        if isinstance(ret, int):
            self._delfunc(ret)
        else:
            ret.sort(reverse=True)
            for i in ret:
                self._delfunc(i)


class SimpleVarList(VarList):
    def __init__(self, modpt, getdimfunc, setdimfunc, getfunc, setfunc):
        VarList.__init__(self)
        self.__modpt = modpt
        self.__getdimfunc = getdimfunc
        self.__setdimfunc = setdimfunc
        self.__setfunc = setfunc
        self.__getfunc = getfunc
    def __len__(self):
        return self.__getdimfunc(self.__modpt)
    def _getfunc(self, indx):
        return self.__getfunc(self.__modpt, indx)
    def _setfunc(self, indx, val):
        self.__setfunc(self.__modpt, indx, val)
    def _setdimfunc(self, num):
        self.__setdimfunc(self.__modpt, num)
