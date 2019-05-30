import _modeller
from modeller.util.modobject import modobject

class pssmdb(modobject):
    """Holds a database of PSSMs"""
    __modpt = None
    env = None

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_pssmdbobj_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_pssmdbobj_new(self)

    def __del__(self):
        if self.__modpt:
            _modeller.mod_pssmdbobj_free(self.__modpt)

    def __get_modpt(self):
        return self.__modpt

    def read(self, pssmdb_name, pssmdb_format):
        return _modeller.mod_pssmdbobj_read(self.__modpt, pssmdb_name,
                                            pssmdb_format)

    modpt = property(__get_modpt)
