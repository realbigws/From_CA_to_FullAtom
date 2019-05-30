import _modeller
from modeller.util.modobject import modobject

class group_restraints(modobject):
    """Holds restraints which act on atom classes/groups"""
    _modpt = None
    env = None
    __class_file = None
    __param_files = None

    def __init__(self, env, classes, parameters=None):
        self.env = env.copy()
        self._modpt = _modeller.mod_group_restraints_new(self)
        self.__read_classes(classes)
        self.__param_files = []
        if parameters:
            self.append(parameters)

    def __del__(self):
        if self._modpt:
            _modeller.mod_group_restraints_free(self._modpt)

    def __setstate__(self, d):
        """Restore internal information from files"""
        self.__dict__.update(d)
        self._modpt = _modeller.mod_group_restraints_new(self)
        self.__read_classes(self.__class_file)
        params = self.__param_files
        self.__param_files = []
        for file in params:
            self.append(file)

    def __read_classes(self, file):
        """Read atom classes from a file"""
        self.__class_file = file
        return _modeller.mod_atom_classes_read(self._modpt, file)

    def append(self, file):
        """Read interaction parameters from a file"""
        self.__param_files.append(file)
        return _modeller.mod_group_restraints_read(self._modpt,
                                                   self.env.libs.modpt, file)
