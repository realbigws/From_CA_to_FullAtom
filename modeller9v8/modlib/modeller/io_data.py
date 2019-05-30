import _modeller
from modeller.util.modobject import modobject
from modeller.util.logger import log

class io_data(modobject):
    """Controls reading from/writing to atom files"""

    hetatm = False
    hydrogen = False
    water = False
    __atom_files_directory = ['']
    __modpt = None

    def __init__(self, copy=None, **kwargs):
        self.__modpt = _modeller.mod_io_data_new()
        if copy:
            for member in copy.__dict__:
                if 'io_data' not in member:
                    self.__dict__[member] = copy.__dict__[member]
            # Copy the list, rather than making another reference
            self.atom_files_directory = copy.atom_files_directory[:]
        for key in kwargs:
            if key in dir(io_data):
                exec("self."+key+"="+str(kwargs[key]))
            else:
                raise KeyError(str(key))

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_io_data_new()

    def __del__(self):
        if self.__modpt:
            _modeller.mod_io_data_free(self.__modpt)

    def __get_atfil(self):
        return self.__atom_files_directory
    def __set_atfil(self, val):
        if isinstance(val, tuple):
            self.__atom_files_directory = list(val)
        elif isinstance(val, list):
            self.__atom_files_directory = val
        elif isinstance(val, str):
            spl = val.split(':')
            log.warning("io_data",
"""Setting io.atom_files_directory to a colon-delimited string is
              deprecated, as it is not robust on Windows systems. Set it to
              a list of directories instead. For example:
              env.io.atom_files_directory = """ + str(spl))
            self.__atom_files_directory = spl
        else:
            raise TypeError("atom_files_directory should be a list")
    def __get_modpt(self):
        modpt = self.__modpt
        _modeller.mod_io_data_set(modpt, self.hydrogen, self.hetatm, self.water,
                                  self.__atom_files_directory)
        return modpt

    atom_files_directory = property(__get_atfil, __set_atfil,
                                    doc="List of paths to search " + \
                                        "for atom files")
    modpt = property(__get_modpt)
