import _modeller

class Information(object):
    """Allows access to information about the current Modeller build"""

    def time_mark(self):
        """Returns the current date, time, and CPU time used"""
        return _modeller.mod_time_mark()

    def __get_version(self):
        return _modeller.mod_long_version_get()
    def __get_version_info(self):
        mod_version = _modeller.mod_short_version_get()
        try:
            return (int(mod_version[0]), int(mod_version[2]))
        except ValueError:
            return mod_version
    def __get_build_date(self):
        return _modeller.mod_build_date_get()
    def __get_exe_type(self):
        return _modeller.mod_exe_type_get()
    def __get_debug(self):
        return _modeller.mod_debug_get()
    def __get_bindir(self):
        return _modeller.mod_bindir_get()
    def __get_jobname(self):
        return _modeller.mod_jobname_get()
    def __set_jobname(self, val):
        return _modeller.mod_jobname_set(val)

    version = property(__get_version, doc="The Modeller version, as a string")
    version_info = property(__get_version_info,
                            doc="The Modeller major, minor version numbers")
    build_date = property(__get_build_date,
                          doc="The date this Modeller binary was built")
    exe_type = property(__get_exe_type,
                        doc="The executable type of this Modeller binary")
    debug = property(__get_debug, doc="Whether this is a debugging build")
    bindir = property(__get_bindir,
                      doc="Directory containing Modeller binaries")
    jobname = property(__get_jobname, __set_jobname,
                       doc="Name of the current Modeller job")


# Global information object
info = Information()
