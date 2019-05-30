"""Classes for file handling."""

__docformat__ = "epytext en"

import _modeller

def delete(file):
    """Delete a file"""
    return _modeller.mod_file_delete(file)


def inquire(file):
    """Check if file exists"""
    return _modeller.mod_file_exists(file)


def default(root_name='undf', file_id='X', id1=1, id2=1, file_ext=''):
    """Generate a default Modeller-style file name"""
    return "%s%s%04d%04d%s" % (root_name, file_id, id1, id2, file_ext)

class File(object):
    """Simple wrapper around Modeller's internal file-reading routines"""

    def __init__(self, filename, mode='r'):
        self._file = _modeller.mod_file_open(filename, mode)

    def __del__(self):
        if hasattr(self, "_file"):
            _modeller.mod_file_close(self._file)

    def close(self):
        """Close the file handle"""
        _modeller.mod_file_close(self.file_pointer)
        del self._file

    def write(self, str):
        """Write a string to the file"""
        _modeller.mod_file_write_buffer(self.file_pointer, str)

    def read(self, len):
        """Read a fixed-length string from the file"""
        return _modeller.mod_file_read_buffer(self.file_pointer, len)[1]

    def __get_file_pointer(self):
        if hasattr(self, "_file"):
            return self._file
        else:
            raise ValueError, "I/O operation on closed file"
    file_pointer = property(__get_file_pointer)


class _CStringIOWrapper(File):
    """Simple wrapper to allow Modeller code to use Python
       cStringIO objects as mod_file objects."""
    def __init__(self, str):
        self._file = _modeller.mod_file_open_cstring(str)


def _check_filehandle(fh):
    """Convert a Python cStringIO object into a Modeller filehandle
       if necessary"""
    if isinstance(fh, File):
        return fh
    else:
        return _CStringIOWrapper(fh)
