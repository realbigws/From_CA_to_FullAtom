"""Utility classes to provide a portable way to spawn a child process and
   communicate with its stdin and combined stdout/stderr."""

__all__ = ["MyPopen4"]

__docformat__ = "epytext en"

import sys
import os
try:
    import subprocess
except ImportError:
    from modeller.python_library import subprocess

def _remove_library_paths(env):
    """Make a new environment, with library path variables removed. This is
       to prevent architecture-dependent variables from being passed to
       processes which may run on a different architecture machine."""
    env = env.copy()
    for key in ('LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'LIBPATH',
                'PYTHONPATH'):
        if key in env:
            del env[key]
    # Restore user-set PYTHONPATH if present:
    if 'ORIGPYPATH' in env:
        env['PYTHONPATH'] = env.pop('ORIGPYPATH')
    return env


class MyPopen4(subprocess.Popen):
    """Utility class to provide a portable way to spawn a child process and
       communicate with its stdin and combined stdout/stderr."""

    def __init__(self, cmd):
        # shell isn't needed on Win32, and may not be found under wine anyway
        shell = (sys.platform != "win32")
        env = _remove_library_paths(os.environ)
        subprocess.Popen.__init__(self, cmd, shell=shell, stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, env=env)

    def require_clean_exit(self):
        """Make sure the child exited with a zero return code"""
        r = self.wait()
        if r != 0:
            raise IOError("Process failed with exit status %d" % r)
