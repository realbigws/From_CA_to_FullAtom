import sys
try:
    import subprocess
except ImportError:
    from modeller.python_library import subprocess

__docformat__ = "epytext en"

def myspawn(cmd, output):
    """Run C{cmd} in the background, and direct stdout and stderr to
       C{output}."""

    fp = file(output, "w")
    print "%s >& %s" % (cmd, output)
    if sys.platform == 'win32':
        _myspawn_win32(cmd, fp)
    else:
        _myspawn_unix(cmd, fp)

def _myspawn_unix(cmd, fp):
    p = subprocess.Popen(cmd, shell=True, stdout=fp, stderr=subprocess.STDOUT)

def _myspawn_win32(cmd, fp):
    try:
        # shell isn't needed on Win32, and may not be found under wine anyway
        p = subprocess.Popen(cmd, shell=False, stdout=fp,
                             stderr=subprocess.STDOUT)
    # Ignore Windows "file not found" errors, so that behavior is consistent
    # between Unix and Windows
    except WindowsError, detail:
        print "WindowsError: %s (ignored)" % detail
