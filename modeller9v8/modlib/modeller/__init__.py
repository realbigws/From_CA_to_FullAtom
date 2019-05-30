"""MODELLER, a package for protein structure modeling

   See U{http://salilab.org/modeller/} for further details.

   @author: Andrej Sali
   @copyright: 1989-2010 Andrej Sali
"""

__all__ = ['energy_data', 'io_data', 'environ', 'group_restraints',
           'ModellerError', 'FileFormatError', 'StatisticsError',
           'SequenceMismatchError', 'model', 'alignment',
           'sequence_db', 'profile', 'saxsdata', 'density', 'pssmdb',
           'excluded_pair', 'rigid_body', 'symmetry', 'selection', 'log',
           'info', 'pseudo_atom', 'virtual_atom', 'modfile', 'features',
           'forms', 'secondary_structure', 'terms', 'physical']

__docformat__ = "epytext en"

# Version check
import sys
if sys.version_info[0] < 2 \
   or (sys.version_info[0] == 2 and sys.version_info[1] < 3):
    raise ImportError("This module requires Python 2.3 or later")

try:
    import config
except ImportError:
    config = None

def get_python_api_ver():
    """Get the Python version at which the C API last changed."""
    ver = sys.version_info[0:2]
    if ver == (2,6):
        return (2,5)
    else:
        return ver

# Special processing on Windows to find _modeller.pyd and Modeller DLLs:
if hasattr(config, 'install_dir') and hasattr(sys, 'dllhandle'):
    dpath = config.install_dir + '\\modlib\\python%d.%d' % sys.version_info[:2]
    if dpath not in sys.path:
        # Insert *after* first entry, so as not to break parallel module's
        # propagation of the first entry (that containing the running script's
        # directory) to slaves
        sys.path.insert(1, dpath)
    try:
        import os
        dpath = config.install_dir + '\\bin;'
        if dpath not in os.environ['PATH']:
            os.environ['PATH'] = dpath + os.environ['PATH']
        del os
    except ImportError:
        pass
    del dpath
# Add Python version-specific directory to search path:
elif hasattr(config, 'install_dir'):
    try:
        import os.path, re, sys
        srch = re.compile("%s/*lib/[^/]+/?" % config.install_dir)
        for (n, pathcomp) in enumerate(sys.path):
            if srch.match(pathcomp):
                modpath = os.path.join(pathcomp,
                                       'python%d.%d' % get_python_api_ver())
                if modpath not in sys.path:
                    sys.path.insert(n, modpath)
                break
        del re, n, pathcomp, os, srch
    except ImportError:
        pass

# Set Modeller install location and license
import _modeller
if hasattr(config, 'license'):
    _modeller.mod_license_key_set(config.license)
if hasattr(config, 'install_dir'):
    _modeller.mod_install_dir_set(config.install_dir)

_modeller.mod_start()

from energy_data import energy_data
from io_data import io_data
from environ import environ
from group_restraints import group_restraints
from error import *
from model import model
from alignment import alignment
from sequence_db import sequence_db
from profile import profile
from saxsdata import saxsdata
from density import density
from pssmdb import pssmdb
from excluded_pair import excluded_pair
from rigid_body import rigid_body
from symmetry import symmetry
from selection import selection
from util.logger import log
from information import info
import pseudo_atom
import virtual_atom
import modfile
import features
import forms
import secondary_structure
import terms
import physical

# Load in readline, if available, to make interactive use easier
try:
    if len(sys.argv) > 0 and sys.argv[0] == '-' and sys.stdin.isatty():
        import readline
except (ImportError, AttributeError):
    pass

# Set job name
if len(sys.argv) > 0 and sys.argv[0] != '-':
    nam = sys.argv[0]
    if nam.endswith('.py'):
        nam = nam[:-3]
    info.jobname = nam
    del nam

del sys, _modeller, config
