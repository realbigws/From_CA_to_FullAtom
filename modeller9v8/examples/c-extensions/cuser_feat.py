from modeller import *
from modeller.scripts import complete_pdb

# Load the C extension module; this needs to be compiled with commands like
# the following (these work for most Linux systems, where 'modXXX' is your
# Modeller binary):
#
# swig -python -noproxy cuser_feat.i
# gcc -shared -Wall -fPIC `modXXX --cflags --libs` \
#     `pkg-config --cflags glib-2.0` \
#     -I/usr/include/python2.4 cuser_feat.c cuser_feat_wrap.c \
#     -o _cuser_feat.so -lm
#
# For a typical Mac system, something like the following should work:
#
# swig -python -noproxy cuser_feat.i
# gcc -bundle -flat_namespace -undefined suppress -Wall `modXXX --cflags` \
#     `pkg-config --cflags --libs glib-2.0` cuser_feat.c cuser_feat_wrap.c \
#     -o _cuser_feat.so -I /System/Library/Frameworks/Python.framework/Headers/
#
# On AIX, the following should work:
#
# swig -python -noproxy cuser_feat.i
# cc -qmkshrobj `modXXX --cflags --libs` `pkg-config --cflags --libs glib-2.0` \
#    -I/usr/local/include/python2.3 cuser_feat.c cuser_feat_wrap.c \
#    -o _cuser_feat.so -lm -bI:/usr/local/lib/python2.3/config/python.exp
import _cuser_feat

env = environ()

env.io.atom_files_directory = ['../atom_files']
log.verbose()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

class MyDist(features.feature):
    """An implementation of Modeller's distance feature (type 1) as a
       C extension module"""

    numatoms = 2
    _builtin_index = _cuser_feat.myfeat_create()


mdl = complete_pdb(env, "1fdn")
sel = selection(mdl)
rsr = mdl.restraints
at = mdl.atoms
rsr.add(forms.gaussian(group=physical.bond,
                       feature=MyDist(at['CA:1'], at['C:1']),
                       mean=1.5380, stdev=0.0364))
sel.energy()
