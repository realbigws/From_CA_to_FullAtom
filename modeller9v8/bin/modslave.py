#!/usr/bin/python

"""A simple script to start up a Modeller slave process. This is useful if,
   for example, you want to use your system's native Python version rather
   than the Python 2.3 interpreter compiled into Modeller.

   Use by setting the 'modeller_path' variable to the name of this script when
   creating a 'job' object. (Note that if your Modeller Python modules are not
   in the default search path, you may have to write a wrapper script to set
   PYTHONPATH and/or LD_LIBRARY_PATH or similar so that the modules can be
   located.)
"""

import sys

def usage():
    print "Usage: %s -slave masterspec" % sys.argv[0]
    print
    sys.exit(1)

if len(sys.argv) != 3 or sys.argv[1] != '-slave':
    usage()
else:
    from modeller.parallel import *
    # Ensure that the current directory is in the search path:
    sys.path.insert(0, '')
    slaveloop(sys.argv[2])
