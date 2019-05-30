The files in this directory are used to build the Modeller Python interface
from the C header files (in the ../include/ directory).

To rebuild the interface:

1. Use SWIG (http://www.swig.org/ - you will need version 1.3.31 or later) to
   build a wrapper file from the .i files, using a command similar to the
   following:
   % swig -python -keyword -nodefaultctor -nodefaultdtor -noproxy modeller.i

2. Edit the file setup.py in this directory to put in your Modeller executable
   type - e.g. i386-intel8 on 32-bit Linux boxes.

3. Build the Python extension with a command similar to the following
   (you will need to have a full Python installation, including 'development'
   files):
   % python setup.py build
