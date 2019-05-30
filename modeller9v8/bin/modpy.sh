#!/bin/sh

# Install location
#MODINSTALL9v8="/home/RaptorX/bin/modeller9v8"
#MODINSTALL9v8="/mnt/data/RaptorXCommon/RaptorX-Threading/BuildModel_Package/modeller9v8"
if [ -z "${MODINSTALL9v8}" ]; then
	echo "Please set the environmental variable MODINSTALL9v8 to the installation folder of the modeller9v8 module"
	exit 1
fi
#TOPDIR="/home/RaptorX/bin/modeller9v8"
TOPDIR=$MODINSTALL9v8
EXETYPE=x86_64-intel8

LLP=${TOPDIR}/lib/${EXETYPE}
if test -z "${LD_LIBRARY_PATH}"; then
  LD_LIBRARY_PATH=${LLP}
else
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LLP}
fi
if test -z "${DYLD_LIBRARY_PATH}"; then
  DYLD_LIBRARY_PATH=${LLP}
else
  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${LLP}
fi
if test -z "${LIBPATH}"; then
  LIBPATH=${LLP}
else
  LIBPATH=${LIBPATH}:${LLP}
fi

if test "x${EXETYPE}" = "xi386-w32"; then
  PP=${TOPDIR}/modlib
else
  PP=${TOPDIR}/lib/${EXETYPE}:${TOPDIR}/modlib
fi
if test -z "${PYTHONPATH}"; then
  PYTHONPATH=${PP}
else
  ORIGPYPATH="${PYTHONPATH}"
  if test "x${EXETYPE}" = "xi386-w32"; then
    PYTHONPATH="${PYTHONPATH};${PP}"
  else
    PYTHONPATH=${PYTHONPATH}:${PP}
  fi
fi
export LD_LIBRARY_PATH DYLD_LIBRARY_PATH PYTHONPATH LIBPATH ORIGPYPATH
exec "$@"
