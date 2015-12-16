#!/usr/bin/env bash
unset PYTHONPATH

echo "setting up science tools"
#VERSION=09-34-03
VERSION=10-00-03 # this one handles gaps in the FT2!

# BLDARCH is defined for convenience
if [ $(fs sysname | grep -c amd64_rhel60) -gt 0 ]; then
    # 64-bit redhat6 at SLAC.
    export BLDARCH=redhat6-x86_64-64bit-gcc44
elif [ $(fs sysname | grep -c amd64_rhel50) -gt 0 ]; then
    # 64-bit redhat5 at SLAC.
    export BLDARCH=redhat5-x86_64-64bit-gcc41
else
    # 32-bit redhat5 at SLAC.
    export BLDARCH=redhat5-i686-32bit-gcc41
fi

export GLAST_EXT=/afs/slac/g/glast/ground/GLAST_EXT/${BLDARCH}

export PFILES=".;" # Clear PFILES

export BUILDS=/nfs/farm/g/glast/u35/ReleaseManagerBuild
export BLDTYPE=Optimized

if [ "$2" = "HEAD" ]; then
    export INST_DIR=/u/gl/mdwood/ki10/ScienceTools-scons
    source ${INST_DIR}/bin/${BLDARCH}-Debug-Optimized/_setup.sh
elif [ -n "$2" ]; then
    export INST_DIR=${BUILDS}/${BLDARCH}/${BLDTYPE}/ScienceTools/$2
    source ${INST_DIR}/bin/${BLDARCH}-${BLDTYPE}/_setup.sh
else
    export INST_DIR=${BUILDS}/${BLDARCH}/${BLDTYPE}/ScienceTools/$VERSION
    source ${INST_DIR}/bin/${BLDARCH}-${BLDTYPE}/_setup.sh
fi

echo "add pointlike pointlike HEAD"
export PYTHONPATH=${INST_DIR}/pointlike/python:$PYTHONPATH

echo "REMEMBER TO SET DATADIR!"