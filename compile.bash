#!/bin/bash

flags="-fcheck=all -fbacktrace -fmax-errors=10"
IXX="-I${ECCODES_INSTALL_ROOT}/include"
LXX="-L${ECCODES_INSTALL_ROOT}/lib"
source_folder="$1"
exe="$2"

test -d binaries || mkdir -p binaries/
cd binaries/

gfortran -g $flags -o $exe $source_folder/*.f90 $IXX $LXX -Wl,-rpath=${ECCODES_INSTALL_ROOT}/lib -leccodes_f90
cp $exe ../$exe
