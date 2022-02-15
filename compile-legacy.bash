#!/bin/bash

flags="-fcheck=all -fbacktrace -fmax-errors=10"
IXX="-I${ECCODES_INSTALL_ROOT}/include"
LXX="-L${ECCODES_INSTALL_ROOT}/lib"
file="$1"
exe="$2"

gfortran -g $flags -o $exe $file $IXX $LXX -Wl,-rpath=${ECCODES_INSTALL_ROOT}/lib -leccodes_f90
