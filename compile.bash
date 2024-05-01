#!/bin/bash


# to install ecCodes, get help here:
# https://confluence.ecmwf.int/download/attachments/73011814/eccodes_bufr_deploy.pdf?api=v2
export ECCODES_INSTALL_ROOT=$YOUR_ECCODES_INSTALLATION_DIR

mkdir bin && cd bin
make -f ../sources/ltt-modular/Makefile clean
make -f ../sources/ltt-modular/Makefile
cp *.exe ..
