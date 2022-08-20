#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
-D CMAKE_INSTALL_PREFIX=/mnt/c/share/calc_wss \
..

make && make install