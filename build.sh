#!/bin/sh

set -x
g++-3.0 -static -DNDEBUG -O3 -g -gstabs -W -Wall -I/nfs_exported/mpich/include -L/nfs_exported/mpich/lib -fverbose-asm -save-temps -o tdg tdg.cxx -lmpich
