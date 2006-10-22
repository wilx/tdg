#!/bin/sh

set -x
g++-3.0 -DNDEBUG -O3 -static -g3 -gstabs -W -Wall -I/nfs_exported/mpich/include -L/nfs_exported/mpich/lib -L/nfs_exported/gm/lib -fverbose-asm -save-temps -o tdg-gm tdg.cxx -lmpich -lgm
