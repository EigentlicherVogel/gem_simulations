#!/bin/bash

export GARFIELD_HOME=/garfieldpp
export GARFIELD_INSTALL=$GARFIELD_HOME/install
export HEED_DATABASE=$GARFIELD_HOME/Heed/heed++/database
export CMAKE_PREFIX_PATH=$GARFIELD_HOME/install:${CMAKE_PREFIX_PATH}

source /usr/local/root/bin/thisroot.sh

#export LD_LIBRARY_PATH=/garfieldpp/install/lib:/usr/lib:/lib:/usr/local/root/lib
export LD_LIBRARY_PATH=/garfieldpp/install/lib:/usr/local/root/lib

echo $ROOTSYS
pwd
#root

/home/zoua1/gem_simulations/sim_newver/build/rdo_iter $1
