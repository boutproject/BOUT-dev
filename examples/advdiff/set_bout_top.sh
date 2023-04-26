#!/usr/bin/env bash
#export BOUT_TOP=$HOME/BOUT/bout++/
export BOUT_TOP=/Users/hong/soft/BOUT/
rm -f pdb2idl.so
ln -s $BOUT_TOP/../pdb2idl/pdb2idl.so
