#!/bin/bash

make

mpirun ./elm_pb -d test


echo "Finished ELM-pb test case"

