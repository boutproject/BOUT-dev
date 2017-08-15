#!/bin/bash
TARGET_DIR=$HOME/miniconda

export PATH="${TARGET_DIR}/bin:$PATH"

#Check to see if the target directory already exists (thanks to caching)
#if it does then we can skip the install of conda etc.
# if [ ! -e ${TARGET_DIR} ]
# then
#     wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# else
#     echo "Using cached version of conda"
# fi

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ${TARGET_DIR}
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Useful for debugging any issues with conda
conda info -a
conda create -q -n test-environment python=3.5 numpy scipy netcdf4
source activate test-environment
