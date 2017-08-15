#!/bin/bash
TARGET_DIR=$HOME/miniconda

export PATH="${TARGET_DIR}/bin:$PATH"

mkdir -p download

#Check to see if the target directory already exists (thanks to caching)
#if it does then we can skip the install of conda etc.
if [[ ! -e download/miniconda.sh ]]
then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O download/miniconda.sh
else
    echo "Using cached version of miniconda"
fi

bash download/miniconda.sh -b -p ${TARGET_DIR} -f
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Useful for debugging any issues with conda
conda info -a
conda create -q -n test-environment python=3.5 numpy scipy netcdf4
ls
echo "------------"
ls ${TARGET_DIR}/bin
source activate test-environment
