#!/bin/bash
TARGET_DIR=$HOME/miniconda
DOWNLOAD_TARGET=${HOME}/download
export PATH="${TARGET_DIR}/bin:$PATH"

mkdir -p ${DOWNLOAD_TARGET}

#Fetch the miniconda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${DOWNLOAD_TARGET}/miniconda.sh

#Install conda
bash ${DOWNLOAD_TARGET}/miniconda.sh -b -p ${TARGET_DIR} -f
hash -r
time conda config --set always_yes yes --set changeps1 no
time conda update -q conda
# Useful for debugging any issues with conda
conda info -a
time conda create -q -n test-environment nomkl "python>=3.5" numpy scipy netcdf4
source activate test-environment
