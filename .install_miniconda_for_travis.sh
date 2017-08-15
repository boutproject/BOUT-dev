#!/bin/bash
TARGET_DIR=$HOME/miniconda
DOWNLOAD_TARGET=${HOME}/download
export PATH="${TARGET_DIR}/bin:$PATH"

mkdir -p ${DOWNLOAD_TARGET}

#Check to see if the miniconda script already exists (thanks to caching)
#if it does then we can skip the download of conda
if [[ ! -e ${DOWNLOAD_TARGET}/miniconda.sh ]]
then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${DOWNLOAD_TARGET}/miniconda.sh
else
    echo "Using cached version of miniconda.sh"
fi

#Check to see if the miniconda build already exists (thanks to caching)
#if it does then we can skip the install of conda etc.
if [[ ! -e ${TARGET_DIR} ]]
then
    bash ${DOWNLOAD_TARGET}/miniconda.sh -b -p ${TARGET_DIR} -f
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    # Useful for debugging any issues with conda
    conda info -a
    conda create -q -n test-environment python=3.5 numpy scipy netcdf4
else
    echo "Using cached version of minconda test-environment"
fi
ls -l ${TARGET_DIR}/bin
source ${TARGET_DIR}/bin/activate test-environment
