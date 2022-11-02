#!/usr/bin/env bash

# install dependencies
apt-get update
apt-get install -y mpich2 libmpich2-dev 
apt-get install -y libfftw3-dev libnetcdf-dev
apt-get install -y g++ make
apt-get install -y python-scipy

# get the bout++ code
apt-get install -y git
rm -rf BOUT
git clone https://github.com/bendudson/BOUT.git

# environment variables for both this session and the vagrant user
echo "export IDL_PATH=$(pwd)/BOUT/tools/idllib" >> $(pwd)/.bashrc
echo "export PYTHONPATH=$(pwd)/BOUT/tools/pylib/:$PYTHONPATH" >> $(pwd)/.bashrc
source $(pwd)/.bashrc
export IDL_PATH=$(pwd)/BOUT/tools/idllib
export PYTHONPATH=$(pwd)/BOUT/tools/pylib/:$PYTHONPATH

# configure bout++
cd BOUT
./configure
make
cd examples
./test_suite

# set owner to vagrant user so compilation works
cd ../..
chown -R vagrant:vagrant BOUT
