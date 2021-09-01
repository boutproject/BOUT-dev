#!/usr/bin/env bash
module purge
export SPACK_ROOT=/usr/WS2/BOUT-GPU/lassen/env/spack
. $SPACK_ROOT/share/spack/setup-env.sh
module refresh
module --ignore-cache load gcc/8.3.1
module --ignore-cache load cmake/3.14.5
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cuda/10.1.243
module --ignore-cache load totalview/2020.1.13
echo "Loading BOUT-GPU Group Environment"
source /usr/WS2/BOUT-GPU/lassen/env/spack/var/spack/environments/bout/loads

module list


