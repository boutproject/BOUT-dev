#!/bin/sh
## Copied script from arpack-ng. Please submit fixes also to
## arpack-ng (if applicable)
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex
## Should we init a container?
if [ ".$1" = .setup ] || [ ".$1" = .podman ]
then
  # fedora
  #   note: when you PR, docker-cp provides, in the container, the branch associated with the PR (not master where there's nothing new)
  #         1. docker create --name mobydick IMAGE CMD        <=> create a container (= instance of image) but container is NOT yet started
  #         2. docker cp -a ${TRAVIS_BUILD_DIR} mobydick:/tmp <=> copy git repository (CI worker, checkout-ed on PR branch) into the container
  #                                                               note: docker-cp works only if copy from/to containers (not images)
  #         3. docker start -a mobydick                       <=> start to run the container (initialized with docker-cp)
    if test $1 = podman ; then
	cmd=podman
	# For the use of testing
	git submodule update --init # --recursive
    else
	cmd="sudo docker"
    fi
    test . != ".$2" && mpi="$2" || mpi=openmpi
    time $cmd pull ghcr.io/dschwoerer/bout-container-base:ci-fedora
    time $cmd create --cap-add=SYS_PTRACE --security-opt seccomp=unconfined \
	 --shm-size 256M \
         --name mobydick ghcr.io/dschwoerer/bout-container-base:ci-fedora \
	     /tmp/BOUT-dev/.ci_fedora.sh $mpi
    time $cmd cp ${TRAVIS_BUILD_DIR:-$(pwd)} mobydick:/tmp/BOUT-dev
    time $cmd start -a mobydick
    exit 0
fi

test . != ".$1" && mpi="$1" || mpi=openmpi

## If we are called as normal user, run test
    cp -a /tmp/BOUT-dev /home/test/
    . /etc/profile.d/modules.sh
    module load mpi/${1}-x86_64
    export OMPI_MCA_rmaps_base_oversubscribe=yes
    export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
    export TRAVIS=true
    # Try limiting openmp threads
    export FLEXIBLAS=NETLIB
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    cd
    cd BOUT-dev
    echo "starting configure"
    time cmake -S . -B build -DBOUT_USE_PETSC=ON \
	 -DBOUT_UPDATE_GIT_SUBMODULE=OFF \
	 -DBOUT_USE_SYSTEM_FMT=ON \
	 -DBOUT_USE_SYSTEM_MPARK_VARIANT=ON \
	 -DBOUT_USE_SUNDIALS=ON

    time make -C build build-check -j 2
    time make -C build check
