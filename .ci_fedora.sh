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
    test . != ".$3" && version="$3" || version=rawhide
    time $cmd pull registry.fedoraproject.org/fedora:$version
    time $cmd create --cap-add=SYS_PTRACE --security-opt seccomp=unconfined \
	 --shm-size 256M \
         --name mobydick registry.fedoraproject.org/fedora:$version \
	     /tmp/BOUT-dev/.ci_fedora.sh $mpi
    time $cmd cp ${TRAVIS_BUILD_DIR:-$(pwd)} mobydick:/tmp/BOUT-dev
    time $cmd start -a mobydick
    exit 0
fi

test . != ".$1" && mpi="$1" || mpi=openmpi

## If we are called as root, setup everything
if [ $UID -eq 0 ]
then
    cat /etc/os-release
    # Ignore weak depencies
    echo "install_weak_deps=False" >> /etc/dnf/dnf.conf
    echo "minrate=10M" >> /etc/dnf/dnf.conf
    export FORCE_COLUMNS=200
    time dnf -y install dnf5
    time dnf5 -y install dnf5-plugins cmake python3-zoidberg python3-natsort
    # Allow to override packages - see #2073
    time dnf5 copr enable -y davidsch/fixes4bout || :
    time dnf5 -y upgrade
    time dnf5 -y builddep bout++
    useradd test
    cp -a /tmp/BOUT-dev /home/test/
    chown -R test /home/test
    chmod u+rwX /home/test -R
    su - test -c "${0/\/tmp/\/home\/test} $mpi"
## If we are called as normal user, run test
else
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
fi
