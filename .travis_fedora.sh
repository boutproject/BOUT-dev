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
    else
	cmd="sudo docker"
    fi
    test . != ".$2" && mpi="$2" || mpi=openmpi
    test . != ".$3" && version="$3" || version=latest
    time $cmd pull registry.fedoraproject.org/fedora:$version
    time $cmd create --name mobydick registry.fedoraproject.org/fedora:$version \
	 /tmp/BOUT-dev/.travis_fedora.sh $mpi
    time $cmd cp ${TRAVIS_BUILD_DIR} mobydick:/tmp
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
    time dnf -y install dnf-plugins-core {petsc,hdf5}-${mpi}-devel /usr/lib/rpm/redhat/redhat-hardened-cc1
    # Allow to override packages - see #2073
    time dnf copr enable -y davidsch/fixes4bout || :
    time dnf -y upgrade
    time dnf -y builddep bout++
    useradd test
    cp -a /tmp/BOUT-dev /home/test/
    chown -R test /home/test
    chmod u+rwX /home/test -R
    sudo -u test ${0/\/tmp/\/home\/test} $mpi
## If we are called as normal user, run test
else
    . /etc/profile.d/modules.sh
    module load mpi/${1}-x86_64
    export OMPI_MCA_rmaps_base_oversubscribe=yes
    export TRAVIS=true
    export FLEXIBLAS=NETLIB
    cd
    cd BOUT-dev
    echo "starting configure"
    time ./configure --with-petsc --enable-shared || cat config.log
    sed -e "s|-L/usr/lib64 ||g" -i make.config
    for f in tests/requirements/*[^y] ; do
	echo -n "$f: "
	$f && echo yes || echo no
    done
    time make build-check -j 2
    time make check
fi
