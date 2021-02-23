#!/usr/bin/env bash

cmake_config=${CMAKE_CONFIG:-""}
if [[ -z ${cmake_config} ]]
then
    arch=$(uname -m)
    compiler=gcc
    cmake_config="tests/gitlab/${arch}-${compiler}.cmake"
else
    arch_comp=${cmake_config##*/}
    arch_comp=${arch_comp%.cmake}
    compiler=${arch_comp#*-}
    arch=${arch_comp%-*}
    echo "$arch and $compiler"
fi

env_prefix=/g/g91/bernede1/Projects/BOUT-dev

build_base_dir=${env_prefix}
build_prefix=${build_base_dir}/build/${arch}-${compiler}/
install_prefix=${build_base_dir}/install/${arch}-${compiler}/
source_prefix=${env_prefix}/

pkg=BOUT-dev

build_dir=${build_prefix}/${pkg}
install_dir=${install_prefix}/${pkg}
source_dir=${source_prefix} 

echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}

if [ "$pkg" == "BOUT-dev" ]; then
    echo 'enter BOUT-dev script'
    mkdir -p $build_dir && cd $build_dir
    cmake -C ${source_dir}/${cmake_config} \
        -DCMAKE_INSTALL_PREFIX=$install_dir \
        $source_dir
fi
