# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install bout-dev
#
# You can edit this file again by typing:
#
#     spack edit bout-dev
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack.package import *


class BoutDev(CMakePackage):
    """
    BOUT++ is a framework for writing fluid and plasma simulations in
    curvilinear geometry. It is intended to be quite modular, with a variety of
    numerical methods and time-integration solvers available. BOUT++ is
    primarily designed and tested with reduced plasma fluid models in mind, but
    it can evolve any number of equations, with equations appearing in a
    readable form.
    """

    homepage = "https://boutproject.github.io/"
    url      = "https://github.com/boutproject/BOUT-dev/releases/download/v4.4.2/BOUT++-v4.4.2.tar.gz"

    maintainers = ['JosephThomasParker', 'ZedThree', 'bendudson']

    version('4.4.2', sha256='bb448c74f3c0b8cb5568aa9ac641d708065552bcd1a623f9fc9649cec2648e63')
    version('4.4.1', sha256='a2bdd70597a68edd8beb5259c9d3fa697bafa0e94b232134a0c79fade43a3ac2')
    version('4.4.0', sha256='a71125af728686450d0f9d65081d2b3112d113dcdacb22ebd776cd84f5f720fe')
    version('4.3.3', sha256='61912db1f09a7df0437b388b9fe9777c3652626a197d1899ad140c49f08ec813')
    version('4.3.2', sha256='ee10aa9ce90fdec4d3ac0cc11bfdedadcbda19b3d0a74851f94633ddc80c5751')
    version('4.3.1', sha256='7763fb4be96dd89551a0bb3a9b435dc09ebac4ef26d80b5edaa0f7c013d1b558')
    version('4.3.0', sha256='db50a66a62edf87f04b8cb6637838811bb9307726e748a9c1979fb1cbf250cd9')
    version('4.2.3', sha256='321192c8297855be65fad8080ed880eb50338db7fff680614d430f7fef750f42')
    version('4.2.2', sha256='80518ccf8aa4d2e61fd1da96b427725ef03e481ee88d193e39eb846c5f292290')
    version('4.2.1', sha256='2e757fc15d1d52bfb80f94963d1aac95a5e910c60ec586e13ee107694a948bf6')

    # Variants
    variant('shared', default=False, description='Enable building bout++ into an shared object')
    variant('checks', default='none', description='Set run-time checking level', values=('1', '2', '3', 'none'), multi=False)
    variant('optimize', default=False, description='Enable optimization')
    variant('track', default=False, description='Enable variable tracking')
    variant('debug', default=False, description='Enable all debugging flags')
    variant('nls', default=True, description='Enable Native Language support')
    variant('openmp', default=False, description='Enable building with OpenMP support')
    variant('pvode-openmp', default=False, description='Enable building PVODE with OpenMP support')
    variant('python', default=True, description='Activate Python packages for data manipulation')
    variant('sundials', default=False, description='Use SUNDIALS CVODE, IDA, and ARKODE solvers')
    variant('blaslapack', default=False, description='Use BLAS and LAPACK')
    variant('petsc', default=False, description='Use PETSc interface')
    variant('slepc', default=False, description='Use SLEPc interface')
    variant('mumps', default=False, description='Use MUMPS library for direct matrix inversions')
    variant('scorep', default=False, description='Enable support for scorep based instrumentation')

    # Mandatory dependencies
    depends_on('mpi@3:', type=('build', 'link', 'run'))
    depends_on('fftw@3:')

    # Optional dependencies
    depends_on('hdf5 +mpi')
    depends_on('netcdf-c +mpi')
    depends_on('netcdf-cxx4')
    depends_on('python@3:', when='+python', type='run')
    depends_on('py-netcdf4', when='+python', type='run')
    depends_on('sundials +CVODE+IDA+ARKODE', when='+sundials')
    depends_on('netlib-lapack~external-blas', when='+blaslapack')
    depends_on('petsc', when='+petsc')
    depends_on('slepc', when='+slepc')
    depends_on('mumps', when='+mumps')
    depends_on('scorep', when='+scorep')

    def cmake_args(self):
        # FIXME: Add arguments other than
        # FIXME: CMAKE_INSTALL_PREFIX and CMAKE_BUILD_TYPE
        # FIXME: If not needed delete this function

        spec = self.spec
        args = []

        if spec.satisfies('+shared'):
            args.append('-DBUILD_SHARED_LIBS=on')

        if spec.satisfies('checks=none'):
            args.append('-DCHECK=0')
        else:
            args.append('-DCHECK={}'.format(spec.variants['checks'].value))

        if spec.satisfies('+optimize'):
            args.append('-DCMAKE_BUILD_TYPE=Release')

        # Profiling and debugging
        if spec.satisfies('+track'):
            args.append('-DBOUT_ENABLE_TRACK=on')

        if spec.satisfies('+debug'):
            #args.append('-DBOUT_ENABLE_OUTPUT_DEBUG=on')
            args.append('-DCMAKE_BUILD_TYPE=Debug')

        if spec.satisfies('+scorep'):
            args.append('-DBOUT_USE_SCOREP={}'.format(spec['scorep'].prefix.bin))

        # Native Language support
        if spec.satisfies('~nls'):
            args.append('-DBOUT_USE_NLS=OFF')

        # OpenMP support
        if spec.satisfies('+openmp'):
            args.append('-DBOUT_ENABLE_OPENMP=ON')
        else:
            args.append('-DBOUT_ENABLE_OPENMP=OFF')

        if spec.satisfies('+pvode-openmp'):
            args.append('-DBOUT_ENABLE_OPENMP=ON')
            args.append('-DBOUT_USE_PVODE=ON')

#        # File format (only the parallel versions)
#        #args.append('--with-hdf5={}/h5cc'.format(spec['hdf5'].prefix.bin))
#        args.append('--with-parallelhdf5={}/h5pcc'.format(spec['hdf5'].prefix.bin))
        args.append('-DBOUT_USE_NETCDF={}'.format(spec['netcdf-cxx4'].prefix))
#        #args.append('--with-pnetcdf={}'.format(spec['netcdf-cxx4'].prefix))
#
        # Math libraries
        args.append('-DBOUT_USE_FFTW={}'.format(spec['fftw'].prefix))

        if spec.satisfies('+sundials'):
            args.append('-DBOUT_USE_SUNDIALS={}'.format(spec['sundials'].prefix))

        if spec.satisfies('+blaslapack'):
            args.append('-DBOUT_USE_LAPACK={}'.format(spec['netlib-lapack'].prefix))

        if spec.satisfies('+petsc'):
            args.append('-DBOUT_USE_PETSC={}'.format(spec['petsc'].prefix))

        if spec.satisfies('+slepc'):
            args.append('-DBOUT_USE_SLEPC={}'.format(spec['slepc'].prefix))

        #if spec.satisfies('+mumps'):
        #    args.append('--with-mumps={}'.format(spec['mumps'].prefix))

        return args
