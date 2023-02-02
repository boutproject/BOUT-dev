// Copied from mpi4py, mpi4py/src/lib-mpi/compat/openmpi.h and
// mpi4py/src/dynload.h (at the time of writing on 14/5/2018, mpi4py was at
// revision 18c52a9):
//     https://bitbucket.org/mpi4py/mpi4py
//     Author:  Lisandro Dalcin
//     Contact: dalcinl@gmail.com
//
//     Copyright (c) 2017, Lisandro Dalcin. All rights reserved.
//
//     Redistribution and use in source and binary forms, with or without modification,
//     are permitted provided that the following conditions are met:
//
//     Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer. Redistributions in binary form must
//     reproduce the above copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided with the
//     distribution.
//
//     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS"
//     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//     ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//     INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//     CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//     ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//     POSSIBILITY OF SUCH DAMAGE.

#ifndef PyMPI_DYNLOAD_H
#define PyMPI_DYNLOAD_H

#if defined(OMPI_MAJOR_VERSION)

bool is_openmpi = true;

#if HAVE_DLFCN_H
#include <dlfcn.h>
#else
#if defined(__linux) || defined(__linux__)
#define RTLD_LAZY 0x00001
#define RTLD_NOW 0x00002
#define RTLD_LOCAL 0x00000
#define RTLD_GLOBAL 0x00100
#define RTLD_NOLOAD 0x00004
#define RTLD_NODELETE 0x01000
#define RTLD_DEEPBIND 0x00008
#elif defined(__sun) || defined(__sun__)
#define RTLD_LAZY 0x00001
#define RTLD_NOW 0x00002
#define RTLD_LOCAL 0x00000
#define RTLD_GLOBAL 0x00100
#define RTLD_NOLOAD 0x00004
#define RTLD_NODELETE 0x01000
#define RTLD_FIRST 0x02000
#elif defined(__APPLE__)
#define RTLD_LAZY 0x1
#define RTLD_NOW 0x2
#define RTLD_LOCAL 0x4
#define RTLD_GLOBAL 0x8
#define RTLD_NOLOAD 0x10
#define RTLD_NODELETE 0x80
#define RTLD_FIRST 0x100
#elif defined(__CYGWIN__)
#define RTLD_LAZY 1
#define RTLD_NOW 2
#define RTLD_LOCAL 0
#define RTLD_GLOBAL 4
#else
#error do not recognize system type
#endif
extern "C" {
extern void* dlopen(const char*, int);
extern void* dlsym(void*, const char*);
extern int dlclose(void*);
extern char* dlerror(void);
}
#endif

#ifndef RTLD_LAZY
#define RTLD_LAZY 1
#endif
#ifndef RTLD_NOW
#define RTLD_NOW RTLD_LAZY
#endif
#ifndef RTLD_LOCAL
#define RTLD_LOCAL 0
#endif
#ifndef RTLD_GLOBAL
#define RTLD_GLOBAL RTLD_LOCAL
#endif

static void PyMPI_OPENMPI_dlopen_libmpi(void) {
  void* handle = 0;
  int mode = RTLD_NOW | RTLD_GLOBAL;
#if defined(__APPLE__)
/* macOS */
#ifdef RTLD_NOLOAD
  mode |= RTLD_NOLOAD;
#endif
#if defined(OMPI_MAJOR_VERSION)
#if OMPI_MAJOR_VERSION == 3
  if (!handle)
    handle = dlopen("libmpi.40.dylib", mode);
#elif OMPI_MAJOR_VERSION == 2
  if (!handle)
    handle = dlopen("libmpi.20.dylib", mode);
#elif OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 10
  if (!handle)
    handle = dlopen("libmpi.12.dylib", mode);
#elif OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 6
  if (!handle)
    handle = dlopen("libmpi.1.dylib", mode);
#elif OMPI_MAJOR_VERSION == 1
  if (!handle)
    handle = dlopen("libmpi.0.dylib", mode);
#endif
#endif
  if (!handle)
    handle = dlopen("libmpi.dylib", mode);
#else
/* GNU/Linux and others */
#ifdef RTLD_NOLOAD
  mode |= RTLD_NOLOAD;
#endif
#if defined(OMPI_MAJOR_VERSION)
#if OMPI_MAJOR_VERSION >= 10 /* IBM Spectrum MPI */
  if (!handle)
    handle = dlopen("libmpi_ibm.so.2", mode);
  if (!handle)
    handle = dlopen("libmpi_ibm.so.1", mode);
  if (!handle)
    handle = dlopen("libmpi_ibm.so", mode);
#elif OMPI_MAJOR_VERSION == 3
  if (!handle)
    handle = dlopen("libmpi.so.40", mode);
#elif OMPI_MAJOR_VERSION == 2
  if (!handle)
    handle = dlopen("libmpi.so.20", mode);
#elif OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 10
  if (!handle)
    handle = dlopen("libmpi.so.12", mode);
#elif OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 6
  if (!handle)
    handle = dlopen("libmpi.so.1", mode);
#elif OMPI_MAJOR_VERSION == 1
  if (!handle)
    handle = dlopen("libmpi.so.0", mode);
#endif
#endif
  if (!handle)
    handle = dlopen("libmpi.so", mode);
#endif
}

#else /* !defined(OMPI_MAJOR_VERSION) */

bool is_openmpi = false;
static void PyMPI_OPENMPI_dlopen_libmpi(void) {}

#endif /* defined(OMPI_MAJOR_VERSION) */

#endif /* !PyMPI_DYNLOAD_H */
