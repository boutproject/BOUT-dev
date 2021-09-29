Backtrace Example
=================

This demonstrates what the exception backtrace looks like when something goes
wrong in a physics model or in BOUT++. Requires both backtrace to be enabled
(done by default) and debug symbols (`--enable-debug` with `configure` or
`-DCMAKE_BUILD_TYPE=Debug` with `CMake` in both BOUT++ _and_ this example).

The output should look something like:

```
...
c is inf
Error encountered
====== Exception path ======
[bt] #10 ./backtrace() [0x40a27e]
_start at /home/abuild/rpmbuild/BUILD/glibc-2.33/csu/../sysdeps/x86_64/start.S:122
[bt] #9 /lib64/libc.so.6(__libc_start_main+0xd5) [0x7fecbfa28b25]
__libc_start_main at /usr/src/debug/glibc-2.33-4.1.x86_64/csu/../csu/libc-start.c:332
[bt] #8 ./backtrace() [0x40a467]
main at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:32 (discriminator 9)
[bt] #7 /path/to/BOUT-dev/build/libbout++.so(_ZN6Solver8setModelEP12PhysicsModel+0xb5) [0x7fecc0ca2e93]
Solver::setModel(PhysicsModel*) at /path/to/BOUT-dev/build/../src/solver/solver.cxx:94
[bt] #6 /path/to/BOUT-dev/build/libbout++.so(_ZN12PhysicsModel10initialiseEP6Solver+0xc0) [0x7fecc0cad594]
PhysicsModel::initialise(Solver*) at /path/to/BOUT-dev/build/../include/bout/physicsmodel.hxx:93 (discriminator 5)
[bt] #5 ./backtrace() [0x40a986]
Backtrace::init(bool) at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:27
[bt] #4 ./backtrace() [0x40a3cf]
f3() at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:19
[bt] #3 ./backtrace() [0x40a3be]
f2(int) at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:15
[bt] #2 ./backtrace() [0x40a386]
f1() at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:13 (discriminator 2)
[bt] #1 ./backtrace(_ZN13BoutExceptionC1IA19_cJEEERKT_DpRKT0_+0xba) [0x40ae16]
BoutException::BoutException<char [19]>(char const (&) [19]) at /path/to/BOUT-dev/build/../include/bout/../boutexception.hxx:28 (discriminator 2)
====== Back trace ======

====== Exception thrown ======
Tomatoes are red?
```
