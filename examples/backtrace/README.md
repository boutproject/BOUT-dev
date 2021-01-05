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
[bt] #11 ./build/backtrace() [0x42ee5e]
_start at ??:?
[bt] #10 /lib64/libc.so.6(__libc_start_main+0xf3) [0x7f7d8bfc11a3]
__libc_start_main at ??:?
[bt] #9 ./build/backtrace() [0x42f08f]
main at /path/to/BOUT++/include/boutmain.hxx:91 (discriminator 9)
[bt] #8 /path/to/BOUT++/build/libbout++.so(_ZN6Solver8setModelEP12PhysicsModel+0xb5) [0x7f7d8d1ce121]
Solver::setModel(PhysicsModel*) at /path/to/BOUT++/build/../src/solver/solver.cxx:92
[bt] #7 /path/to/BOUT++/build/libbout++.so(_ZN12PhysicsModel10initialiseEP6Solver+0xbe) [0x7f7d8d1d6bf6]
PhysicsModel::initialise(Solver*) at /path/to/BOUT++/build/../include/bout/physicsmodel.hxx:79 (discriminator 5)
[bt] #6 ./build/backtrace() [0x433ceb]
LegacyModel::init(bool) at /path/to/BOUT++/include/boutmain.hxx:63
[bt] #5 ./build/backtrace() [0x42f2f6]
physics_init(bool) at /path/to/BOUT++/examples/backtrace/backtrace.cxx:28
[bt] #4 ./build/backtrace() [0x42f2dd]
f3() at /path/to/BOUT++/examples/backtrace/backtrace.cxx:22
[bt] #3 ./build/backtrace() [0x42f2cc]
f2(int) at /path/to/BOUT++/examples/backtrace/backtrace.cxx:18
[bt] #2 ./build/backtrace() [0x42f294]
f1() at /path/to/BOUT++/examples/backtrace/backtrace.cxx:14 (discriminator 2)
[bt] #1 ./build/backtrace(_ZN13BoutExceptionC1IA19_cJEEERKT_DpRKT0_+0x37) [0x4343c5]
BoutException::BoutException<char [19]>(char const (&) [19]) at /path/to/BOUT++/include/boutexception.hxx:26 (discriminator 2)
====== Back trace ======

====== Exception thrown ======
Tomatoes are red?
```
