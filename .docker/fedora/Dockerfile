ARG BASE
FROM ghcr.io/dschwoerer/bout-container-base:$BASE

# ----------------------------------------------------------------
# Build and install BOUT++
# ----------------------------------------------------------------
# user: boutuser
# password: boutforever

ARG MPI
ARG URL
ARG COMMIT
ARG CMAKE_OPTIONS

RUN sudo ls

# Checkout submodules now so configure later is fast, and iterating on
# it less painful
RUN git clone $URL \
 && chown -R boutuser /home/boutuser/BOUT-dev \
 && cd BOUT-dev \
 && git checkout $COMMIT \
 && git submodule update --init --recursive


WORKDIR /home/boutuser/BOUT-dev

RUN cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/bout++/ \
          -DBOUT_GENERATE_FIELDOPS=OFF \
          -DBOUT_USE_PETSC=ON -DPETSc_ROOT=/opt/petsc \
          -DBOUT_ENABLE_PYTHON=ON \
          -DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=/usr/lib64/$MPI/ -DSUNDIALS_INCLUDE_DIR=/usr/include/$MPI-x86_64/sundials/ \
	  $CMAKE_OPTIONS || (cat /home/boutuser/BOUT-dev/build/CMakeFiles/CMake{Output,Error}.log  ; exit 1)


RUN make -C build -j 2
RUN sudo make -C build install

RUN find /opt/bout++
