cat <<EOF
FROM docker.io/oi4ai/bout-base:openmp-$MPI-mini

# ----------------------------------------------------------------
# Build and install BOUT++
# ----------------------------------------------------------------
RUN useradd boutuser -p boutforever -G wheel
RUN mkdir /opt/bout++ && chown boutuser /opt/bout++

USER boutuser
WORKDIR /home/boutuser

# Checkout submodules now so configure later is fast, and iterating on
# it less painful
RUN git clone https://github.com/dschwoerer/BOUT-dev.git \
 && chown -R boutuser /home/boutuser/BOUT-dev \
 && cd BOUT-dev \
 && git checkout $COMMIT \
 && git -c submodule.recurse=false submodule update --init --recursive


WORKDIR /home/boutuser/BOUT-dev

ENV MPI_BIN=/usr/lib64/$MPI/bin \
    MPI_SYSCONFIG=/etc/$MPI-x86_64 \
    MPI_FORTRAN_MOD_DIR=/usr/lib64/gfortran/modules/$MPI \
    MPI_INCLUDE=/usr/include/$MPI-x86_64 \
    MPI_LIB=/usr/lib64/$MPI/lib \
    MPI_MAN=/usr/share/man/$MPI-x86_64 \
    MPI_PYTHON_SITEARCH=/usr/lib64/python3.9/site-packages/$MPI \
    MPI_PYTHON3_SITEARCH=/usr/lib64/python3.9/site-packages/$MPI \
    MPI_COMPILER=$MPI-x86_64 \
    MPI_SUFFIX=_$MPI \
    MPI_HOME=/usr/lib64/$MPI


RUN export PATH=\$MPI_BIN:\$PATH LD_LIBRARY_PATH=\$MPI_LIB:\$LD_LIBRARY_PATH \
 && export PKG_CONFIG_PATH=\$MPI_LIB/pkgconfig:\$PKG_CONFIG_PATH MANPATH=\$MPI_MAN:\$MANPATH \
 && cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/bout++/ \
          -DBOUT_GENERATE_FIELDOPS=OFF \
          -DBOUT_USE_PETSC=ON -DPETSc_ROOT=/opt/petsc \
          -DBOUT_ENABLE_METRIC_3D=ON -DBOUT_ENABLE_OPENMP=ON -DBOUT_ENABLE_PYTHON=ON \
          -DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=/usr/lib64/$MPI/ -DSUNDIALS_INCLUDE_DIR=/usr/include/$MPI-x86_64/sundials/ \
	  $CMAKE_OPTIONS || (cat /home/boutuser/BOUT-dev/build/CMakeFiles/CMake{Output,Error}.log  ; exit 1)


RUN make -C build -j 2
RUN make -C build install

RUN find /opt/bout++
EOF
