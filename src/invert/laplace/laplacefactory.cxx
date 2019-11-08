#include "laplacefactory.hxx"

#include "impls/serial_tri/serial_tri.hxx"
#include "impls/serial_band/serial_band.hxx"
#include "impls/pdd/pdd.hxx"
#include "impls/spt/spt.hxx"
#include "impls/petsc/petsc_laplace.hxx"
#include "impls/mumps/mumps_laplace.hxx"
#include "impls/cyclic/cyclic_laplace.hxx"
#include "impls/shoot/shoot_laplace.hxx"
#include "impls/multigrid/multigrid_laplace.hxx"
#include "impls/naulin/naulin_laplace.hxx"

std::string StandardFactoryTraits<Laplacian>::getDefaultType() {
  return LAPLACE_CYCLIC;
}
