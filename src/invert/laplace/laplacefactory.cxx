
#include <globals.hxx>
#include <boutexception.hxx>
#include <strings.h>

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
#include "impls/petscamg/petscamg.hxx"
#include "impls/petsc3damg/petsc3damg.hxx"

#define LAPLACE_SPT  "spt"
#define LAPLACE_PDD  "pdd"
#define LAPLACE_TRI  "tri"
#define LAPLACE_BAND "band"
#define LAPLACE_PETSC "petsc"
#define LAPLACE_MUMPS "mumps"
#define LAPLACE_CYCLIC "cyclic"
#define LAPLACE_SHOOT "shoot"
#define LAPLACE_MULTIGRID "multigrid"
#define LAPLACE_NAULIN "naulin"
#define LAPLACE_PETSCAMG "petscamg"
#define LAPLACE_PETSC3DAMG "petsc3damg"

LaplaceFactory *LaplaceFactory::instance = nullptr;

LaplaceFactory* LaplaceFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new LaplaceFactory();
  }
  return instance;
}

Laplacian* LaplaceFactory::createLaplacian(Options *options, const CELL_LOC loc, Mesh *mesh_in) {
  if (options == nullptr)
    options = Options::getRoot()->getSection("laplace");

  if (mesh_in == nullptr) {
    mesh_in = bout::globals::mesh;
  }

  std::string type;

  if(mesh_in->firstX() && mesh_in->lastX()) {
    // Can use serial algorithm

    options->get("type", type, LAPLACE_CYCLIC);

    if(strcasecmp(type.c_str(), LAPLACE_TRI) == 0) {
      return new LaplaceSerialTri(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_BAND) == 0) {
      return new LaplaceSerialBand(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_SPT) == 0) {
      return new LaplaceSPT(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_PETSC) == 0) {
      return new LaplacePetsc(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_MUMPS) == 0) {
      return new LaplaceMumps(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_CYCLIC) == 0) {
      return new LaplaceCyclic(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_SHOOT) == 0) {
      return new LaplaceShoot(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_MULTIGRID) == 0) {
      return new LaplaceMultigrid(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_NAULIN) == 0) {
      return new LaplaceNaulin(options, loc, mesh_in);
    }else if(strcasecmp(type.c_str(), LAPLACE_PETSCAMG) == 0) {
      return new LaplacePetscAmg(options);
    }else if(strcasecmp(type.c_str(), LAPLACE_PETSC3DAMG) == 0) {
      return new LaplacePetsc3dAmg(options);
    }else {
      throw BoutException("Unknown serial Laplacian solver type '%s'", type.c_str());
    }
  }

  options->get("type", type, LAPLACE_CYCLIC);

  // Parallel algorithm
  if(strcasecmp(type.c_str(), LAPLACE_PDD) == 0) {
    return new LaplacePDD(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_SPT) == 0) {
    return new LaplaceSPT(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_PETSC) == 0) {
    return new LaplacePetsc(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_MUMPS) == 0) {
    return new LaplaceMumps(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_CYCLIC) == 0) {
    return new LaplaceCyclic(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_SHOOT) == 0) {
    return new LaplaceShoot(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_MULTIGRID) == 0) {
      return new LaplaceMultigrid(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_NAULIN) == 0) {
    return new LaplaceNaulin(options, loc, mesh_in);
  }else if(strcasecmp(type.c_str(), LAPLACE_PETSCAMG) == 0) {
      return new LaplacePetscAmg(options);
  }else if(strcasecmp(type.c_str(), LAPLACE_PETSC3DAMG) == 0) {
      return new LaplacePetsc3dAmg(options);
  }else {
    throw BoutException("Unknown parallel Laplacian solver type '%s'", type.c_str());
  }
}

