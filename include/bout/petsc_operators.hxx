/// Represent operators using PETSc matrices
///

#include "bout/mesh.hxx"
#include <bout/petsclib.hxx>

#include <petscmat.h>

class PetscOperators {
public:
  PetscOperators(Mesh* mesh);

  // Static functions for unit testing
  static Region<Ind3D> create_region(Field3D cell_number);
  static Region<Ind3D> create_region_xin(Field3D cell_number);
  static Region<Ind3D> create_region_xout(Field3D cell_number);

private:
  PetscLib lib{}; // Initialize and finalize PETSc

  Region<Ind3D> evolving_region; ///< Evolving cells
  Region<Ind3D> xin_region;      ///< X boundary cells in inner boundary
  Region<Ind3D> xout_region;     ///< X boundary cells in outer boundary
  Region<Ind3D> yup_region;      ///< Y boundary cells in forward direction
  Region<Ind3D> ydown_region;    ///< Y boundary cells in backward direction

  Mat mat_mesh_to_petsc;

  Mat mat_yforward;
};
