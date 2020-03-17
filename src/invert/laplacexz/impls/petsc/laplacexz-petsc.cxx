/*
 * LaplaceXZ implementation using PETSc
 *
 * Ben Dudson, October 2015
 *
 * Based on:
 * PETSc example of reusing preconditioner matrix (in Fortran):
 * http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex6f.F.html
 */

#include "laplacexz-petsc.hxx"

#ifdef BOUT_HAS_PETSC  // Requires PETSc

#include <bout/assert.hxx>
#include <bout/sys/timer.hxx>

#include <msg_stack.hxx>
#include <output.hxx>

LaplaceXZpetsc::LaplaceXZpetsc(Mesh *m, Options *opt, const CELL_LOC loc)
  : LaplaceXZ(m, opt, loc), coefs_set(false) {
  /* Constructor: LaplaceXZpetsc
   * Purpose:     - Setting inversion solver options
   *              - Setting the solver method
   *              - Setting the preconditioner
   *              - Allocating memory to the matrix A (not expanded in memory)
   *              - Allocating memory to the vectors x and b in Ax=b
   */
  /* A note about the boundary condition:
   * ====================================
   * NOTE: All boundaries are currently implemented half between grid points
   *
   * We set the inner boundary first in the x-array in Ax=b, and the outer
   * boundary last in the x-array (z is periodic)
   *
   * The problem at the inner boundaries will be set by
   *
   * A_{0,j}x_{0,j} + A_{1,j}x_{1,j} = b_{0,j}
   *
   * Where the indices denote the indicies in the physical grid, and NOT in
   * the matrix. I.e. the first index refers to the x-index in the grid, and
   * the second index refers to the z-index in the grid.
   * The problem at the outer boundaries will be set by
   *
   * A_{xend,j}x_{xend,j} + A_{xend+1,j}x_{xend+1,j} = b_{xend_1,j}
   *
   * Determining the boundary conditions:
   * ====================================
   * Using x_BC at the value at the boundary, and k as the x index
   *
   * 2nd order dirichlet:
   * --------------------
   * 0.5*x_{k,j} + 0.5*x_{k+1,j} = x_BC
   * x_BC will be calculated from the grid points of the input "in" (given to
   * LaplaceXZpetsc::solve) so
   * 0.5*x_{k,j} + 0.5*x_{k+1,j} = 0.5*in_{k,j} + 0.5*in_{k+1,j}
   * which we will write as
   *
   * x_{k,j} + x_{k+1,j} = in_{k,j} + in_{k+1,j}
   *
   * 2nd order neumann:
   * --------------------
   * -(1/dx)*(x_{k,j}) + (1/dx)*x_{k+1,j} = x_BC
   * x_BC will be calculated from the grid points of the input "in" (given to
   * LaplaceXZpetsc::solve) so
   *   -(1/dx)*(x_{k,j}) + (1/dx)*x_{k+1,j}
   * = -(1/dx)*(in_{k,j}) + (1/dx)*in_{k+1,j}
   * which we will write as
   *
   * x_{0,j} - x_{1,j} = in_{0,j} - in_{1,j}
   * - x_{end-1,j} + x_{end,j} = - in_{end-1,j} + in_{end,j}
   *
   * Notice the change of sign for the inner boundary, which ensures that the
   * boundary part of main diagonal in the matrix A remains positive
   *
   * Specification through ghost point:
   * ----------------------------------
   * If options are set to use the ghost points of the input, we simply use
   * x_{k,j} = in_{k,j}
   *
   * Defualt behaviour:
   * =================-
   * If the boundary flags are not:
   * - xin & xout will be set to 2nd order dirichlet as described above
   *
   * Current implementations:
   * ========================
   * INVERT_AC_GRAD - As "2nd order neumann", but with 0 at the RHS
   * INVERT_SET     - As "Specifiaction through ghost points", where "in" is
   *                  given by x0
   * INVERT_RHS     - As "Specifiaction through ghost points", where "in" is
   *                  given by b
   */

  TRACE("LaplaceXZpetsc::LaplaceXZpetsc");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root())["laplacexz"];
  }

  // Getting the boundary flags
  OPTION(opt, inner_boundary_flags, 0);
  OPTION(opt, outer_boundary_flags, 0);
#if CHECK > 0
  // Checking flags are set to something which is not implemented
  // This is done binary (which is possible as each flag is a power of 2)
  if (inner_boundary_flags & ~implemented_boundary_flags) {
    throw BoutException("Attempted to set LaplaceXZ inversion boundary flag that is not "
                        "implemented in petsc_laplace.cxx");
  }
  if (outer_boundary_flags & ~implemented_boundary_flags) {
    throw BoutException("Attempted to set LaplaceXZ inversion boundary flag that is not "
                        "implemented in petsc_laplace.cxx");
  }
  if (localmesh->periodicX) {
    throw BoutException("LaplacePetsc does not work with periodicity in the x direction "
                        "(localmesh->PeriodicX == true). Change boundary conditions or "
                        "use serial-tri or cyclic solver instead");
  }
#endif

  reuse_limit = (*opt)["reuse_limit"]
                    .doc("How many solves can the preconditioner be reused?")
                    .withDefault(100);
  reuse_count = reuse_limit + 1; // So re-calculates first time

  // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
  // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.

  const BoutReal rtol = (*opt)["rtol"].doc("Relative tolerance").withDefault(1e-5);
  const BoutReal atol = (*opt)["atol"]
          .doc("Absolute tolerance. The solution is considered converged if |Ax-b| "
               "< max( rtol * |b| , atol )")
          .withDefault(1e-10);
  const BoutReal dtol = (*opt)["dtol"]
                      .doc("The solution is considered diverged if |Ax-b| > dtol * |b|")
                      .withDefault(1e3);
  const int maxits = (*opt)["maxits"].doc("Maximum iterations").withDefault(100000);

  // Get KSP Solver Type
  const std::string ksptype = (*opt)["ksptype"].doc("KSP solver type").withDefault("gmres");
  
  // Get PC type
  const std::string pctype = (*opt)["pctype"].doc("Preconditioner type").withDefault("none");

  const std::string factor_package = (*opt)["factor_package"]
          .doc("Package to use in preconditioner. Passed to PCFactorSetMatSolver")
          .withDefault("petsc");

  // Get MPI communicator
  MPI_Comm comm = localmesh->getXcomm();

  // Local size
  int localN = (localmesh->xend - localmesh->xstart + 1) * (localmesh->LocalNz);
  if (localmesh->firstX()) {
    localN += localmesh->LocalNz;
  }
  if (localmesh->lastX()) {
    localN += localmesh->LocalNz;
  }

  // Create Vectors
  VecCreate( comm, &xs );
  VecSetSizes( xs, localN, PETSC_DETERMINE );
  VecSetFromOptions( xs );
  VecDuplicate( xs , &bs );

  for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
    YSlice data;

    data.yindex = y;

    // Set size of Matrix on each processor to localN x localN
    MatCreate( comm, &data.MatA );
    MatSetSizes( data.MatA, localN, localN, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions(data.MatA);

    //////////////////////////////////////////////////
    // Pre-allocate space for matrix elements

    PetscInt *d_nnz, *o_nnz;
    PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
    PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );

    for (int i=0;i<localN;i++) {
      // Non-zero elements on this processor
      d_nnz[i] = 9; 
      // Non-zero elements on neighboring processor
      o_nnz[i] = 3;
    }

    // X boundaries
    if (localmesh->firstX()) {
      for (int z = 0; z < localmesh->LocalNz; z++) {
        d_nnz[z] = 2;
      }
    } else {
      // One point on another processor
      for (int z = 0; z < localmesh->LocalNz; z++) {
        d_nnz[z] -= 1;
        o_nnz[z] += 1;
      }
    }

    if (localmesh->lastX()) {
      for (int z = 0; z < localmesh->LocalNz; z++) {
        int ind = localN - (localmesh->LocalNz) + z;
        d_nnz[ind] = 2;
      }
    } else {
      // One point on another processor
      for (int z = 0; z < localmesh->LocalNz; z++) {
        int ind = localN - (localmesh->LocalNz) + z;
        d_nnz[ind] -= 1;
        o_nnz[ind] += 1;
      }
    }

    MatMPIAIJSetPreallocation( data.MatA, 0, d_nnz, 0, o_nnz );
    MatSetUp(data.MatA);
    PetscFree( d_nnz );
    PetscFree( o_nnz );

    //////////////////////////////////////////////////
    // Declare KSP Context
    KSPCreate( comm, &data.ksp );

    // Set KSP type
    KSPSetType( data.ksp, ksptype.c_str() );
    // Set KSP tolerances
    KSPSetTolerances( data.ksp, rtol, atol, dtol, maxits );

    // Set KSP preconditioner
    PC pc;
    KSPGetPC(data.ksp,&pc);
    PCSetType(pc, pctype.c_str());
#if PETSC_VERSION_GE(3,9,0)
    PCFactorSetMatSolverType(pc,factor_package.c_str());
#else
    PCFactorSetMatSolverPackage(pc,factor_package.c_str());
#endif

    KSPSetFromOptions( data.ksp );

    /// Add to slice vector
    slice.push_back(data);
  }

}

LaplaceXZpetsc::~LaplaceXZpetsc() {

  TRACE("LaplaceXZpetsc::~LaplaceXZpetsc");

  PetscBool petsc_is_finalised;
  PetscFinalized(&petsc_is_finalised);

  for (auto &it : slice) {
    MatDestroy(&it.MatA);
    MatDestroy(&it.MatP);

    if (!petsc_is_finalised) {
      // PetscFinalize may already have destroyed this object
      KSPDestroy(&it.ksp);
    }
  }

  VecDestroy(&bs);
  VecDestroy(&xs);
}

void LaplaceXZpetsc::setCoefs(const Field3D &Ain, const Field3D &Bin) {
  /* Function: LaplaceXZpetsc::setCoefs
   * Purpose:  - Set the matrix coefficients in the matrix MatA (member data)
   *             in Ax=b
   *             NOTE: Do not confuse the matrix A with the input coefficient A
   *                   referring to A in div(A grad_perp(B)) + Bf = b
   *
   * Input
   * Ain       - The A coefficient in div(A grad_perp(B)) + Bf = b
   * Bin       - The B coefficient in div(A grad_perp(B)) + Bf = b
   */

  TRACE("LaplaceXZpetsc::setCoefs");

  ASSERT1(Ain.getMesh() == localmesh);
  ASSERT1(Bin.getMesh() == localmesh);
  ASSERT1(Ain.getLocation() == location);
  ASSERT1(Bin.getLocation() == location);

  #if CHECK > 0
    // Checking flags are set to something which is not implemented
    // This is done binary (which is possible as each flag is a power of 2)
    if ( inner_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set LaplaceXZ inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( outer_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set LaplaceXZ inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
  #endif
  Timer timer("invert");
  // Set coefficients

  Field3D A = Ain;
  Field3D B = Bin;

  // Each Y slice is handled as a separate set of matrices and KSP context
  for (auto &it : slice) {
    // Get Y index
    int y = it.yindex;

    int Istart, Iend;
    MatGetOwnershipRange(it.MatA, &Istart, &Iend);

    ////////////////////////////////////////////////
    // Inner X boundary (see note about BC in LaplaceXZ constructor)
    int row = Istart;
    if(localmesh->firstX()) {
      if (inner_boundary_flags & INVERT_AC_GRAD){
        // Neumann 0
        /* NOTE: Sign of the elements are opposite of what one might expect,
         *       see note about BC in LaplaceXZ constructor for more details
         */
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row + (localmesh->LocalNz); // +1 in X
          val = -1.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else if(inner_boundary_flags & INVERT_SET){
        // Setting BC from x0
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row + (localmesh->LocalNz); // +1 in X
          val = 0.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else if(inner_boundary_flags & INVERT_RHS){
        // Setting BC from b
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row + (localmesh->LocalNz); // +1 in X
          val = 0.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else{
        // Default: Dirichlet on inner x boundary
        /* NOTE: Sign of the elements are opposite of what one might expect,
         *       see note about BC in LaplaceXZ constructor for more details
         */
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 0.5;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row + (localmesh->LocalNz); // +1 in X
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
    }

    ////////////////////////////////////////////////
    // Set matrix elements
    //
    // (1/J) d/dx ( A * J * g11 d/dx ) + (1/J) d/dz ( A * J * g33 d/dz ) + B
    Coordinates *coords = localmesh->getCoordinates(location);
    
    for(int x=localmesh->xstart; x <= localmesh->xend; x++) {
      for(int z=0; z < localmesh->LocalNz; z++) {
	// stencil entries
	PetscScalar c, xm, xp, zm, zp, xpzp, xpzm, xmzp, xmzm;
	PetscScalar zpxp, zpxm, zmxp, zmxm;
	  
	// XX component
	  
	// Metrics on x+1/2 boundary
	BoutReal J = 0.5*(coords->J(x,y,z) + coords->J(x+1,y,z));
	BoutReal g11 = 0.5*(coords->g11(x,y,z) + coords->g11(x+1,y,z));
	BoutReal dx = 0.5*(coords->dx(x,y,z) + coords->dx(x+1,y,z));
	BoutReal Acoef = 0.5*(A(x,y,z) + A(x+1,y,z));
	  
	BoutReal val = Acoef * J * g11 / (coords->J(x,y,z) * dx * coords->dx(x,y,z));
	xp = val;
	c  = -val;
	  
	// Metrics on x-1/2 boundary
	J = 0.5*(coords->J(x,y,z) + coords->J(x-1,y,z));
	g11 = 0.5*(coords->g11(x,y,z) + coords->g11(x-1,y,z));
	dx = 0.5*(coords->dx(x,y,z) + coords->dx(x-1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x-1,y,z));
	  
	val = Acoef * J * g11 / (coords->J(x,y,z) * dx * coords->dx(x,y,z));
	xm = val;
	c  -= val;
	  
	// ZZ component
	// Wrap around z-1 and z+1 indices
	int zminus = (z - 1 + (localmesh->LocalNz)) % (localmesh->LocalNz);
	int zplus = (z + 1) % (localmesh->LocalNz);
	  
	J = 0.5*(coords->J(x,y,z) + coords->J(x,y,zplus));
	BoutReal g33 = 0.5*(coords->g33(x,y,z) + coords->g33(x,y,zplus));
	BoutReal dz = coords->dz;
	// Metrics on z+1/2 boundary
	Acoef = 0.5*(A(x,y,z) + A(x,y,zplus));
	  
	val = Acoef * J * g33 / (coords->J(x,y,z) * dz * dz);
	zp = val;
	c  -= val;
	  
	// Metrics on z-1/2 boundary
	J = 0.5*(coords->J(x,y,z) + coords->J(x,y,zminus));
	g33 = 0.5*(coords->g33(x,y,z) + coords->g33(x,y,zminus));
	Acoef = 0.5*(A(x,y,z) + A(x,y,zminus));
	  
	val = Acoef * J * g33 / (coords->J(x,y,z) * dz * dz);
	zm = val;
	c  -= val;
	  
	// XZ components

	// x+1/2, z+1/2 
	J =  0.5*(coords->J(x,y,z) + coords->J(x+1,y,z));
	BoutReal g13 = 0.5*(coords->g13(x,y,z) + coords->g13(x+1,y,z));
	dz = 4.0*coords->dz;
	Acoef = 0.5*(A(x,y,z) + A(x,y,zplus));

	val = Acoef * J * g13 / (coords->J(x,y,z) * dz * coords->dx(x,y,z));
	xpzp = val;
	c -= val;

	// x+1/2, z-1/2
	// J =  0.5*(coords->J(x,y,z) + coords->J(x+1,y,z));
	// g13 = 0.5*(coords->gxz(x,y,z) + coords->gxz(x+1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x,y,zminus));

	val = - Acoef * J * g13 / (coords->J(x,y,z) * dz * coords->dx(x,y,z));
	xpzm = val;
	c -= val;

	// x-1/2, z+1/2
	J =  0.5*(coords->J(x,y,z) + coords->J(x-1,y,z));
	g13 = 0.5*(coords->g13(x,y,z) + coords->g13(x-1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x,y,zplus));

	val = - Acoef * J * g13 / (coords->J(x,y,z) * dz * coords->dx(x,y,z));
	xmzp = val;
	c -= val;

	// x-1/2, z-1/2
	Acoef = 0.5*(A(x,y,z) + A(x,y,zminus));

	val = Acoef * J * g13 / (coords->J(x,y,z) * dz * coords->dx(x,y,z));
	xmzm = val;
	c -= val;


	// ZX components
	// z+1/2, x+1/2
	J = 0.5*(coords->J(x,y,z) + coords->J(x,y,zplus));
	g13 = 0.5*(coords->g13(x,y,z) + coords->g13(x,y,zplus));
	dx = 2.0*(coords->dx(x,y,z) + coords->dx(x+1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x+1,y,z));

	val = Acoef * J * g13 / (coords->J(x,y,z) * dx * coords->dz);
	zpxp = val;
	c -= val;

	//z+1/2, x-1/2
	dx = 2.0*(coords->dx(x,y,z) + coords->dx(x-1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x-1,y,z));

	val = - Acoef * J * g13 / (coords->J(x,y,z) * dx * coords->dz);
	zpxm = val;
	c -= val;

	// z-1/2, x+1/2
	J = 0.5*(coords->J(x,y,z) + coords->J(x,y,zminus));
	g13 = 0.5*(coords->g13(x,y,z) + coords->g13(x,y,zminus));
	dx = 2.0*(coords->dx(x,y,z) + coords->dx(x+1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x+1,y,z));

	val = - Acoef * J * g13 / (coords->J(x,y,z) * dx * coords->dz);
	zmxp = val;
	c -= val;

	// z-1/2, x-1/2
	dx = 2.0*(coords->dx(x,y,z) + coords->dx(x-1,y,z));
	Acoef = 0.5*(A(x,y,z) + A(x-1,y,z));

	val = Acoef * J * g13 / (coords->J(x,y,z) * dx * coords->dz);
	zmxm = val;
	c -= val;

 	  
	///////////// OLDER CODE /////////////////
	// Note that because metrics are constant in Z many terms cancel
	  
	// // Wrap around z-1 and z+1 indices
	// // int zminus = (z - 1 + (localmesh->LocalNz)) % (localmesh->LocalNz);
	// // int zplus = (z + 1) % (localmesh->LocalNz);
	  
	// // Metrics on z+1/2 boundary
	// Acoef = 0.5*(A(x,y,z) + A(x,y,z+1));
	  
	// val = Acoef * coords->gzz(x,y,z) / (coords->dz(x,y,z)*coords->dz(x,y,z));
	// zp = val;
	// c -= val;
	  
	// // Metrics on z-1/2 boundary
	// Acoef = 0.5*(A(x,y,z) + A(x,y,z-1));
	  
	// val = Acoef * coords->gzz(x,y,z) / (coords->dz(x,y,z)*coords->dz(x,y,z));
	// zm = val;
	// c -= val;
	//////////////////////////////////////////
	  
	  
	// B term
	c += B(x,y,z);
	  
	/////////////////////////////////////////////////
	// Now have a 5-point stencil for the Laplacian
	  
	// Set the centre (diagonal)
	MatSetValues(it.MatA,1,&row,1,&row,&c,INSERT_VALUES);
	  
	// X + 1
	int col = row + (localmesh->LocalNz);
	MatSetValues(it.MatA,1,&row,1,&col,&xp,INSERT_VALUES);
	  
	// X - 1
	col = row - (localmesh->LocalNz);
	MatSetValues(it.MatA,1,&row,1,&col,&xm,INSERT_VALUES);
	  
	// Z + 1
	col = row + 1;
	if(z == localmesh->LocalNz-1) {
	  col -= localmesh->LocalNz;  // Wrap around
	}
	MatSetValues(it.MatA,1,&row,1,&col,&zp,INSERT_VALUES);

	int xzcol = col + (localmesh->LocalNz); // X+1
	MatSetValues(it.MatA,1,&row,1,&xzcol,&xpzp,INSERT_VALUES);
	MatSetValues(it.MatA,1,&row,1,&xzcol,&zpxp,INSERT_VALUES);

	xzcol = col - (localmesh->LocalNz); // X-1
	MatSetValues(it.MatA,1,&row,1,&xzcol,&xmzp,INSERT_VALUES);
	MatSetValues(it.MatA,1,&row,1,&xzcol,&zpxm,INSERT_VALUES);
	  
	// Z - 1
	col = row - 1;
	if(z == 0) {
	  col += localmesh->LocalNz;  // Wrap around
	}
	MatSetValues(it.MatA,1,&row,1,&col,&zm,INSERT_VALUES);

	xzcol = col + (localmesh->LocalNz); // X+1
	MatSetValues(it.MatA,1,&row,1,&xzcol,&xpzm,INSERT_VALUES);
	MatSetValues(it.MatA,1,&row,1,&xzcol,&zmxp,INSERT_VALUES);

	xzcol = col - (localmesh->LocalNz); // X-1
	MatSetValues(it.MatA,1,&row,1,&xzcol,&xmzm,INSERT_VALUES);
	MatSetValues(it.MatA,1,&row,1,&xzcol,&zmxm,INSERT_VALUES);
	  
	row++;
      }
    }
    
    ////////////////////////////////////////////////
    // Outer X boundary (see note about BC in LaplaceXZ constructor)
    if(localmesh->lastX()) {
      if (outer_boundary_flags & INVERT_AC_GRAD){
        // Neumann 0
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row - (localmesh->LocalNz); // -1 in X
          val = -1.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else if (outer_boundary_flags & INVERT_SET){
        // Setting BC from x0
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row - (localmesh->LocalNz); // -1 in X
          val = 0.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else if (outer_boundary_flags & INVERT_RHS){
        // Setting BC from b
        for(int z=0; z < localmesh->LocalNz; z++) {
          PetscScalar val = 1.0;
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row - (localmesh->LocalNz); // -1 in X
          val = 0.0;
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
      else{
        //Default: Dirichlet on outer X boundary
        PetscScalar val = 0.5;

        for(int z=0; z < localmesh->LocalNz; z++) {
          MatSetValues(it.MatA, 1, &row, 1, &row, &val, INSERT_VALUES);

          int col = row - (localmesh->LocalNz); // -1 in X
          MatSetValues(it.MatA, 1, &row, 1, &col, &val, INSERT_VALUES);

          row++;
        }
      }
    }

    ASSERT1(row == Iend); // Check that row is currently on the last row

    // Assemble Matrix
    MatAssemblyBegin(it.MatA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(it.MatA, MAT_FINAL_ASSEMBLY);
  }

  // Increase reuse count
  reuse_count++;
  if(reuse_count > reuse_limit) {
    // Reuse limit exceeded. Reset count
    reuse_count = 0;

    // Modifying preconditioner matrix
    for (auto &it : slice) {
      // Copy matrix into preconditioner
      if(coefs_set) {
        // Preconditioner already set
        MatDestroy(&it.MatP);
      }
      MatConvert(it.MatA, MATSAME, MAT_INITIAL_MATRIX, &it.MatP);
    }

    // Set operators
    for (auto &it : slice) {

      // Note: This is a hack to force update of the preconditioner matrix
#if PETSC_VERSION_GE(3,5,0)
      KSPSetOperators(it.ksp, it.MatA, it.MatP);
#else
      KSPSetOperators(it.ksp, it.MatA, it.MatP, SAME_NONZERO_PATTERN);
#endif
    }
  }else {
    for (auto &it : slice) {
      /// Reuse the preconditioner, even if the operator changes

#if PETSC_VERSION_GE(3,5,0)
      KSPSetReusePreconditioner(it.ksp, PETSC_TRUE);
#else
      KSPSetOperators(it.ksp, it.MatA, it.MatP, SAME_PRECONDITIONER);
#endif
    }
  }

  coefs_set = true;
}

Field3D LaplaceXZpetsc::solve(const Field3D &bin, const Field3D &x0in) {
  /* Function: LaplaceXZpetsc::solve
   * Purpose:  - Set the values of b in  Ax=b
   *           - Set the initial guess x0, and use this for x in  Ax=b
   *           - Solve Ax=b for x
   *           - Recast x to a Field3D
   *
   * Input
   * bin       - The b to be used in Ax=b
   * x0in      - The initial guess x0 to be used in Ax=b
   *
   * Output
   * result    - The solved x (returned as a Field3D) in the matrix problem Ax=b
   */

  TRACE("LaplaceXZpetsc::solve");

  ASSERT1(bin.getMesh() == localmesh);
  ASSERT1(x0in.getMesh() == localmesh);
  ASSERT1(bin.getLocation() == location);
  ASSERT1(x0in.getLocation() == location);

  if(!coefs_set) {
    throw BoutException("LaplaceXZpetsc: solve called before setCoefs");
  }

  Timer timer("invert");

  Field3D b = bin;
  Field3D x0 = x0in;

  Field3D result{emptyFrom(bin)};

  for (auto &it : slice) {
    /// Get y index
    int y = it.yindex;

    /// Specify non-zero starting guess for solution (from input x0)
    KSPSetInitialGuessNonzero(it.ksp, PETSC_TRUE);

    //////////////////////////
    // Load initial guess x0 into xs and rhs into bs
    int Istart, Iend;
    MatGetOwnershipRange(it.MatA, &Istart, &Iend);

    // Starting index
    int ind = Istart;

    // Inner X boundary (see note about BC in LaplaceXZ constructor)
    if(localmesh->firstX()) {
      if (inner_boundary_flags & INVERT_AC_GRAD){
        // Neumann 0
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xstart-1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = 0.0;
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
          ind++;
        }
      }
      else if (inner_boundary_flags & INVERT_SET){
        // Setting BC from x0
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xstart-1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = x0(localmesh->xstart,y,z);
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
          ind++;
        }
      }
      else if (inner_boundary_flags & INVERT_RHS){
        // Setting BC from b
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xstart-1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = b(localmesh->xstart,y,z);
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
          ind++;
        }
      }
      else{
        // Default: Neumann on inner x boundary
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xstart-1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = x0(localmesh->xstart-1,y,z) - x0(localmesh->xstart,y,z);
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
          ind++;
        }
      }
    }

    // Set the inner points
    for(int x=localmesh->xstart;x<= localmesh->xend;x++) {
      for(int z=0; z < localmesh->LocalNz; z++) {
	PetscScalar val = x0(x,y,z);
	VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
	  
	val = b(x,y,z);
	VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
	ind++;
      }
    }

    // Outer X boundary (see note about BC in LaplaceXZ constructor)
    if(localmesh->lastX()) {
      if (outer_boundary_flags & INVERT_AC_GRAD){
        // Neumann 0
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xend+1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = 0.0;
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );

          ind++;
        }
      }
      else if (outer_boundary_flags & INVERT_SET){
        // Setting BC from x0
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xend+1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = x0(localmesh->xend+1,y,z);
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );

          ind++;
        }
      }
      else if (outer_boundary_flags & INVERT_RHS){
        // Setting BC from b
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xend+1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = b(localmesh->xend+1,y,z);
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );

          ind++;
        }
      }
      else{
        //Default: Dirichlet on outer X boundary
        for(int z=0; z < localmesh->LocalNz; z++) {
          // Setting the initial guess x0
          PetscScalar val = x0(localmesh->xend+1,y,z);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );

          // Setting the solution b
          val = 0.5*(x0(localmesh->xend,y,z) + x0(localmesh->xend+1,y,z));
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );

          ind++;
        }
      }
    }

    ASSERT1(ind == Iend); // Reached end of range

    // Assemble RHS Vector
    VecAssemblyBegin(bs);
    VecAssemblyEnd(bs);

    // Assemble Trial Solution Vector
    VecAssemblyBegin(xs);
    VecAssemblyEnd(xs);

    //////////////////////////
    // Solve the system

    KSPSolve(it.ksp, bs, xs);

    // Check if the solve converged
    KSPConvergedReason reason;
    KSPGetConvergedReason(it.ksp, &reason);

    if(reason <= 0) {
      throw BoutException("LaplaceXZ failed to converge. Reason %d", reason);
    }

    //////////////////////////
    // Copy data into result

    ind = Istart;
    // Inner X boundary
    if(localmesh->firstX()) {
      for(int z=0; z < localmesh->LocalNz; z++) {
        PetscScalar val;
        VecGetValues(xs, 1, &ind, &val );
        result(localmesh->xstart-1,y,z) = val;
        ind++;
      }
    }

    for(int x=localmesh->xstart;x<= localmesh->xend;x++) {
      for(int z=0; z < localmesh->LocalNz; z++) {
        PetscScalar val;
        VecGetValues(xs, 1, &ind, &val );
        result(x,y,z) = val;
        ind++;
      }
    }

    // Outer X boundary
    if(localmesh->lastX()) {
      for(int z=0; z < localmesh->LocalNz; z++) {
        PetscScalar val;
        VecGetValues(xs, 1, &ind, &val );
        result(localmesh->xend+1,y,z) = val;
        ind++;
      }
    }
    ASSERT1(ind == Iend); // Reached end of range
  }

  return result;
}

#endif // BOUT_HAS_PETSC
