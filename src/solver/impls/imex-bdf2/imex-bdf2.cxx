#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC

#include "imex-bdf2.hxx"

#include <bout/mesh.hxx>
#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/assert.hxx>

#include <cmath>

#include <output.hxx>

#include "petscsnes.h"
#include "petscmat.h"

// Redundent definition because < C++17
constexpr int IMEXBDF2::MAX_SUPPORTED_ORDER;

IMEXBDF2::IMEXBDF2(Options *opt)
    : Solver(opt), snes_f(nullptr), snes_x(nullptr), snes(nullptr), snesAlt(nullptr),
      snesUse(nullptr), Jmf(nullptr) {

  has_constraints = true; ///< This solver can handle constraints
}

IMEXBDF2::~IMEXBDF2() {
  if (snes_f != nullptr) {
    VecDestroy(&snes_f);
  }
  if (snes_x != nullptr) {
    VecDestroy(&snes_x);
  }
}

/*!
 * PETSc callback function, which evaluates the nonlinear
 * function to be solved by SNES.
 *
 * This function assumes the context void pointer is a pointer
 * to an IMEXBDF2 object.
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
static PetscErrorCode FormFunction(SNES UNUSED(snes), Vec x, Vec f, void *ctx) {
  return static_cast<IMEXBDF2 *>(ctx)->snes_function(x, f, false);
}

/*!
 * PETSc callback function for forming Jacobian
 *
 * This function can be a linearised form of FormFunction
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionForDifferencing"
static PetscErrorCode FormFunctionForDifferencing(void* ctx, Vec x, Vec f) {
  return static_cast<IMEXBDF2*>(ctx)->snes_function(x, f, true);
}

/*!
 * SNES callback for forming Jacobian with coloring
 *
 * This can be a linearised and simplified form of FormFunction
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionForColoring"
static PetscErrorCode FormFunctionForColoring(SNES UNUSED(snes), Vec x, Vec f,
                                              void *ctx) {
  return static_cast<IMEXBDF2 *>(ctx)->snes_function(x, f, true);
}

#undef __FUNCT__
#define __FUNCT__ "imexbdf2PCapply"
static PetscErrorCode imexbdf2PCapply(PC pc,Vec x,Vec y) {
  int ierr;

  // Get the context
  IMEXBDF2 *s;
  ierr = PCShellGetContext(pc, reinterpret_cast<void**>(&s));CHKERRQ(ierr);

  PetscFunctionReturn(s->precon(x, y));
}

/*!
 * Initialisation routine. Called once before solve.
 *
 */
int IMEXBDF2::init(int nout, BoutReal tstep) {

  TRACE("Initialising IMEX-BDF2 solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  output << "\n\tIMEX-BDF2 time-integration solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  int ntmp;
  if (bout::globals::mpi->MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
               n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Check if there are any constraints
  have_constraints = false;

  for(int i=0;i<n2Dvars();i++) {
    if(f2d[i].constraint) {
      have_constraints = true;
      break;
    }
  }
  for(int i=0;i<n3Dvars();i++) {
    if(f3d[i].constraint) {
      have_constraints = true;
      break;
    }
  }

  if (have_constraints) {
    is_dae.reallocate(nlocal);
    // Call the Solver function, which sets the array
    // to zero when not a constraint, one for constraint
    set_id(std::begin(is_dae));
  }

  // Get options
  OPTION(options, timestep, tstep); // Internal timestep
  OPTION(options, mxstep, 100000); //Maximum number of internal iterations
  ninternal = static_cast<int>(out_timestep / timestep);
  if((ninternal == 0) || (out_timestep / ninternal > timestep))
    ++ninternal;
  if(ninternal>mxstep){
    throw BoutException("Error: Number of internal timesteps ({:d}) exceeds mxstep ({:d})", ninternal, mxstep);
  };

  timestep = out_timestep / ninternal;
  output.write("\tUsing timestep = {:e}, {:d} internal steps per output\n", timestep, ninternal);


  OPTION(options, maxOrder, 2); //Maximum order of the scheme (1/2/3)
  if(maxOrder > MAX_SUPPORTED_ORDER){
    throw BoutException("Requested maxOrder greater than MAX_SUPPORTED_ORDER ({:d})",MAX_SUPPORTED_ORDER);
  }

  // Allocate memory and initialise structures
  u.reallocate(nlocal);
  for(int i=0;i<maxOrder;i++){
    uV.emplace_back(Array<BoutReal>{nlocal});
    fV.emplace_back(Array<BoutReal>{nlocal});
    timesteps.push_back(timestep);
    uFac.push_back(0.0);
    fFac.push_back(0.0);
    gFac.push_back(0.0);
  }

  rhs.reallocate(nlocal);

  OPTION(options, adaptive, true); //Do we try to estimate the error?
  OPTION(options, nadapt, 4); //How often do we check the error
  OPTION(options, dtMinFatal, 1.0e-10);
  OPTION(options, dtMax, out_timestep);
  OPTION(options, dtMin, dtMinFatal);
  if(adaptive){
    err.reallocate(nlocal);
    OPTION(options, adaptRtol, 1.0e-3); //Target relative error
    OPTION(options, mxstepAdapt, mxstep); //Maximum no. consecutive times we try to reduce timestep
    OPTION(options, scaleCushUp, 1.5);
    OPTION(options, scaleCushDown, 1.0);
  }

  // Put starting values into u
  saveVars(std::begin(u));
  for(int i=0; i<nlocal; i++){
    for (auto& u_: uV) {
      u_[i] = u[i];
    }
  }

  // Initialise PETSc components
  int ierr;

  // Vectors
  ierr = VecCreate(BoutComm::get(), &snes_x);CHKERRQ(ierr);
  ierr = VecSetSizes(snes_x, nlocal, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(snes_x);CHKERRQ(ierr);
  VecDuplicate(snes_x,&snes_f);

  //The SNES solver object(s)
  constructSNES(&snes);
  if(adaptive) constructSNES(&snesAlt);

  return 0;
}

//Set up a snes object stored at the specified location
void IMEXBDF2::constructSNES(SNES *snesIn){

  // Use global mesh for now
  Mesh* mesh = bout::globals::mesh;

  // Nonlinear solver interface (SNES)
  SNESCreate(BoutComm::get(),snesIn);

  // Set the callback function
  SNESSetFunction(*snesIn,snes_f,FormFunction,this);

  /////////////////////////////////////////////////////
  // Set up the Jacobian

  bool matrix_free;
  OPTION(options, matrix_free, true); // Default is matrix free
  if(matrix_free) {
    /*!
      PETSc SNES matrix free Jacobian, using a different
      operator for differencing.

      See PETSc examples
      http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tests/ex7.c.html
      and this thread:
      http://lists.mcs.anl.gov/pipermail/petsc-users/2014-January/020075.html

     */
    MatCreateSNESMF(*snesIn,&Jmf);

    // Set a function to be called for differencing
    // This can be a linearised form of the SNES function
    MatMFFDSetFunction(Jmf,FormFunctionForDifferencing,this);

    // Calculate Jacobian matrix free using FormFunctionForDifferencing
    SNESSetJacobian(*snesIn,Jmf,Jmf,MatMFFDComputeJacobian,this);
  }else {
    /*!
     * Calculate Jacobian using finite differences.
     * 
     */

    bool use_coloring;
    OPTION(options, use_coloring, true);
    if(use_coloring) {
      // Use matrix coloring to calculate Jacobian

      //////////////////////////////////////////////////
      // Get the local indices by starting at 0
      Field3D index = globalIndex(0);

      //////////////////////////////////////////////////
      // Pre-allocate PETSc storage

      int localN = getLocalN(); // Number of rows on this processor
      int n2d = f2d.size();
      int n3d = f3d.size();

      // Set size of Matrix on each processor to localN x localN
      MatCreate( BoutComm::get(), &Jmf );                                
      MatSetSizes( Jmf, localN, localN, PETSC_DETERMINE, PETSC_DETERMINE );
      MatSetFromOptions(Jmf);
      
      PetscInt *d_nnz, *o_nnz;
      PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
      PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );

      // Set values for most points
      if(mesh->LocalNz > 1) {
        // A 3D mesh, so need points in Z

        for(int i=0;i<localN;i++) {
          // Non-zero elements on this processor
          d_nnz[i] = 7*n3d + 5*n2d; // Star pattern in 3D
          // Non-zero elements on neighboring processor
          o_nnz[i] = 0;
        }
      }else {
        // Only one point in Z
        
        for(int i=0;i<localN;i++) {
          // Non-zero elements on this processor
          d_nnz[i] = 5*(n3d+n2d); // Star pattern in 2D
          // Non-zero elements on neighboring processor
          o_nnz[i] = 0;
        }
      }

      // X boundaries
      if(mesh->firstX()) {
        // Lower X boundary
        for(int y=mesh->ystart;y<=mesh->yend;y++) {
          for(int z=0;z<mesh->LocalNz;z++) {
            int localIndex = ROUND(index(mesh->xstart, y, z));
            ASSERT2( (localIndex >= 0) && (localIndex < localN) );
            if(z == 0) {
              // All 2D and 3D fields
              for(int i=0;i<n2d+n3d;i++)
                d_nnz[localIndex + i] -= (n3d + n2d);
            }else {
              // Only 3D fields
              for(int i=0;i<n3d;i++)
                d_nnz[localIndex + i] -= (n3d + n2d);
            }
          }
        }
      }else {
        // On another processor
        for(int y=mesh->ystart;y<=mesh->yend;y++) {
          for(int z=0;z<mesh->LocalNz;z++) {
            int localIndex = ROUND(index(mesh->xstart, y, z));
            ASSERT2( (localIndex >= 0) && (localIndex < localN) );
            if(z == 0) {
              // All 2D and 3D fields
              for(int i=0;i<n2d+n3d;i++) {
                d_nnz[localIndex+i] -= (n3d + n2d);
                o_nnz[localIndex+i] += (n3d + n2d);
              }
            }else {
              // Only 3D fields
              for(int i=0;i<n3d;i++) {
                d_nnz[localIndex+i] -= (n3d + n2d);
                o_nnz[localIndex+i] += (n3d + n2d);
              }
            }
          }
        }
      }

      if(mesh->lastX()) {
        // Upper X boundary
        for(int y=mesh->ystart;y<=mesh->yend;y++) {
          for(int z=0;z<mesh->LocalNz;z++) {
            int localIndex = ROUND(index(mesh->xend, y, z));
            ASSERT2( (localIndex >= 0) && (localIndex < localN) );
            if(z == 0) {
              // All 2D and 3D fields
              for(int i=0;i<n2d+n3d;i++)
                d_nnz[localIndex + i] -= (n3d + n2d);
            }else {
              // Only 3D fields
              for(int i=0;i<n3d;i++)
                d_nnz[localIndex + i] -= (n3d + n2d);
            }
          }
        }
      }else {
        // On another processor
        for(int y=mesh->ystart;y<=mesh->yend;y++) {
          for(int z=0;z<mesh->LocalNz;z++) {
            int localIndex = ROUND(index(mesh->xend, y, z));
            ASSERT2( (localIndex >= 0) && (localIndex < localN) );
            if(z == 0) {
              // All 2D and 3D fields
              for(int i=0;i<n2d+n3d;i++) {
                d_nnz[localIndex+i] -= (n3d + n2d);
                o_nnz[localIndex+i] += (n3d + n2d);
              }
            }else {
              // Only 3D fields
              for(int i=0;i<n3d;i++) {
                d_nnz[localIndex+i] -= (n3d + n2d);
                o_nnz[localIndex+i] += (n3d + n2d);
              }
            }
          }
        }
      }
      
      // Y boundaries
  
      for(int x=mesh->xstart; x <=mesh->xend; x++) {
        // Default to no boundary
        // NOTE: This assumes that communications in Y are to other
        //   processors. If Y is communicated with this processor (e.g. NYPE=1)
        //   then this will result in PETSc warnings about out of range allocations

        // z = 0 case
        int localIndex = ROUND(index(x, mesh->ystart, 0));
        // All 2D and 3D fields
        for(int i=0;i<n2d+n3d;i++) {
          //d_nnz[localIndex+i] -= (n3d + n2d);
          o_nnz[localIndex+i] += (n3d + n2d);
        }
        
        for(int z=1;z<mesh->LocalNz;z++) {
          localIndex = ROUND(index(x, mesh->ystart, z));
          
          // Only 3D fields
          for(int i=0;i<n3d;i++) {
            //d_nnz[localIndex+i] -= (n3d + n2d);
            o_nnz[localIndex+i] += (n3d + n2d);
          }
        }

        // z = 0 case
        localIndex = ROUND(index(x, mesh->yend, 0));
        // All 2D and 3D fields
        for(int i=0;i<n2d+n3d;i++) {
          //d_nnz[localIndex+i] -= (n3d + n2d);
          o_nnz[localIndex+i] += (n3d + n2d);
        }
        
        for(int z=1;z<mesh->LocalNz;z++) {
          localIndex = ROUND(index(x, mesh->yend, z));
          
          // Only 3D fields
          for(int i=0;i<n3d;i++) {
            //d_nnz[localIndex+i] -= (n3d + n2d);
            o_nnz[localIndex+i] += (n3d + n2d);
          }
        }
      }

      for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
        // A boundary, so no communication

        // z = 0 case
        int localIndex = ROUND(index(it.ind, mesh->ystart, 0));
        // All 2D and 3D fields
        for(int i=0;i<n2d+n3d;i++) {
          o_nnz[localIndex+i] -= (n3d + n2d);
        }
        
        for(int z=1;z<mesh->LocalNz;z++) {
          int localIndex = ROUND(index(it.ind, mesh->ystart, z));
          
          // Only 3D fields
          for(int i=0;i<n3d;i++) {
            o_nnz[localIndex+i] -= (n3d + n2d);
          }
        }
      }

      for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
        // A boundary, so no communication

        // z = 0 case
        int localIndex = ROUND(index(it.ind, mesh->yend, 0));
        // All 2D and 3D fields
        for(int i=0;i<n2d+n3d;i++) {
          o_nnz[localIndex+i] -= (n3d + n2d);
        }
        
        for(int z=1;z<mesh->LocalNz;z++) {
          int localIndex = ROUND(index(it.ind, mesh->yend, z));
          
          // Only 3D fields
          for(int i=0;i<n3d;i++) {
            o_nnz[localIndex+i] -= (n3d + n2d);
          }
        }
      }
      
      // Pre-allocate
      MatMPIAIJSetPreallocation( Jmf, 0, d_nnz, 0, o_nnz );
      MatSetUp(Jmf); 
      MatSetOption(Jmf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);      
      PetscFree( d_nnz );
      PetscFree( o_nnz );
      
      // Determine which row/columns of the matrix are locally owned
      int Istart, Iend;
      MatGetOwnershipRange( Jmf, &Istart, &Iend );
      
      // Convert local into global indices
      index += Istart;
      
      // Now communicate to fill guard cells
      mesh->communicate(index);

      //////////////////////////////////////////////////
      // Mark non-zero entries

      
      // Offsets for a 5-point pattern
      const int xoffset[5] = {0,-1, 1, 0, 0};
      const int yoffset[5] = {0, 0, 0,-1, 1};
      
      PetscScalar val = 1.0;
      
      for(int x=mesh->xstart; x <= mesh->xend; x++) {
        for(int y=mesh->ystart;y<=mesh->yend;y++) {
          
          int ind0 = ROUND(index(x,y,0));

          // 2D fields
          for(int i=0;i<n2d;i++) {
            PetscInt row = ind0 + i;

            // Loop through each point in the 5-point stencil
            for(int c=0;c<5;c++) {
              int xi = x + xoffset[c];
              int yi = y + yoffset[c];
                
              if( (xi < 0) || (yi < 0) ||
                  (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy) )
                continue;
              
              int ind2 = ROUND(index(xi, yi, 0));
              
              if(ind2 < 0)
                continue; // A boundary point
              
              // Depends on all variables on this cell
              for(int j=0;j<n2d;j++) {
                PetscInt col = ind2 + j;

                //output.write("SETTING 1: {:d}, {:d}\n", row, col);
                MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
              }
            }
          }
          
          // 3D fields
          for(int z=0;z<mesh->LocalNz;z++) {
            
            int ind = ROUND(index(x,y,z));
            
            for(int i=0;i<n3d;i++) {
              PetscInt row = ind + i;
              if(z == 0)
                row += n2d;
              
              // Depends on 2D fields
              for(int j=0;j<n2d;j++) {
                PetscInt col = ind0 + j;
                //output.write("SETTING 2: {:d}, {:d}\n", row, col);
                MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
              }
              
              // 5 point star pattern
              for(int c=0;c<5;c++) {
                int xi = x + xoffset[c];
                int yi = y + yoffset[c];
                
                if( (xi < 0) || (yi < 0) ||
                    (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy) )
                  continue;
                
                int ind2 = ROUND(index(xi, yi, z));
                if(ind2 < 0)
                  continue; // Boundary point
                
                if(z == 0)
                  ind2 += n2d;
                
                // 3D fields on this cell
                for(int j=0;j<n3d;j++) {
                  PetscInt col = ind2 + j;
                  //output.write("SETTING 3: {:d}, {:d}\n", row, col);
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }
              }

              int nz = mesh->LocalNz;
              if(nz > 1) {
                // Multiple points in z
                
                int zp = (z + 1) % nz;

                int ind2 = ROUND(index(x, y, zp));
                if(zp == 0)
                  ind2 += n2d;
                for(int j=0;j<n3d;j++) {
                  PetscInt col = ind2 + j;
                  //output.write("SETTING 4: {:d}, {:d}\n", row, col);
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }

                int zm = (z - 1 + nz) % nz;
                ind2 = ROUND(index(x, y, zm));
                if(zm == 0)
                  ind2 += n2d;
                for(int j=0;j<n3d;j++) {
                  PetscInt col = ind2 + j;
                  //output.write("SETTING 5: {:d}, {:d}\n", row, col);
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }
                
              }
              
            }
          }
        }
      }
      // Finished marking non-zero entries
      
      // Assemble Matrix
      MatAssemblyBegin( Jmf, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( Jmf, MAT_FINAL_ASSEMBLY );
      
      
      ISColoring iscoloring;
      
#if PETSC_VERSION_GE(3,5,0)
      MatColoring coloring; // This new in PETSc 3.5
      MatColoringCreate(Jmf,&coloring);
      MatColoringSetType(coloring,MATCOLORINGSL);
      MatColoringSetFromOptions(coloring);
      // Calculate index sets
      MatColoringApply(coloring,&iscoloring);
      MatColoringDestroy(&coloring);
#else
      // Pre-3.5
      MatGetColoring(Jmf,MATCOLORINGSL,&iscoloring);
#endif

      // Create data structure for SNESComputeJacobianDefaultColor
      MatFDColoringCreate(Jmf,iscoloring,&fdcoloring);
      ISColoringDestroy(&iscoloring);
      // Set the function to difference
      //MatFDColoringSetFunction(fdcoloring,(PetscErrorCode (*)(void))FormFunctionForDifferencing,this);
      MatFDColoringSetFunction(
          fdcoloring, reinterpret_cast<PetscErrorCode(*)()>(FormFunctionForColoring),
          this);
      MatFDColoringSetFromOptions(fdcoloring);
      //MatFDColoringSetUp(Jmf,iscoloring,fdcoloring);
      
#if PETSC_VERSION_GE(3,4,0)
      SNESSetJacobian(*snesIn,Jmf,Jmf,SNESComputeJacobianDefault,fdcoloring);
#else
      // Before 3.4
      SNESSetJacobian(*snesIn,Jmf,Jmf,SNESDefaultComputeJacobian,fdcoloring);
#endif

      // Re-use Jacobian
      int lag_jacobian;
      OPTION(options, lag_jacobian,   4);
      SNESSetLagJacobian(*snesIn,lag_jacobian);
      
      //MatView(Jmf, PETSC_VIEWER_DRAW_WORLD);
      //MatView(Jmf, PETSC_VIEWER_STDOUT_WORLD);
    }else {
      // Brute force calculation
      // NOTE: Slow!
      
      MatCreateAIJ(BoutComm::get(),
                 nlocal,nlocal,  // Local sizes
                 PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
                 3,   // Number of nonzero entries in diagonal portion of local submatrix
                 PETSC_NULL,
                 0,   // Number of nonzeros per row in off-diagonal portion of local submatrix
                 PETSC_NULL,
                 &Jmf);
      
#if PETSC_VERSION_GE(3,4,0)
    SNESSetJacobian(*snesIn,Jmf,Jmf,SNESComputeJacobianDefault,this);
#else
    // Before 3.4
    SNESSetJacobian(*snesIn,Jmf,Jmf,SNESDefaultComputeJacobian,this); 
#endif
    
      MatSetOption(Jmf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    }
  }
  
  /////////////////////////////////////////////////////
  // Set tolerances
  BoutReal atol, rtol; // Tolerances for SNES solver
  options->get("atol", atol, 1e-16);
  options->get("rtol", rtol, 1e-10);
  int max_nonlinear_it = (*options)["max_nonlinear_iterations"]
                         .doc("Maximum number of nonlinear iterations per SNES solve")
                         .withDefault(5);
  SNESSetTolerances(*snesIn,atol,rtol,PETSC_DEFAULT,max_nonlinear_it,PETSC_DEFAULT);

  /////////////////////////////////////////////////////
  // Predictor method
  options->get("predictor", predictor, 1);

  /////////////////////////////////////////////////////
  // Preconditioner

  bool use_precon;
  OPTION(options, use_precon,   false);
 
  // Get KSP context from SNES
  KSP ksp;
  SNESGetKSP(*snesIn, &ksp);
  
  bool kspsetinitialguessnonzero;
  options->get("kspsetinitialguessnonzero", kspsetinitialguessnonzero, false);
  if(kspsetinitialguessnonzero) {
    //Set the initial guess to be non-zero
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  }

  int maxl; // Maximum number of linear iterations
  OPTION(options, maxl, 20);
  KSPSetTolerances(ksp, 
                   PETSC_DEFAULT,  // rtol
                   PETSC_DEFAULT,  // abstol
                   PETSC_DEFAULT,  // dtol (divergence tolerance)
                   maxl);  // Maximum number of iterations

  // Get PC context from KSP
  PC pc;
  KSPGetPC(ksp,&pc);

  if (use_precon && hasPreconditioner()) {
    output.write("\tUsing user-supplied preconditioner\n");

    // Set a Shell (matrix-free) preconditioner type
    PCSetType(pc, PCSHELL);

    // Specify the preconditioner function
    PCShellSetApply(pc,imexbdf2PCapply);
    // Context used to supply object pointer
    PCShellSetContext(pc,this);
  } else if (matrix_free) {
    PCSetType(pc, PCNONE);
  }

  /////////////////////////////////////////////////////
  // diagnostics
  
  OPTION(options, diagnose, false); // Print diagnostics
  OPTION(options, verbose, false); // More outputs at each timestep
  
  /////////////////////////////////////////////////////
  // Get runtime options
  SNESSetFromOptions(*snesIn);

  if(diagnose) {
    //Some reporting
    PCType pctype; PCGetType(pc, &pctype);
    KSPType ksptype; KSPGetType(ksp, &ksptype);
    SNESType snestype; SNESGetType(*snesIn, &snestype);
    output<<"SNES Type : "<<snestype<<endl;
    output<<"KSP Type : "<<ksptype<<endl;
    output<<"PC Type : "<<pctype<<endl;
  }

}

int IMEXBDF2::run() {
  TRACE("IMEXBDF2::run()");

  // Multi-step scheme, so first steps are different
  int order = 1;
  int lastOrder = -1;
  BoutReal dt = timestep;
  std::vector<BoutReal> lastTimesteps = timesteps;
  BoutReal dtNext = dt; //Timestep to try for next internal iteration
  
  //By default use the main snes object.
  snesUse = snes;

  int internalCounter=0; //Cumulative number of successful internal iterations

  for(int s=0;s<nsteps;s++) {
    BoutReal cumulativeTime = 0.;
    int counter = 0; //How many iterations in this output step
    
    // Reset linear and nonlinear fail counts
    linear_fails = 0;
    nonlinear_fails = 0;
    while(cumulativeTime<out_timestep){
      //Move state history along one stage (i.e. u_2-->u_3,u_1-->u_2, u-->u_1 etc.)
      //Note: This sets the current timestep to be the same as the last timestep.
      shuffleState();

      //First part of time step -- Run the convective part to find f_1
      // Calculate time-derivative of u_1, put into f_1
      /*
	Would be more efficient to include this loadVars at the end of the inner loop to avoid
	need to loadVars(u) before doing the monitors. Would need a setup call outside
	loops however and savings probably minimal.
      */
      loadVars(std::begin(uV[0]));
      run_convective(simtime);
      saveDerivs(std::begin(fV[0]));

      bool running = true;
      bool checkingErr = adaptive && (internalCounter%nadapt) ==0 && order>1;
      int adaptCounter=0;
      int failCounter = 0; // Number of failed steps
      while(running){
	running = false;

	//Validate our desired next timestep
	if(dtNext<dtMinFatal){ //Don't allow the timestep to go below requested fatal min
	  throw BoutException("Aborting: Timestep ({:f}) tried to go below minimum allowed",dtNext);
	}else if(dtNext<dtMin){ //Don't allow timestep below requested min
	  dtNext = dtMin;
	}else if(dtNext>dtMax){ //Don't allow timestep above request max
	  dtNext = dtMax;
	}else{ //Timestep is fine so don't do anything
	};
	
	//Check if we will go past the target time (i.e. past the output step).
	//If so we want to limit the timestep.
	//There's potential for this to confuse the adaptive calculation so
	//we'll set a flag to alert us to this forced change. 
	bool artificalLimit = false;
        BoutReal dtNoLimit = dtNext; // What dt would have been without artificial limit
	if(cumulativeTime+dtNext > out_timestep){
	  artificalLimit = true;
	  dtNext = out_timestep - cumulativeTime;
	}

        if(verbose) {
          output << endl << "At t=" << cumulativeTime << " attempting internal step "<<counter<<" (attempt "<<adaptCounter<<")"<<endl;
          output << "Using dt = "<<dtNext<<endl;
        }

	//Set the current timestep to try -- Has to be before calculateCoeffs call
	timesteps[0] = dtNext;

	//If we're checking the error at this point do the low order solution now
	if(checkingErr){
	  //First find coefficients for use with lower order scheme
	  calculateCoeffs(order-1);  

	  //Use alternative snes solver for low order scheme
	  snesUse = snesAlt;

	  //Solve
          try {
            take_step(simtime, timesteps[0], order-1);
          }catch (const BoutException &e) {
            // An error occurred. If adaptive, reduce timestep
            if(!adaptive)
              throw;

            failCounter++;
            if(failCounter > 10) {
              throw BoutException("Too many failed steps\n");
            }
            
            running = true; // Retry
            dtNext = timesteps[0]*0.5;
            continue;
          }
	  
	  //Store this solution in err
	  for(int i=0;i<nlocal;i++){
	    err[i] = u[i];
	  };

	  //Go back to using the main snes object
	  snesUse = snes;
	}

	//Now we get the coefficients if the order has changed *or* any of the
	//timesteps in the history are different *or* we had to calc coefficients
	//for lower order scheme.
	if( (order!=lastOrder) || (lastTimesteps != timesteps) || checkingErr){
	  calculateCoeffs(order);  
	}

	//Now we complete the timestep by constructing rhs and solving the implicit part
        try {
          take_step(simtime, timesteps[0], order);
        }catch (const BoutException &e) {
          // An error occurred. If adaptive, reduce timestep
          if(!adaptive)
            throw;

          failCounter++;
          if(failCounter > 10) {
            throw BoutException("Too many failed steps\n");
          }

          running = true; // Retry
          dtNext = timesteps[0]*0.5;
          continue;
        }

	//Now we can calculate the error and decide what we want to do
	if(checkingErr){
	  //Now we want to find the actual (abs) error
	  BoutReal errTot[3] = {0,0,0}; 
	  BoutReal errGlobTot[3] = {0,0,0};

	  //Find local data
	  for(int i=0;i<nlocal;i++){
	    errTot[0] += std::abs(err[i]-u[i]);
	    errTot[1] += std::abs(u[i]);
	    errTot[2] += std::abs(err[i]);
	  };

	  //Now reduce across procs
          bout::globals::mpi->MPI_Allreduce(&errTot, &errGlobTot, 3, MPI_DOUBLE, MPI_SUM,
                                            BoutComm::get());

          BoutReal aRtol = errGlobTot[0]/errGlobTot[1];
          if(verbose) {
            output<<"The average errors are aerr = "<<errGlobTot[0]<<" and rerr = "<<aRtol<<endl;
            output<<"The err mag is "<<errGlobTot[2]<<" and the sol mag is "<<errGlobTot[1]<<endl;
          }

	  /*
	   * The following is how we argue the timestep should be scaled (s) 
	   * to achieve the target error (adaptRtol).
	   * U_{N-1} = T + O(dt^{N-1})
	   * U_{N} = T + O(dt^N)
	   * Next we assume U_{N} == T
	   * => aRtol ~ C dt^{N-1}   {A}
	   * adaptRtol = C (s*dt)^{N-1} = delta*aRtol {B}
	   * {B}/{A} = delta = s^{N-1}
	   * => s is the {N-1}th root of delta, where delta = adaptRtol/aRtol
	   * If s<scaleCushDown then we recommend a reduction in the timestep.
	   * If scaleCushDown<=s<scaleCushUp timestep change possible but not worth it
	   * If s>=scaleCushUp we want to increase the timestep next time
	   * It should be noted that this may not be a good approximation, particularly
	   * when we use a low order scheme. In addition we argue U_N = T but then assume
	   * aRtol represents the error on the U_N as this is what we follow.
	   */
	  BoutReal delta = adaptRtol/aRtol;
	  BoutReal s = pow(delta, 1.0/(order-1.0));

	  //Work out if we need to change the timestep and repeat this step
	  if(s<scaleCushDown){
	    running = true;
	    dtNext = timesteps[0]*s; 
	  }else if( (s>=scaleCushUp) && (adaptCounter==0) && (failCounter == 0) ){ 
	    // Here we decide to increase the timestep
	    // but note we only allow this if this is the first attempt at this step.
            // and if there have been no failed steps
	    // This is designed to prevent oscillation in timestep.
            
            s = BOUTMIN(s, 1.25); // Limit increase 
	    dtNext = timesteps[0]*s;
	  }else{ //No change to the timestep
	    dtNext = timesteps[0];
	  }

          if(verbose) {
            output << "Error ratio is "<<delta<<" so scaling factor is "<<s<<" and dtNext is "<<dtNext<<endl;
          }
	
	  adaptCounter++;
	  if(adaptCounter>mxstepAdapt){
	    throw BoutException("Aborting: Maximum number of adapative iterations ({:d}) exceeded", mxstepAdapt);
	  }
	}else if(artificalLimit) {
          // Reset dtNext if it was artificially limited
          // to the value it would have been without artificial limit
          dtNext = dtNoLimit;
        }
      }//End of running -- Done a single internal step

      //Update record of what was used to complete this step
      lastOrder = order;
      lastTimesteps = timesteps;

      //Increment order if we're not at the maximum requested
      if(order<maxOrder) order++;

      //Update simulation time and record of how far through this output step we are.
      simtime += timesteps[0];
      cumulativeTime += timesteps[0];

      call_timestep_monitors(simtime, timesteps[0]);
      
      //Increment internal counter to keep track of number of internal iterations
      internalCounter++;

      //Increment iteration counter to ensure we don't get an infinite loop
      counter++;
      if(counter>mxstep){
	throw BoutException("Aborting: Maximum number of internal iterations ({:d}) exceeded", mxstep);
      };
    }

    if(diagnose) {
      output.write("\n   Last dt = {:e}, order = {:d}\n", timesteps[0], lastOrder);
      output.write("   Linear fails = {:d}, nonlinear fails = {:d}\n", linear_fails, nonlinear_fails);
    }

    loadVars(std::begin(u));// Put result into variables
    run_rhs(simtime); // Run RHS to calculate auxilliary variables

    iteration++; // Advance iteration number

    /// Call the monitor function

    if(call_monitors(simtime, s, nsteps)) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

/*
 * Calculate the coefficients required for this order calculation
 * See: http://summit.sfu.ca/item/9862 for more details
 * Note that reference uses increasing index to indicate increasing time
 * i.e. u_n is one stage older than u_(n+1). This is opposite to the order 
 * we use here where uV[0] is newer than uV[1].
 * Note: uFac corresponds to -alpha. fFrac corresponds to beta and gFac corresponds to C
 * The BDF schemes have gamma=1, theta,c = 0
 * Currently only implemented schemes with gFac=0 (these are BDF schemes) could look at
 * changing this, for order>2 TVB schemes should be better, for details see
 * http://homepages.cwi.nl/~willem/DOCART/JCP07.pdf
 */
void IMEXBDF2::calculateCoeffs(int order){
  BoutReal uCurrFac;

  switch (order){
  case 1: {
    uCurrFac = 1.0;
    uFac[0] = 1.0;
    fFac[0] = timesteps[0];
    dtImp = timesteps[0];
    break;
  }
  case 2: {
    BoutReal omega1 = timesteps[0]/timesteps[1];
    uCurrFac = (1+2*omega1)/(1+omega1);
    uFac[0] = (1+omega1);
    uFac[1] = -pow(omega1,2)/(1+omega1);
    fFac[0] = timesteps[0]*(1+omega1);
    fFac[1] = -timesteps[0]*omega1;
    dtImp = timesteps[0];
    break;
  }
  case 3: {
    BoutReal omega1 = timesteps[1]/timesteps[2];
    BoutReal omega2 = timesteps[0]/timesteps[1];
    uCurrFac = 1 + omega2/(1+omega2) + omega1*omega2/(1+omega1*(1+omega2));
    uFac[0] = 1 + omega2 + omega1*omega2*(1+omega2)/(1+omega1);
    uFac[1] = -pow(omega2,2)*(omega1+1/(1+omega2));
    uFac[2] = pow(omega1,3)*pow(omega2,2)*(1+omega2)/((1+omega1)*(1+omega1+omega1*omega2));
    fFac[0] = timesteps[0]*(1+omega2)*(1+omega1*(1+omega2))/(1+omega1);
    fFac[1] = -timesteps[0]*omega2*(1+omega1*(1+omega2));
    fFac[2] = timesteps[0]*pow(omega1,2)*omega2*(1+omega2)/(1+omega1);
    dtImp = timesteps[0];
    break;
  }
  case 4: {
    BoutReal omega1 = timesteps[2]/timesteps[3];
    BoutReal omega2 = timesteps[1]/timesteps[2];
    BoutReal omega3 = timesteps[0]/timesteps[1];
    BoutReal A1 = 1+omega1*(1+omega2);
    BoutReal A2 = 1+omega2*(1+omega3);
    BoutReal A3 = 1+omega1*A2;
    uCurrFac = 1 + omega3/(1+omega3) + omega2*omega3/A2 + omega1*omega2*omega3/A3;
    uFac[0] = 1 + omega3*(1+omega2*(1+omega3)*(1+omega1*A2/A1)/(1+omega2));
    uFac[1] = -omega3*(omega3/(1+omega3) + omega2*omega3*(A3+omega1)/(1+omega1));
    uFac[2] = pow(omega2,3)*pow(omega3,2)*(1+omega3)*A3/((1+omega2)*A2);
    uFac[3] = -((1+omega3)/(1+omega1))*(A2/A1)*pow(omega1,4)*pow(omega2,3)*pow(omega3,2)/A3;
    fFac[0] =  timesteps[0]*
      (omega2*(1+omega3)/(1+omega2))*
      ((1+omega3)*(A3+omega1)+(1+omega1)/omega2)/A1;
    fFac[1] = -timesteps[0]*A2*A3*omega3/(1+omega1);
    fFac[2] =  timesteps[0]*pow(omega2,2)*omega3*A3*(1+omega3)/(1+omega2);
    fFac[3] = -timesteps[0]*pow(omega1,3)*pow(omega2,2)*omega3*(A2/A1)*(1+omega3)/(1+omega1);
    dtImp = timesteps[0];
    break;
  }
  default:{
    throw BoutException("Invalid order supplied in IMEXBDF2::calculateCoeffs");
  }
  };

  //Scale the factors by uCurrFac
  for(int i=0;i<order;i++){
    uFac[i] /= uCurrFac;
    fFac[i] /= uCurrFac;
    gFac[i] /= uCurrFac;
  }
  dtImp /= uCurrFac;
}

/*!
 * Take a full IMEX-BDF step of order "order". Note that this assumes
 * that enough time points are already available (in u and f).
 *
 * Inputs:
 * u*   - Solution history
 * f*   - Non-stiff component history
 *
 * Outputs:
 * u   - Latest Solution
 * f1  - Non-stiff time derivative at current time
 */
void IMEXBDF2::take_step(BoutReal curtime, BoutReal UNUSED(dt), int order) {

  //First zero out rhs
  std::fill(std::begin(rhs), std::end(rhs), 0.0);

  //Now add the contribution to rhs from each history step
  for(int j=0;j<order;j++){
    for(int i=0;i<nlocal;i++){
      rhs[i] += uV[j][i]*uFac[j] + fV[j][i]*fFac[j];
    }
  }

  // Now need to solve u - dtImp*G(u) = rhs
  solve_implicit(curtime+timesteps[0], dtImp);
}

/* 
 *  Moves the solution histories along one step. Also handles the timesteps.
 */
void IMEXBDF2::shuffleState(){

  //Note: std::rotate takes the start and end of a range and a third value (2nd arg)
  //which says rotate the elements of the vector such that this element is first.
  //Here we want the last point in history to become the work array for the first
  //as we are losing this last point from our records.

  //Shuffle stashed values along a step
  //Non-stiff solutions
  std::rotate(fV.begin(),fV.end()-1,fV.end());

  // Rotate u -> u_1, u_1 -> u_2, u_2 -> u . U later overwritten
  std::rotate(uV.begin(),uV.end()-1,uV.end()); //Rotate
  //Slight extra handling required as the current state "u" is held externally 
  //from the history vector *for reasons*
  std::swap(u, uV[0]);

  //Timesteps used
  std::rotate(timesteps.begin(),timesteps.end()-1,timesteps.end());
  //Note -- timesteps[0] is currently not correct in general. Must set it
  //before we take_step. This is a bit unpleasent but is somewhat necessary
  //in order to allow us to trial different step sizes. Really it's probably
  //nicer if we at least leave timesteps in a sensible default state. Hence
  //here we say lets just use the same timestep as last time by default.
  //That way we only need to fiddle with timesteps if we're adapting.
  timesteps[0] = timesteps[1];
}

/*
 * Solves u - gamma*G(u) = rhs
 *
 * where u is the result, G(u) is the stiff part of the rhs (run_diffusive)
 * and gamma is a factor depending on the time-step and method
 *
 * Inputs:
 * rhsvec
 *
 *
 */
PetscErrorCode IMEXBDF2::solve_implicit(BoutReal curtime, BoutReal gamma) {
  implicit_curtime = curtime;
  implicit_gamma = gamma;

  // Set initial guess at the solution
  BoutReal *xdata;
  int ierr;
  ierr = VecGetArray(snes_x,&xdata);CHKERRQ(ierr);

  switch(predictor) {
  case 0: {
    // Constant, so next step is same as last step
    for(int i=0;i<nlocal;i++) {
      xdata[i] = uV[0][i];     // Use previous solution
    }
    break;
  }
  case 1: {
    // Linear extrapolation from last two steps
    for(int i=0;i<nlocal;i++) {
      xdata[i] = 2.*uV[0][i] - uV[1][i];
    }
    break;
  }
  case 2: {
    // Quadratic extrapolation.
    for(int i=0;i<nlocal;i++) {
      xdata[i] = 3.*uV[0][i] - 3.*uV[1][i] + uV[2][i];
    }
    break;
  }
  //Could add a cubic extrapolation here
  //
  default: {
    // Assume that there is no non-linear solve, so G = 0
    for(int i=0;i<nlocal;i++) {
      xdata[i] = rhs[i];   // If G = 0
    }
  }
  }

  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);

  SNESSolve(snesUse, nullptr, snes_x);

  // Find out if converged
  SNESConvergedReason reason;
  SNESGetConvergedReason(snesUse,&reason);
  if(reason < 0) {
    // Diverged
    KSP ksp;
    SNESGetKSP(snesUse, &ksp);
    KSPConvergedReason kreason;
    KSPGetConvergedReason(ksp,&kreason);
    if(kreason<0){
      if(verbose) {
        output<<"KSP Failed to converge with reason "<<kreason<<endl;
      }
      linear_fails++;
    }else{
      nonlinear_fails++;
      if(verbose) {
        output << "KSP Succeeded with reason "<<kreason<<endl;
      }
    };
    if(verbose) {
      output << "SNES failed to converge with reason " << reason << endl;
    }
    throw BoutException("SNES failed to converge. Reason: {:d}\n", reason);
  }

  int its;
  SNESGetIterationNumber(snesUse,&its);

  if(verbose) {
    output << "Number of SNES iterations: " << its << endl;
  }
    
  // Put the result into u
  ierr = VecGetArray(snes_x,&xdata);CHKERRQ(ierr);

  for(int i=0;i<nlocal;i++)
    u[i] = xdata[i];
  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);

  return 0;
}

// f = (x - gamma*G(x)) - rhs
PetscErrorCode IMEXBDF2::snes_function(Vec x, Vec f, bool linear) {
  const BoutReal *xdata;
  BoutReal *fdata;
  int ierr;

  // Get data from PETSc into BOUT++ fields
  ierr = VecGetArrayRead(x,&xdata);CHKERRQ(ierr);

  // Hacky
  loadVars(const_cast<BoutReal*>(xdata));

  // Call RHS function
  run_diffusive(implicit_curtime, linear);

  // Copy derivatives back
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  saveDerivs(fdata);

  // G(x) now in fdata
  
  if(!have_constraints) {
    // No constraints, so simple loop over all variables
    
    for(int i=0;i<nlocal;i++) {
      fdata[i] = xdata[i] - implicit_gamma * fdata[i] - rhs[i];
    }
  }else {
    // Some constraints
    for(int i=0;i<nlocal;i++) {
      if(is_dae[i] > 0.5) { // 1 -> differential, 0 -> algebraic
        fdata[i] = xdata[i] - implicit_gamma * fdata[i] - rhs[i];
      }
      // Otherwise the constraint is that fdata[i] = 0.0
    }
  }

  // Restore data arrays to PETSc
  ierr = VecRestoreArray(f,&fdata);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(x,&xdata);CHKERRQ(ierr);

  return 0;
}

/*
 * Preconditioner function
 */
PetscErrorCode IMEXBDF2::precon(Vec x, Vec f) {
  if (!hasPreconditioner()) {
    // No user preconditioner
    throw BoutException("No user preconditioner");
  }

  int ierr;

  // Get data from PETSc into BOUT++ fields
  Vec solution;
  SNESGetSolution(snes, &solution);
  BoutReal *soldata;
  ierr = VecGetArray(x,&soldata);CHKERRQ(ierr);
  load_vars(soldata);
  ierr = VecRestoreArray(solution,&soldata);CHKERRQ(ierr);

  // Load vector to be inverted into ddt() variables
  BoutReal *xdata;
  ierr = VecGetArray(x,&xdata);CHKERRQ(ierr);
  load_derivs(xdata);
  ierr = VecRestoreArray(x,&xdata);CHKERRQ(ierr);

  // Run the preconditioner
  runPreconditioner(implicit_curtime, implicit_gamma, 0.0);

  // Save the solution from F_vars
  BoutReal *fdata;
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  save_derivs(fdata);
  ierr = VecRestoreArray(f,&fdata);CHKERRQ(ierr);

  return 0;
}

/*!
 * Loop over arrays, using template parameter
 * to specify the operation to be performed at each point
 *
 */
template< class Op >
void IMEXBDF2::loopVars(BoutReal *u) {
  // Use global mesh for now
  Mesh* mesh = bout::globals::mesh;

  // Loop over 2D variables
  for(auto & it : f2d) {
    Op op(it.var, it.F_var); // Initialise the operator

    if(it.evolve_bndry) {
      // Include boundary regions

      // Inner X
      if(mesh->firstX() && !mesh->periodicX) {
        for(int jx=0;jx<mesh->xstart;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy) {
            op.run(jx, jy, u); ++u;
          }
      }

      // Outer X
      if(mesh->lastX() && !mesh->periodicX) {
        for(int jx=mesh->xend+1;jx<mesh->LocalNx;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy) {
            op.run(jx, jy, u); ++u;
          }
      }
      // Lower Y
      for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); ++xi) {
        for(int jy=0;jy<mesh->ystart;++jy) {
          op.run(*xi, jy, u); ++u;
        }
      }

      // Upper Y
      for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); ++xi) {
        for(int jy=mesh->yend+1;jy<mesh->LocalNy;++jy) {
          op.run(*xi, jy, u); ++u;
        }
      }
    }

    // Bulk of points
    for(int jx=mesh->xstart; jx <= mesh->xend; ++jx)
      for(int jy=mesh->ystart; jy <= mesh->yend; ++jy) {
        op.run(jx, jy, u); ++u;
      }
  }

  // Loop over 3D variables
  for(auto & it : f3d) {
    Op op(it.var, it.F_var); // Initialise the operator
    if(it.evolve_bndry) {
      // Include boundary regions

      // Inner X
      if(mesh->firstX() && !mesh->periodicX) {
        for(int jx=0;jx<mesh->xstart;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy)
            for(int jz=0; jz < mesh->LocalNz; ++jz) {
              op.run(jx, jy, jz, u); ++u;
            }
      }

      // Outer X
      if(mesh->lastX() && !mesh->periodicX) {
        for(int jx=mesh->xend+1;jx<mesh->LocalNx;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy)
            for(int jz=0; jz < mesh->LocalNz; ++jz) {
              op.run(jx, jy, jz, u); ++u;
            }
      }
      // Lower Y
      for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); ++xi) {
        for(int jy=0;jy<mesh->ystart;++jy)
          for(int jz=0; jz < mesh->LocalNz; ++jz) {
            op.run(*xi, jy, jz, u); ++u;
          }
      }

      // Upper Y
      for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); ++xi) {
        for(int jy=mesh->yend+1;jy<mesh->LocalNy;++jy)
          for(int jz=0; jz < mesh->LocalNz; ++jz) {
            op.run(*xi, jy, jz, u); ++u;
          }
      }
    }

    // Bulk of points
    for(int jx=mesh->xstart; jx <= mesh->xend; ++jx)
      for(int jy=mesh->ystart; jy <= mesh->yend; ++jy)
        for(int jz=0; jz < mesh->LocalNz; ++jz) {
          op.run(jx, jy, jz, u); ++u;
        }
  }
}

///////////////////////////////////////////////////////////////////

class SaveVarOp {
public:
  // Initialise with a Field2D iterator
  SaveVarOp(Field2D *var, Field2D *UNUSED(F_var)) : var2D(var) {}
  // Initialise with a Field3D iterator
  SaveVarOp(Field3D *var, Field3D *UNUSED(F_var)) : var3D(var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    *u = (*var2D)(jx,jy);
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    *u = (*var3D)(jx,jy,jz);
  }
private:
  Field2D *var2D;
  Field3D *var3D;
};

/*!
 * Copy data from fields into array
 */
void IMEXBDF2::saveVars(BoutReal *u) {
  //loopVars<SaveVarOp>(u);
  save_vars(u);
}

///////////////////////////////////////////////////////////////////

class LoadVarOp {
public:
  // Initialise with a Field2D iterator
  LoadVarOp(Field2D *var, Field2D *UNUSED(F_var)) : var2D(var) {}
  // Initialise with a Field3D iterator
  LoadVarOp(Field3D *var, Field3D *UNUSED(F_var)) : var3D(var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    (*var2D)(jx,jy) = *u;
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    (*var3D)(jx,jy,jz) = *u;
  }
private:
  Field2D *var2D;
  Field3D *var3D;
};

/*!
 * Copy data from array into fields
 */
void IMEXBDF2::loadVars(BoutReal *u) {
  //loopVars<LoadVarOp>(u);
  load_vars(u);
}

///////////////////////////////////////////////////////////////////

class SaveDerivsOp {
public:
  // Initialise with a Field2D iterator
  SaveDerivsOp(Field2D *UNUSED(var), Field2D *F_var) : F_var2D(F_var) {}
  // Initialise with a Field3D iterator
  SaveDerivsOp(Field3D *UNUSED(var), Field3D *F_var) : F_var3D(F_var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    *u = (*F_var2D)(jx,jy);
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    *u = (*F_var3D)(jx,jy,jz);
  }
private:
  Field2D *F_var2D;
  Field3D *F_var3D;
};

/*!
 * Copy time derivatives from fields into array
 */
void IMEXBDF2::saveDerivs(BoutReal *u) {
  //loopVars<SaveDerivsOp>(u);
  save_derivs(u);
}

#endif // BOUT_HAS_PETSC
