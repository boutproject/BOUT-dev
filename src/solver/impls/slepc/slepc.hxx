/**************************************************************************
 * Interface to SLEPc solver
 * NOTE: This class needs tidying, generalising to use FieldData interface
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson, D. Dickinson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __SLEPC_SOLVER_H__
#define __SLEPC_SOLVER_H__

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_SLEPC

namespace {
RegisterUnavailableSolver registerunavailableslepc("slepc",
                                                   "BOUT++ was not configured with SLEPc");
}

#else

class SlepcSolver;

#include <slepc.h>
// PETSc creates macros for MPI calls, which interfere with the MpiWrapper class
#undef MPI_Allreduce
#undef MPI_Gatherv
#undef MPI_Irecv
#undef MPI_Isend
#undef MPI_Recv
#undef MPI_Scatterv
#undef MPI_Send
#undef MPI_Wait
#undef MPI_Waitall
#undef MPI_Waitany

#include <field2d.hxx>
#include <field3d.hxx>
#include <utils.hxx>
#include <vector2d.hxx>
#include <vector3d.hxx>

#include <bout/petsclib.hxx>
#include <bout/slepclib.hxx>
#include <vector>

#define OPT_SIZE 40

// Define a name to use with SolverType to indicate SlepcSolver
// is in charge of advancing fields
#define SOLVERSLEPCSELF "self"

namespace {
RegisterSolver<SlepcSolver> registersolverslepc("slepc");
}

class SlepcSolver : public Solver {
public:
  SlepcSolver(Options *options);
  ~SlepcSolver();

  int advanceStep(Mat &matOperator, Vec &inData, Vec &outData);
  int compareEigs(PetscScalar ar, PetscScalar ai, PetscScalar br, PetscScalar bi);
  void monitor(PetscInt its, PetscInt nconv, PetscScalar eigr[], PetscScalar eigi[],
               PetscReal errest[], PetscInt nest);

  // These contain slepc specific code and call the advanceSolver code
  int init(int NOUT, BoutReal TIMESTEP) override;
  int run() override;

  ////////////////////////////////////////
  /// OVERRIDE
  ///      Here we override *all* other
  ///      virtual functions in order to
  ///      pass through control to the
  ///      actual solver (advanceSolver)
  ///      This is only required if allow
  ///      use of additional solver
  ////////////////////////////////////////

  void setModel(PhysicsModel *model) override { // New API
    Solver::setModel(model);
    if (!selfSolve) {
      advanceSolver->setModel(model);
    }
  }

  void setRHS(rhsfunc f) override { // Old API
    Solver::setRHS(f);
    if (!selfSolve) {
      advanceSolver->setRHS(f);
    }
  }

  //////Following overrides all just pass through to advanceSolver

  // Override virtual add functions in order to pass through to advanceSolver
  void add(Field2D& v, const std::string& name) override {
    Solver::add(v, name);
    if (!selfSolve) {
      advanceSolver->add(v, name);
    }
  }
  void add(Field3D& v, const std::string& name) override {
    Solver::add(v, name);
    if (!selfSolve) {
      advanceSolver->add(v, name);
    }
  }
  void add(Vector2D& v, const std::string& name) override {
    Solver::add(v, name);
    if (!selfSolve) {
      advanceSolver->add(v, name);
    }
  }
  void add(Vector3D& v, const std::string& name) override {
    Solver::add(v, name);
    if (!selfSolve) {
      advanceSolver->add(v, name);
    }
  }

  // Set operations
  void setJacobian(Jacobian j) override {
    if (!selfSolve) {
      advanceSolver->setJacobian(j);
    }
  }
  void setSplitOperator(rhsfunc fC, rhsfunc fD) override {
    if (selfSolve) {
      Solver::setSplitOperator(fC, fD);
    } else {
      advanceSolver->setSplitOperator(fC, fD);
    }
  }

  // Constraints
  bool constraints() override {
    if (selfSolve) {
      return false;
    } else {
      return advanceSolver->constraints();
    }
  }
  void constraint(Field2D& v, Field2D& C_v, std::string name) override {
    if (!selfSolve) {
      advanceSolver->constraint(v, C_v, std::move(name));
    }
  }
  void constraint(Field3D& v, Field3D& C_v, std::string name) override {
    if (!selfSolve) {
      advanceSolver->constraint(v, C_v, std::move(name));
    }
  }
  void constraint(Vector2D& v, Vector2D& C_v, std::string name) override {
    if (!selfSolve) {
      advanceSolver->constraint(v, C_v, std::move(name));
    }
  }
  void constraint(Vector3D& v, Vector3D& C_v, std::string name) override {
    if (!selfSolve) {
      advanceSolver->constraint(v, C_v, std::move(name));
    }
  }

  // Override count operations
  int n2Dvars() const override {
    if (selfSolve) {
      return Solver::n2Dvars();
    } else {
      return advanceSolver->n2Dvars();
    }
  }
  int n3Dvars() const override {
    if (selfSolve) {
      return Solver::n3Dvars();
    } else {
      return advanceSolver->n3Dvars();
    }
  }
  // Time steps
  void setMaxTimestep(BoutReal dt) override {
    if (selfSolve) {
      Solver::setMaxTimestep(dt);
    } else {
      advanceSolver->setMaxTimestep(dt);
    }
  }
  BoutReal getCurrentTimestep() override {
    if (selfSolve) {
      return Solver::max_dt;
    }
    { return advanceSolver->getCurrentTimestep(); }
  }

  int compareState;

  void slepcToBout(PetscScalar &reEigIn, PetscScalar &imEigIn, BoutReal &reEigOut,
                   BoutReal &imEigOut, bool force = false);

  Mat shellMat; //"Shell" matrix operator
private:
  MPI_Comm comm;
  EPS eps; // Slepc solver handle

  ST st;               // Spectral transform object
  PetscBool stIsShell; // Is the ST a shell object?

  std::unique_ptr<Solver> advanceSolver{nullptr}; // Pointer to actual solver used to advance fields

  void vecToFields(Vec &inVec);
  void fieldsToVec(Vec &outVec);

  void createShellMat();
  void createEPS();

  void analyseResults();
  void boutToSlepc(BoutReal &reEigIn, BoutReal &imEigIn, PetscScalar &reEigOut,
                   PetscScalar &imEigOut, bool force = false);

  SlepcLib slib;  // Handles initialize / finalize
  bool ddtMode;   // If true then slepc deals with the ddt operator
  bool selfSolve; // If true then we don't create an advanceSolver
  bool eigenValOnly;

  // For selfSolve=true
  Array<BoutReal> f0, f1;

  // Timestep details
  int nout;
  BoutReal tstep;

  // Used for SLEPc options
  int nEig, maxIt;
  int mpd; // Maximum projected dimension
  PetscReal tol, target;
  BoutReal targRe, targIm;
  bool userWhich;

  // Generic options
  bool useInitial; // If true then set the first vector in subspace to be initial
                   // conditions
  bool debugMonitor;

  // Grid size
  PetscInt localSize;
};

#endif // BOUT_HAS_SLEPC

#endif // __SLEPC_SOLVER_H__
