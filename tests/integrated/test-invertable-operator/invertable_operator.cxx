#include <invert_laplace.hxx>
#include <msg_stack.hxx>

#include <bout/invertable_operator.hxx>
#include <bout/physicsmodel.hxx>

class InvertableOperatorTest : public PhysicsModel {
private:
  Field3D n, solutionInv, solutionLap;

  int passVerification = 0;
  BoutReal maxRelErrLaplacians = 0;

  struct myLaplacian {
    Field3D D = 1.0, A = 1.0e-1;

    // Drop C term for now
    Field3D operator()(const Field3D& input) {
      TRACE("myLaplacian::operator()");
      Field3D result = A * input + D * Delp2(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };

  myLaplacian mm;
  bout::inversion::InvertableOperator<Field3D> mySolver;

  class Laplacian* laplacianSolver;

protected:
  int init(bool) override {
    SOLVE_FOR(n);
    SOLVE_FOR(solutionLap);
    SOLVE_FOR(solutionInv);

    SAVE_REPEAT(passVerification);
    SAVE_REPEAT(maxRelErrLaplacians);

    mm.A = 1.0e-1;
    mm.D = 1.0;

    // Note mySolve takes a copy of the passed functor so updates to the local
    // instance won't have any effect, but the function _can_ be changed (currently)
    // through setOperatorFunction

    mySolver.setOperatorFunction(mm);
    mySolver.setup();

    laplacianSolver = Laplacian::create();
    laplacianSolver->setCoefA(mm.A);
    laplacianSolver->setCoefC(1.0);
    laplacianSolver->setCoefD(mm.D);

    n.applyBoundary("dirichlet");

    return 0;
  }

  int rhs(BoutReal) override {
    ddt(n) = 0.;
    ddt(solutionInv) = 0.;
    ddt(solutionLap) = 0.;

    // First run to get the solution with zero initial guess
    solutionInv = mySolver.invert(n, 0.0);
    mesh->communicate(solutionInv);

    passVerification = mySolver.verify(n, 1.e-3);

    solutionLap = laplacianSolver->solve(n);

    maxRelErrLaplacians =
        max(abs(solutionInv - solutionLap), true) / max(abs(solutionLap));

    return 0;
  }

public:
};

// Define a main() function
BOUTMAIN(InvertableOperatorTest);
