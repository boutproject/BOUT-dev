#include <invert_laplace.hxx>
#include <msg_stack.hxx>

#include <bout/invertable_operator.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>
#include <derivs.hxx>

#include <boundary_region.hxx>

Field3D minus(const Field3D &input) { return -1.0 * input; };
Field3D delp(const Field3D &input) { return input + Delp2(input); };

class HW : public PhysicsModel {
private:
  Field3D n, solutionInv, solutionLap;

  struct myOp {
    BoutReal factor = 1.;
    Field3D operator()(const Field3D &input) { return factor * input + Delp2(input); };
  };
  myOp myDelp;

  struct myLaplacian {
    Field3D D = 1.0, C = 1.0, A = 0.0;

    // Drop C term for now
    Field3D operator()(const Field3D &input) {
      TRACE("myLaplacian::operator()");
      Timer timer("invertable_operator_operate");
      Field3D result = A * input + D * Delp2(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };

  struct my3DLaplacian {
    Field3D A = 1.0, B = 1.0, C = 0.0, D = 0.0;
    bool withDiv = false;

    // Drop C term for now
    Field3D operator()(const Field3D &input) {
      TRACE("myLaplacian::operator()");
      Timer timer("invertable_operator_operate");
      Field3D result = A * input + B * Laplace_perp(input);
      if (withDiv) {
        auto tmp = C * Grad_perp(input);
        input.getMesh()->communicate(tmp);
        result += Div(tmp);
      }
      result += D * Laplace(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };

  myLaplacian mm;
  bout::inversion::InvertableOperator<Field3D> mySolver;
  // Above could also be:
  // mySolver(delp);
  // or even a Lambda
  // mySolver([](const Field3D &input) { return input + Delp2(input); });

  class Laplacian* laplacianSolver;

  const int nits = 10;

protected:
  int init(bool restart) {
    SOLVE_FOR(n);
    SOLVE_FOR(solutionLap);
    SOLVE_FOR(solutionInv);

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

  int rhs(BoutReal time) {
    ddt(n) = 0.;
    ddt(solutionInv) = 0.;
    ddt(solutionLap) = 0.;

    // First run to get the solution with zero initial guess
    solutionInv = mySolver.invert(n, 0.0);
    mySolver.reportTime();

    try {
      for (int i = 0; i < nits; i++) {
        //    solutionInv = mySolver.invert(n);
        solutionInv = mySolver.invert(n, solutionInv);
        // mesh->communicate(solutionInv);
      }
    } catch (BoutException &e) {
    };

    mesh->communicate(solutionInv);

    output << std::endl;
    auto pass = mySolver.verify(n) == 0 ? "False" : "True";
    output << "Has test passed ? " << pass << std::endl;
    output << std::endl;

    {
      Timer timer("sol_lap");
      try {
        for (int i = 0; i < nits; i++) {
          solutionLap = laplacianSolver->solve(n);
        }
      } catch (BoutException &e) {
      };
    }

    output << "Max diff undo Invertable is " << max(abs(mm(solutionInv) - n), true)
           << endl;
    output << "MAX DIFF SOL " << max(abs(solutionLap - solutionInv), true) << endl;

    mySolver.reportTime();
    output << "Laplacian time is " << Timer::resetTime("sol_lap") << endl;

    return 0;
  }

public:
};

// Define a main() function
BOUTMAIN(HW);
