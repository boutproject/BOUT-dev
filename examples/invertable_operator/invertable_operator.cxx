#include <invert_laplace.hxx>
#include <msg_stack.hxx>

#include <bout/sys/timer.hxx>
#include <bout/invertable_operator.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

#include <boundary_region.hxx>

Field3D minus(const Field3D &input) { return -1.0 * input; };
Field3D delp(const Field3D &input) { return input + Delp2(input); };

class HW : public PhysicsModel {
private:
  Field3D n, solutionInv, solutionLap;
  
  struct myOp {
    BoutReal factor=1.;
    Field3D operator()(const Field3D &input) {return factor*input + Delp2(input); } ;
  };
  myOp myDelp;

  struct myLaplacian {
    Field3D D=1.0, C=1.0, A=0.0;

    // Drop C term for now
    Field3D operator()(const Field3D &input) {
      TRACE("myLaplacian::operator()");
      Field3D result = A*input + D*Delp2(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };

  
  myLaplacian mm;
  InvertableOperator<Field3D> mySolver;
  // Above could also be:
  // mySolver(delp);
  // or even a Lambda
  // mySolver([](const Field3D &input) { return input + Delp2(input); });

  class Laplacian *phiSolverXZ;
  
  const int nits=1000;    
  
protected:
  int init(bool restart) {
    SOLVE_FOR(n);
    SOLVE_FOR(solutionLap);
    SOLVE_FOR(solutionInv);

    mm.A=1.0e-1;
    mm.D=1.0;

    // Note mySolve takes a copy of the passed functor so updates to the local
    // instance won't have any effect, but the function _can_ be changed (currently)
    // through setOperatorFunction

    mySolver.setOperatorFunction(mm);
    mySolver.setup();


    phiSolverXZ = Laplacian::create();
    phiSolverXZ->setCoefA(mm.A);
    phiSolverXZ->setCoefC(1.0);
    phiSolverXZ->setCoefD(mm.D);

    n.applyBoundary("dirichlet");
    
    return 0;
  }

  int rhs(BoutReal time) {
    ddt(n) = 0.;
    ddt(solutionInv) = 0.;
    ddt(solutionLap) = 0.;

    try{
    for(int i=0; i<nits; i++){
      solutionInv = mySolver.invert(n);
    }
    }catch(BoutException){
    };
    
    output_warn << std::endl;
    auto pass = mySolver.verify(n) == 0 ? "False" : "True";
    output_warn << "Has test passed ? " << pass << std::endl;
    output_warn << std::endl;

    {
      Timer timer("sol_lap");
      try{
      for(int i=0; i<nits; i++){    
    	solutionLap = phiSolverXZ->solve(n);
      }
      }catch(BoutException){
      };
    }

    output_warn<<"Max diff undo Invertable is "<<max(abs(mm(solutionInv)-n),true)<<endl;
    output_warn<<"MAX DIFF SOL "<<max(abs(solutionLap-solutionInv),true)<<endl;

    mySolver.reportTime();
    output_warn<<"Laplacian time is "<<Timer::resetTime("sol_lap")<<endl;

    return 0;
  }

public:
};

// Define a main() function
BOUTMAIN(HW);
