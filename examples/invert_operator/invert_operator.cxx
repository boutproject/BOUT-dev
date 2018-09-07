
#include <bout/invert_operator.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

Field3D minus(const Field3D &input) { return -1.0 * input; };
Field3D delp(const Field3D &input) { return input + Delp2(input); };

class HW : public PhysicsModel {
private:
  Field3D n;
  InvertOperator<Field3D> *mySolver;
  
  struct myOp : public OperatorWrapper {
    BoutReal factor=1.;
    Field3D operator()(const Field3D &input) override {output<<"Factor is "<<factor<<endl; return factor*input + Delp2(input); } ;
  };
  myOp myDelp;

protected:
  int init(bool restart) {

    SOLVE_FOR(n);

    mySolver = new InvertOperator<Field3D>();

    // mySolver->setup(delp);
    mySolver->setup(myDelp);
    // Above could also be:
    // mySolver->setup(delp);
    // or even a Lambda
    // mySolver->setup([](const Field3D &input) { return input + Delp2(input); });

    // Note mySolve takes a copy of the passed functor so updates to the local
    // instance won't have any effect, but the function _can_ be changed (currently)
    // as it is a public member, so the following should work.
    myDelp.factor = 3.0;
    mySolver->func = myDelp;
    Field3D output = 3.0;
    n = -2.0;
    for(int i=0; i<10000; i++){
      output = mySolver->invert(n);
    }
    output_warn << endl;
    output_warn << "Max difference is " << max(abs(output - n), true) << std::endl;
    output_warn << std::endl;
    auto pass = mySolver->verify(n) == 0 ? "False" : "True";
    output_warn << "Has test passed ? " << pass << std::endl;
    output_warn << std::endl;

    mySolver->reportTime();
    return 0;
  }

  int rhs(BoutReal time) {
    ddt(n) = 0.;
    return 0;
  }

public:
  ~HW() { delete mySolver; };
};

// Define a main() function
BOUTMAIN(HW);
