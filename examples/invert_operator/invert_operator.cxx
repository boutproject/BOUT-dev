
#include <bout/invert_operator.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

Field3D minus(const Field3D &input) { return -1.0 * input; };
Field3D delp(const Field3D &input) { return input + Delp2(input); };

class HW : public PhysicsModel {
private:
  Field3D n;
  class InvertOperator<Field3D> *mySolver;

protected:
  int init(bool restart) {

    SOLVE_FOR(n);

    mySolver = new InvertOperator<Field3D>();

    mySolver->setup(delp);

    Field3D output = 3.0;
    n = -2.0;
    output = mySolver->invert(n);

    output_warn << endl;
    output_warn << "Max difference is " << max(abs(output - n), true) << std::endl;
    output_warn << std::endl;
    auto pass = mySolver->verify(n) == 0 ? "False" : "True";
    output_warn << "Has test passed ? " << pass << std::endl;
    output_warn << std::endl;

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
