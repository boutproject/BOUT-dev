/*
 * Tests reading of options
 * 
 */

#include <bout/physicsmodel.hxx>
#include <options.hxx>
#include <field_factory.hxx>

class TestFieldFactory2 : public PhysicsModel {
protected:
  // Initialisation
  int init(bool restarting) {
    // Create a field factory for parsing strings
    FieldFactory f(mesh);
    
    Options *opt = Options::getRoot();

    string astr;
    opt->get("a", astr, "0.0");
    
    Field3D a = f.create3D(astr, opt);
    
    return 1; // Quit
  }
};

// Create a default main()
BOUTMAIN(TestFieldFactory2);
