/*
 * Tests reading of options
 * 
 */

#include <bout/physicsmodel.hxx>
#include <bout/options.hxx>
#include <bout/field_factory.hxx>

class TestFieldFactory2 : public PhysicsModel {
protected:
  // Initialisation
  int init(bool restarting) {
    // Create a field factory for parsing strings
    FieldFactory f(mesh);
    
    Options *opt = Options::getRoot();

    Field3D a = f.create3D("a", opt);
    
    return 1; // Quit
  }
};

// Create a default main()
BOUTMAIN(TestFieldFactory2);
