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
  int init(bool UNUSED(restarting)) {
    // Create a field factory for parsing strings
    FieldFactory f(mesh);
    
    Options *opt = Options::getRoot();

    Field3D a = f.create3D("a", opt);
    
    return 1; // Quit
  }
};

// Create a default main()
BOUTMAIN(TestFieldFactory2);
