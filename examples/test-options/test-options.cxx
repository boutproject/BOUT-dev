
#include <bout.hxx>
#include <assert.h>
#include <math.h>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Options *opt = Options::getRoot();

  assert( opt != nullptr );
  assert( !opt->isSet("test") );

  opt->set("test", "a string", "test source");
  assert( opt->isSet("test") );

  string str;
  opt->get("test", str, "default value");
  assert( str == "a string" );

  BoutReal val;
  opt->get("value", val, 0.0);
  assert( fabs(val) < 1e-10 );

  // Repeat, to test that the value has been set correctly
  opt->get("value", val, 0.0);
  assert( fabs(val) < 1e-10 );

  try { // Try to read with a different default
   opt->get("value", val, 1.0);
   assert(false); // Should never get to here
  } catch (BoutException &e) {
    // success
  }

  BoutFinalise();
  return 0;
}
