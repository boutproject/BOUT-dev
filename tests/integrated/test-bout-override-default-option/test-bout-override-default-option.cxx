#include<options.hxx>

// Use an integrated test for what is effectively a unit test of the
// BOUT_OVERRIDE_DEFAULT_OPTION() macro because the functionality relies on the state of
// the global options instance - in particular it will not work if Options::cleanup() is
// called before the overridden defaults are tested.

BOUT_OVERRIDE_DEFAULT_OPTION("OverrideDefaultValueOptionsMacro_str",
    "macro_override_value");
BOUT_OVERRIDE_DEFAULT_OPTION("OverrideDefaultValueOptionsMacro_int", 42);
BOUT_OVERRIDE_DEFAULT_OPTION("OverrideDefaultValueOptionsMacro_boutreal", 11.);

int main() {
  Options& options = Options::root();

  std::string value_str =
    options["OverrideDefaultValueOptionsMacro_str"].withDefault("macro_default_value");
  int value_int = options["OverrideDefaultValueOptionsMacro_int"].withDefault(1);
  BoutReal value_boutreal =
    options["OverrideDefaultValueOptionsMacro_boutreal"].withDefault(2.);

  bool success = true;
  if (value_str != "macro_override_value") {
    output_error << "value_str=" << value_str << " but should be macro_override_value"
      << endl;
    success = false;
  }
  if (value_int != 42) {
    output_error << "value_int=" << value_int << " but should be " << 42 << endl;
    success = false;
  }
  if (value_boutreal != 11.) {
    output_error << "value_boutreal=" << value_boutreal << " but should be " << 11.
      << endl;
    success = false;
  }

  // Return 0 for success=true, 1 for success=false
  return not success;
}
