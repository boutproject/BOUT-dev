#ifdef BOUT_HAS_PYBIND11
#include <pybind11/embed.h>

namespace py = pybind11;
using namespace py::literals;

// Provide a simple function that returns an int
int myFunc(){
  return 42;
};

int main() {
  {
    // Get access to a python interpreter with lifetime of the current scope
    py::scoped_interpreter guard{};

    // Construct a python dictionary -- note this includes a call
    // to a native C++ method as well as including constants.
    auto kwargs = py::dict("name"_a="BOUT++ user",
                           "number"_a=myFunc(),
                           "floatNumber"_a=1.2345                         
                           );
    // Format the string using the dictionary constructed above
    auto message = "Hello, {name}! The answer is {number} and not {floatNumber}"_s.format(**kwargs);
    // Display message
    py::print(message);
  }

  {
    // Get access to another python interpreter
    py::scoped_interpreter guard{};

    py::module sys = py::module::import("sys");
    py::print("The python path is:");
    py::print(sys.attr("path"));
  }
}
#else
#error Must configure with `--with-pybind11` to be able to build pybind_embed example.
#endif
