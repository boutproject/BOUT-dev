#ifdef BOUT_HAS_PYBIND11
#include <pybind11/numpy.h>
#include <pybind11/embed.h>

#include <iostream>
#include <type_traits>

namespace py = pybind11;
using namespace py::literals;

// Provide a simple function that returns an int
int myFunc(){
  return 42;
};

template<typename In, typename Out>
struct has_cast {
private:
  // Routine to check if the result of the cast of In to Out works --
  // this fails if the types are different but also if In doesn't have
  // a template cast method.  The type of this will either be
  // std::false_type or std::true_type depending on if In and Out are
  // compatible
  template<typename T>
  static constexpr typename std::is_same<decltype( std::declval<T>().template cast<Out>()), Out>::type check(T*);

  // Fall back check type should the above check fail to be
  // available (i.e. because In doesn't have a valid cast method)
  template<typename>
  static constexpr std::false_type check(...);

  // Find the type of invoking the check, this will either be true
  // or false.
  using type = decltype(check<In>(nullptr));

public:
  static constexpr bool value = type::value;
};

template<typename U, typename T>
void getVal(U in, T& out) {
  static_assert(has_cast<U, T>::value,"Error: Can't cast U to T in getVal, did you pass a pybind11 object as the first argument?");
  out = in.template cast<T>();
  return;
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
    py::print();    
  }

  {
    // Get access to another python interpreter
    py::scoped_interpreter guard{};

    py::module sys = py::module::import("sys");
    py::print("The python path is:");
    py::print(sys.attr("path"));
    py::print("");

    py::module np = py::module::import("numpy");
    py::print(np);
    py::print("Numpy version is", np.attr("__version__"));
    auto linspaceTest = np.attr("linspace")(0.0,1.0,20);
    py::print(linspaceTest);
    py::print();    
    py::print(linspaceTest[py::make_tuple(2)]);
    linspaceTest[py::make_tuple(2)] = 2.0;
    py::print(linspaceTest[py::make_tuple(2)]);
    py::print();
    for(int i=0; i<20; i++) {
      linspaceTest[py::make_tuple(i)] = i;
    }
    py::print(linspaceTest);
    py::module matplot = py::module::import("matplotlib.pyplot");
    matplot.attr("ioff")();
    matplot.attr("plot")(linspaceTest,"-o");
    matplot.attr("xlabel")("This is the x axis label");
    matplot.attr("ylabel")("$\\alpha$ is alpha");
    matplot.attr("show")();    
    double theValue = linspaceTest[py::make_tuple(10)].cast<double>();
    std::cout<<"The value is "<<theValue<<std::endl;    
    getVal(linspaceTest[py::make_tuple(5)], theValue);
    std::cout<<"The value is "<<theValue<<std::endl;
    auto tmp = linspaceTest.cast<py::array_t<double>>();
    for (const auto &i: tmp){
      std::cout<<" with value "<< i.cast<double>() << std::endl;
    }

    // Here's another way of indexing the python array
    auto ndArray = np.attr("ones")(py::make_tuple(20,30,40,50));    
    auto buf = ndArray.cast<py::array_t<double>>().request();
    std::cout<<"Element 3 from buffer with shape [ ";
    for(const auto &i: buf.shape){
      std::cout<<i<<", ";
    }
    std::cout<<"] is "<<((double *) buf.ptr)[3]<<std::endl;    
  }

}
#else
#error Must configure with `--with-pybind11` to be able to build pybind_embed example.
#endif
