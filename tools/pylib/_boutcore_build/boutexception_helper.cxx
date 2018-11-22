#include "boutexception.hxx"
#include <Python.h>

void raise_bout_py_error() {
  try {
    throw;
  } catch (const BoutException& e) {
    PyErr_SetString(PyExc_RuntimeError, e.getBacktrace().c_str());
  }
}
    
