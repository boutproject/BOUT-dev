/*!
 *
 */

#include <bout/array.hxx>

#include <bout_types.hxx>
#include <dcomplex.hxx>

template<>
std::map< int, std::vector<Array<double>::ArrayData* > > Array<double>::store = {}; // NB: C++11

template<>
bool Array<double>::use_store = true;

template<>
std::map< int, std::vector<Array<dcomplex>::ArrayData* > > Array<dcomplex>::store = {};

template<>
bool Array<dcomplex>::use_store = true;
