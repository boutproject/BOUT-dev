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


#ifdef UNIT
/*
 * Unit testing of Array class
 * 
 *
 */

#include <iostream>
using std::cout;

#include <assert.h>

int main() {
  Array<double> a(10);

  assert(!a.empty());      // Not empty
  assert(a.size() == 10);  // Correct size
  assert(a.unique());      // Should be unique

  // Set some values
  
  int count = 0;
  for(auto &i : a)
    i = count++;

  assert(a[1] = 1);
  assert(a[9] = 9);

  Array<double> b(a);   // Copy constructor

  assert(!b.empty());
  assert(b.size() == 10);
  assert(b[5] == 5);
  assert(!b.unique());
  assert(!a.unique());

  // Make both b and a unique
  b.ensureUnique();
  
  assert(b.unique());
  assert(a.unique());
  assert(a.size() == 10);
  assert(b.size() == 10);
  // Should have the same values
  for(auto ai = a.begin(), bi = b.begin();
      ai != a.end(); ++ai, ++bi)
    assert(*ai == *bi);

  // Release the data. Should put into store.
  a.clear();
  assert(a.empty());
  assert(!b.empty());
  assert(a.size() == 0);

  // Construct, retrieve from store, and move assign
  a = Array<double>(10);

  assert(!a.empty());
  assert(a.size() == 10);
  assert(a[4] == 4); // Test if reused data from store
  assert(a.unique());

  a.ensureUnique(); // Should have no effect
  
  assert(a.size() == 10);
  assert(a[4] == 4); // Test if reused data from store
  assert(a.unique());

  // Assign
  a = b;
  
  assert(!a.unique());
  assert(!b.unique());
  
  return 0;
}

#endif
