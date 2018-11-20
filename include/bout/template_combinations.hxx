/*!************************************************************************
 * \file template_combinations.hxx
 *
 * Routines and helper types for calling a templated function with all combinations
 * of various sets of types/values.
 *
 **************************************************************************
 * Copyright 2018
 *    D.Dickinson, P.Hill
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __TEMPLATE_COMBINATIONS_H__
#define __TEMPLATE_COMBINATIONS_H__

#include <unused.hxx>

/// Here we define an empty templated struct that can represent
/// a collection of arbitrary types. This is useful for passing
/// template packs (typename...) around whilst being able to
/// distinguish between different template packs.
template <typename... Ts>
struct Set {};

/// Here we provide a container type that can be used to pass around
/// a type without needing to create instances of the specific type
/// (instead we create instances of the container type).
template <typename T>
struct TypeContainer {
  using type = T;
};

/// Define a struct (functor) that we use to build up the final
/// collection of template values. Each time we create one of these
/// objects we provide one more type to the templatePack. We may call
/// the operator() method by providing all the required template parameters
/// for whatever the storedFunc is at that point -- once we have built a
/// complete templatePack we don't need to specify any of these template
/// parameters as they can be deduced/inferred.
template <typename currentFunction, typename currentType>
struct DeferredFunction {
  // Just store the actual function we wish to apply
  explicit DeferredFunction(currentFunction f) : storedFunc(f){};

  // The actual function we wish to apply with the functor
  currentFunction storedFunc;

  // Make the struct a functor by defining operator() we use
  // type inference to populate the templatePack/args
  template <typename... templatePack>
  void operator()(templatePack... args) {
    storedFunc(currentType{}, args...);
  }
};

///////////////////////////////////////////////////////////
/// Now we define routines for dealing with Sets of types
/// We use recursion to unpack each Set such that we can
/// form all combinations of the contained types.
/// As we empty Sets and get individual items we change
/// the number of arguments and hence provide several
/// overloads with different numbers of Sets.
/// Finally we end up with a single item -- this is the
/// point at which we have a unique combination of the
/// template types and can final invoke the DeferredFunction
///
/// Note that we define the routines from the bottom up so
/// that we don't have to pre-declare routines.
///
/// Note we make use of type inferenence with templates to
/// be able to refer to the first item in a Set/pack allowing
/// us to extract items from Sets.
///
///////////////////////////////////////////////////////////

//--------------------------------------------------------
// Routines for extracting single items from Sets and
// adding to DeferredFunction.
//--------------------------------------------------------

/// This is the lowest level routine -- we now have a unique
/// combination of template parameters provided to the
/// DeferredFunction (completed by passing in `item`) and
/// can therefore invoke this functor.
template <typename item, typename theFunction>
void addItemToDeferredFunction(theFunction func, item) {
  // Create a new DeferredFunction building on the passed DeferredFunction (func)
  // and adding the next template parameter, `item`, to the stored pack.
  DeferredFunction<theFunction, item> theFinalFunction(func);

  // As we have no more Sets to process here we have a unique/complete
  // combination of the template items and hence can call the final DeferredFunction
  // This terminates this path and we move onto the next values in the Sets.
  theFinalFunction();
}

/// One Set left to process so no template pack required
template <typename item, typename lastSet, typename theFunction>
void addItemToDeferredFunction(theFunction func, item, lastSet) {
  // Create a new DeferredFunction building on the passed DeferredFunction (func)
  // and adding the next template parameter, `item`, to the stored pack.
  DeferredFunction<theFunction, item> theNextFunction(func);

  // Process the lastSet, passing in the current state of DeferredFunction.
  processSet(theNextFunction, lastSet{});
}

/// More than one Set left to process.
template <typename item, typename nextSet, typename... otherSets, typename theFunction>
void addItemToDeferredFunction(theFunction func, item, nextSet, otherSets...) {
  // Create a new DeferredFunction building on the passed DeferredFunction (func)
  // and adding the next template parameter, `item`, to the stored pack.
  DeferredFunction<theFunction, item> theNextFunction(func);

  // Process the nextSet, passing in the current state of DeferredFunction
  // and the other remaining Sets.
  processSet(theNextFunction, nextSet{}, otherSets{}...);
}

//--------------------------------------------------------
// Routines for extracting Sets from Set packs and
// passing onto the routines for processing items in a Set.
//--------------------------------------------------------

/// Terminal routine -- the current Set is empty
/// so nothing left to do.
template <typename... Sets, typename theFunction>
void processSet(theFunction UNUSED(func), Set<>, Sets...){};

/// Here we use type inference to allow us to refer to the firstItem in the first Set
/// and the otherItems in this Set. We use this to pass the firstItem off to the routines
/// that will add this to the DeferredFunction. Following this we use recursion to call
/// this
/// routine again to process the rest of this Set.
template <typename firstItem, typename... otherItems, typename... otherSets,
          typename theFunction>
void processSet(theFunction func, Set<firstItem, otherItems...>, otherSets... others) {
  // Take the firstItem out of the current (first) Set and add to the DeferredFunction
  addItemToDeferredFunction(func, firstItem{}, others...);

  // Invoke this routine again with the items left in this Set (i.e. without firstItem)
  processSet(func, Set<otherItems...>{}, others...);
}

//--------------------------------------------------------
// Routine(s) for kicking off the process of processing
// the Sets.
//--------------------------------------------------------

/// This is the top level routine that takes the different Sets of types
/// and triggers the construction of instances of theFunction with all the
/// combinations of the template types defined by the Sets.
///
/// A use of this might look like:
/// produceCombinations<
///     Set<typeA, typeB, typeC>,
///     Set<int, double, std::string>
///   >(someFunctionWhichTakesTwoTemplateTypeArguments);
///
/// Note we wrap this in a struct such that by declaring a global variable of this
/// type we trigger the creation of the combinations.
template <typename FirstSet, typename... otherSets>
struct produceCombinations {
  template <typename theFunction>
  explicit produceCombinations(theFunction func) {
    processSet(func, FirstSet{}, otherSets{}...);
  };
};
#endif
