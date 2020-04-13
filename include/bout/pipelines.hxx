
#pragma once
#ifndef PIPELINES_H
#define PIPELINES_H

#include <functional>

#include "../field3d.hxx"

namespace bout {
  namespace experimental {
  
    using GeneratorFunction = std::function<BoutReal(const Ind3D&)>;
    
    GeneratorFunction generate(BoutReal value) {
      return [value](const Ind3D&) { return value; };
    }
    
    GeneratorFunction generate(const Field3D &field) {
      return [&](const Ind3D& ind) { return field[ind]; };
    }

    auto generate(GeneratorFunction func) {
      return func;
    }
    
    // Convert a Field into a generator
    template<typename F>
    auto operator|(const F &f, std::function< Field3D (GeneratorFunction)> &&func) {
      return func(generate(f));
    }

    template<typename F>
    auto operator|(const F &f, std::function< GeneratorFunction (GeneratorFunction)> &&func) {
      return func(generate(f));
    }

    // Bind one generator to the next
    template<typename T>
    T operator|(GeneratorFunction left, std::function< T (GeneratorFunction)> &&right) {
      return right(left);
    }

    // Iterate and evaluate generators
    auto collect() {
      return [](GeneratorFunction func) {
        Field3D result;
        result.allocate();
        BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = func(i); }
        return result;
      };
    }

    template<typename T>
    auto add(const T &value) {
      auto adding = generate(value);
      // Function which takes and returns a generator
      return [adding](GeneratorFunction func) {
               // The function to calculate the combined result
               return [=](const Ind3D& ind) {
                        return func(ind) + adding(ind);
                      };
             };
    }

    
  } // experimental
} // bout
  
#endif // PIPELINES_H
