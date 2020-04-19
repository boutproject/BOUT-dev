
#pragma once
#ifndef PIPELINES_H
#define PIPELINES_H

#include <functional>

#include "../field3d.hxx"

namespace bout {
  namespace experimental {

    /// Expressions are converted into a GeneratorFunction which takes an index,
    /// and returns a value (BoutReal). This could be by looking up the 
    /// value in a field, or performing some calculation.
    ///
    using GeneratorFunction = std::function<BoutReal(const Ind3D&)>;

    /// The generate function converts values into GeneratorFunction objects
    GeneratorFunction generate(BoutReal value) {
      return [value](const Ind3D&) { return value; };
    }
    
    GeneratorFunction generate(const Field3D &field) {
      return [&](const Ind3D& ind) { return field[ind]; };
    }

    // Pass through generator functions
    auto generate(GeneratorFunction func) {
      return func;
    }
    
    /// Convert a Field into a generator
    ///
    /// Note: Variations on
    ///   template<typename F, typename T>
    ///    auto operator|(const F &f, std::function< T (GeneratorFunction)> func)
    /// don't work because T cannot be inferred from a lambda function
    /// 
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

    /// Iterate and evaluate generators
    ///
    /// This creates a function which evaluates a GeneratorFunction over a Field3D
    /// It is used in a pipeline when all values need to be evaluated e.g. at the end
    ///
    /// Example:
    ///
    /// Field3D result = 1.0 | add(2.0) | toField();
    ///     -> result is filled with 3.0
    ///
    /// Field3D result = field | slow_operation()
    ///                        | toField()
    ///                        | differentiate()
    ///                        | toField();
    ///
    /// where `differentiate` is an operation which needs to use stencils.
    /// If the first collect were not used, then `slow_operation` would be done
    /// more often than needed, every time a value was accessed.
    ///
    /// By default all points are evaluated (region "RGN_ALL"), but a region
    /// can be provided:
    ///
    ///   toField(mesh->getRegion3D("RGN_NOBNDRY"))
    /// 
    template<typename FieldType = Field3D>
    auto toField(const Region<typename FieldType::ind_type>& region = bout::globals::mesh->getRegion<FieldType>("RGN_ALL")) {
      return [&](GeneratorFunction func) {
        FieldType result;
        result.allocate();
        BOUT_FOR(i, region) { result[i] = func(i); }
        return result;
      };
    }

    template<typename T>
    auto add(const T &value) {
      // Convert the value into a generator
      auto adding = generate(value);
      // Function which takes and returns a generator
      return [adding](GeneratorFunction func) {
               // The function to calculate the combined result
               return [=](const Ind3D& ind) {
                        return func(ind) + adding(ind);
                      };
             };
    }

    template<typename T>
    auto subtract(const T &value) {
      // Convert the value into a generator
      auto subtracting = generate(value);
      // Function which takes and returns a generator
      return [subtracting](GeneratorFunction func) {
               // The function to calculate the combined result
               return [=](const Ind3D& ind) {
                        return func(ind) - subtracting(ind);
                      };
             };
    }

    template<typename T>
    auto multiply(const T &value) {
      // Convert the value into a generator
      auto multiplying = generate(value);
      // Function which takes and returns a generator
      return [multiplying](GeneratorFunction func) {
               // The function to calculate the combined result
               return [=](const Ind3D& ind) {
                        return func(ind) * multiplying(ind);
                      };
             };
    }

    template<typename T>
    auto divide(const T &value) {
      // Convert the value into a generator
      auto dividing = generate(value);
      // Function which takes and returns a generator
      return [dividing](GeneratorFunction func) {
               // The function to calculate the combined result
               return [=](const Ind3D& ind) {
                        return func(ind) / dividing(ind);
                      };
             };
    }
  } // experimental
} // bout
  
#endif // PIPELINES_H
