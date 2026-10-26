/**************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_NVECTOR_H
#define BOUT_NVECTOR_H

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/build_defines.hxx>
#include <bout/field.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <functional>
#if BOUT_HAS_SUNDIALS_MANYVECTOR
#include <nvector/nvector_manyvector.h>
#endif
#include <sundials/sundials_nvector.h>
#include <utility>

#include "bout/sundials_backports.hxx"

class BoutNVector {
private:
  static double all_reduce(const double local, const MPI_Op op = MPI_SUM) {
    double global;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, op, BoutComm::get());
    return global;
  }

  static double all_reduce(const unsigned int local, const MPI_Op op = MPI_SUM) {
    unsigned int global;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED, op, BoutComm::get());
    return global;
  }

  template <typename T>
  struct Content {
    T& field;
    const bool own;
    const bool evolve_bndry;
    const sunindextype length;

    Content(T& field, const bool own, const bool evolve_bndry)
        : field(field), own(own), evolve_bndry(evolve_bndry),
          length(all_reduce(getRegion().size())) {}

    auto getRegion() const {
      return field.getRegion(evolve_bndry ? RGN_ALL : RGN_NOBNDRY);
    }

    ~Content() {
      if (own) {
        delete &field;
      }
    }
  };

  template <typename T>
  static Content<T>& get_content(const N_Vector v) {
    return *static_cast<Content<T>*>(v->content);
  }

  template <typename T>
  static T& get_field(const N_Vector v) {
    return get_content<T>(v).field;
  }

  template <typename T, typename Reduce, typename... V>
  static BoutReal reduce_field(Reduce&& reduce, const bool allpe, const N_Vector v1,
                               const V... v2) {
    BoutReal result = 0;
    const auto region = get_content<T>(v1).getRegion();
    BOUT_FOR_OMP(i, region, parallel for reduction(+:result)) {
      result += std::invoke(std::forward<Reduce>(reduce), get_field<T>(v1)[i],
                            get_field<T>(v2)[i]...);
    }

    return allpe ? all_reduce(result) : result;
  }

public:
  BoutNVector() = delete; // Enforce static access only

  template <IsField T, typename Ctx>
  static N_Vector create(Ctx&& ctx, T& field, const bool evolve_bndry,
                         const bool own = false) {
    N_Vector v = callWithSUNContext(N_VNewEmpty, std::forward<Ctx>(ctx));
    if (v == nullptr) {
      throw BoutException("N_VNewEmpty failed\n");
    }

    v->content = static_cast<void*>(new Content<T>(field, own, evolve_bndry));

    v->ops->nvgetvectorid = [](N_Vector) { return SUNDIALS_NVEC_CUSTOM; };

    v->ops->nvclone = [](N_Vector x) {
      const Content<T>& content = get_content<T>(x);
      // TODO: ensure no memory leaks
      // T* field_clone = new T(*content.field);
      T* field_clone = new T(content.field.getMesh(), content.field.getLocation(),
                             content.field.getDirections());
      field_clone->allocate();
      field_clone->copyBoundary(content.field);
      return create(x->sunctx, *field_clone, content.evolve_bndry, true);
    };

    v->ops->nvdestroy = [](N_Vector x) {
      if (x == nullptr) {
        return;
      }
      delete &get_content<T>(x);
      N_VFreeEmpty(x);
    };

    v->ops->nvgetlength = [](N_Vector x) { return get_content<T>(x).length; };

    v->ops->nvlinearsum = [](sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                             N_Vector z) {
      get_field<T>(z) = a * get_field<T>(x) + b * get_field<T>(y);
    };

    v->ops->nvconst = [](sunrealtype c, N_Vector x) { get_field<T>(x) = c; };

    v->ops->nvprod = [](N_Vector x, N_Vector y, N_Vector z) {
      get_field<T>(z) = get_field<T>(x) * get_field<T>(y);
    };

    v->ops->nvdiv = [](N_Vector x, N_Vector y, N_Vector z) {
      get_field<T>(z) = get_field<T>(x) / get_field<T>(y);
    };

    v->ops->nvscale = [](sunrealtype a, N_Vector x, N_Vector y) {
      get_field<T>(y) = a * get_field<T>(x);
    };

    v->ops->nvabs = [](N_Vector x, N_Vector y) {
      get_field<T>(y) = abs(get_field<T>(x));
    };

    v->ops->nvinv = [](N_Vector x, N_Vector y) { get_field<T>(y) = 1 / get_field<T>(x); };

    v->ops->nvaddconst = [](N_Vector x, sunrealtype a, N_Vector y) {
      get_field<T>(y) = a + get_field<T>(x);
    };

    v->ops->nvdotprod = [](N_Vector x, N_Vector y) {
      return reduce_field<T>(std::multiplies<BoutReal>(), true, x, y);
    };

    v->ops->nvmaxnorm = [](N_Vector x) { return max(abs(get_field<T>(x)), true); };

    v->ops->nvwrmsnorm = [](N_Vector x, N_Vector y) {
      return N_VWL2Norm(x, y) / std::sqrt(static_cast<BoutReal>(N_VGetLength(x)));
    };

    v->ops->nvmin = [](N_Vector x) { return min(get_field<T>(x), true); };

    v->ops->nvwl2norm = [](N_Vector x, N_Vector y) {
      return std::sqrt(reduce_field<T>(
          [](const auto v1, const auto v2) { return std::pow(v1 * v2, 2); }, true, x, y));
    };

    v->ops->nvl1norm = [](N_Vector x) {
      return reduce_field<T>([](const auto value) { return std::abs(value); }, true, x);
    };

    return v;
  }

  template <IsField T>
  static void swap(const N_Vector v, T& field) {
    field.swapData(get_field<T>(v));
  }

  template <IsField T>
  static T& get(const N_Vector v) {
    return get_field<T>(v);
  }

#if BOUT_HAS_SUNDIALS_MANYVECTOR
  template <typename V>
  static N_Vector create(const sundials::Context& ctx, V& subvectors) {
    const auto v = callWithSUNContext(N_VNew_ManyVector, ctx, std::size(subvectors),
                                      std::data(subvectors));
    ((N_VectorContent_ManyVector)v->content)->own_data = true;
    return v;
  }

  template <IsField T>
  static void swap(const N_Vector v, T& field, std::size_t subvector) {
    return BoutNVector::swap(N_VGetSubvector_ManyVector(v, subvector), field);
  }

  template <IsField T>
  static T& get(const N_Vector v, std::size_t subvector) {
    return BoutNVector::get<T>(N_VGetSubvector_ManyVector(v, subvector));
  }
#endif
};

#endif

/* Other functions to implement (most optional)
N_Vector (*nvcloneempty)(N_Vector); // Probably not needed
void (*nvdestroy)(N_Vector); // Implemented
void (*nvspace)(N_Vector, sunindextype*, sunindextype*); // Not needed
sunrealtype* (*nvgetarraypointer)(N_Vector); // Probably not needed
sunrealtype* (*nvgetdevicearraypointer)(N_Vector); // Probably not needed
void (*nvsetarraypointer)(sunrealtype*, N_Vector); // Probably not needed
SUNComm (*nvgetcommunicator)(N_Vector); // Probably not needed
sunindextype (*nvgetlength)(N_Vector); // Implemented
sunindextype (*nvgetlocallength)(N_Vector);  // Probably not needed
void (*nvlinearsum)(sunrealtype, N_Vector, sunrealtype, N_Vector, N_Vector); // Implemented
void (*nvconst)(sunrealtype, N_Vector); // Implemented
void (*nvprod)(N_Vector, N_Vector, N_Vector); // Implemented
void (*nvdiv)(N_Vector, N_Vector, N_Vector); // Implemented
void (*nvscale)(sunrealtype, N_Vector, N_Vector); // Implemented
void (*nvabs)(N_Vector, N_Vector); // Implemented
void (*nvinv)(N_Vector, N_Vector); // Implemented
void (*nvaddconst)(N_Vector, sunrealtype, N_Vector); // Implemented
sunrealtype (*nvdotprod)(N_Vector, N_Vector); // Implemented
sunrealtype (*nvmaxnorm)(N_Vector); // Implemented
sunrealtype (*nvwrmsnorm)(N_Vector, N_Vector); // Implemented
sunrealtype (*nvwrmsnormmask)(N_Vector, N_Vector, N_Vector); // TODO
sunrealtype (*nvmin)(N_Vector); // Implemented
sunrealtype (*nvwl2norm)(N_Vector, N_Vector); // Implemented
sunrealtype (*nvl1norm)(N_Vector); // Implemented
void (*nvcompare)(sunrealtype, N_Vector, N_Vector); // TODO
sunbooleantype (*nvinvtest)(N_Vector, N_Vector); // TODO
sunbooleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector); // TODO
sunrealtype (*nvminquotient)(N_Vector, N_Vector); // TODO
SUNErrCode (*nvlinearcombination)(int, sunrealtype*, N_Vector*, N_Vector); // TODO
SUNErrCode (*nvscaleaddmulti)(int, sunrealtype*, N_Vector, N_Vector*,
                          N_Vector*); // Probably not needed
SUNErrCode (*nvdotprodmulti)(int, N_Vector, N_Vector*, sunrealtype*); // Probably not needed
SUNErrCode (*nvlinearsumvectorarray)(int, sunrealtype, N_Vector*, sunrealtype,
                                    N_Vector*, N_Vector*); // Probably not needed
SUNErrCode (*nvscalevectorarray)(int, sunrealtype*, N_Vector*, N_Vector*); // Probably not needed
SUNErrCode (*nvconstvectorarray)(int, sunrealtype, N_Vector*); // Probably not needed
SUNErrCode (*nvwrmsnormvectorarray)(int, N_Vector*, N_Vector*, sunrealtype*); // Probably not needed
SUNErrCode (*nvwrmsnormmaskvectorarray)(int, N_Vector*, N_Vector*, N_Vector,
                                      sunrealtype*); // Probably not needed
SUNErrCode (*nvscaleaddmultivectorarray)(int, int, sunrealtype*, N_Vector*,
                                      N_Vector**, N_Vector**); // Probably not needed
SUNErrCode (*nvlinearcombinationvectorarray)(int, int, sunrealtype*,
                                          N_Vector**, N_Vector*); // Probably not needed
sunrealtype (*nvdotprodlocal)(N_Vector, N_Vector); // Probably not needed
sunrealtype (*nvmaxnormlocal)(N_Vector); // Probably not needed
sunrealtype (*nvminlocal)(N_Vector); // Probably not needed
sunrealtype (*nvl1normlocal)(N_Vector); // Probably not needed
sunbooleantype (*nvinvtestlocal)(N_Vector, N_Vector); // Probably not needed
sunbooleantype (*nvconstrmasklocal)(N_Vector, N_Vector, N_Vector); // Probably not needed
sunrealtype (*nvminquotientlocal)(N_Vector, N_Vector); // Probably not needed
sunrealtype (*nvwsqrsumlocal)(N_Vector, N_Vector); // Probably not needed
sunrealtype (*nvwsqrsummasklocal)(N_Vector, N_Vector, N_Vector); // Probably not needed
SUNErrCode (*nvdotprodmultilocal)(int, N_Vector, N_Vector*, sunrealtype*); // Probably not needed
SUNErrCode (*nvdotprodmultiall_reduce)(int, N_Vector, sunrealtype*); // Probably not needed
SUNErrCode (*nvbufsize)(N_Vector, sunindextype*); // Probably not needed
SUNErrCode (*nvbufpack)(N_Vector, void*); // Probably not needed
SUNErrCode (*nvbufunpack)(N_Vector, void*); // Probably not needed
void (*nvprint)(N_Vector); // Nice, but optional
void (*nvprintfile)(N_Vector, FILE*); // Probably not needed
*/
