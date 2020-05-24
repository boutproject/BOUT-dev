/**************************************************************************
 * Perpendicular Laplacian inversion in X-Z
 *
 * Equation solved is:
 *
 * Div( A * Grad_perp(x) ) + B*x = b
 *
 *
 **************************************************************************
 * Copyright 2015 B.D.Dudson
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

#ifndef __LAPLACEXZ_H__
#define __LAPLACEXZ_H__

#include <options.hxx>
#include <field3d.hxx>
#include <bout/mesh.hxx>
#include <unused.hxx>
#include <bout/generic_factory.hxx>

class LaplaceXZ;

class LaplaceXZFactory
    : public Factory<
          LaplaceXZ, LaplaceXZFactory,
          std::function<std::unique_ptr<LaplaceXZ>(Mesh*, Options*, CELL_LOC)>> {
public:
  static constexpr auto type_name = "LaplaceXZ";
  static constexpr auto section_name = "laplacexz";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "cyclic";

  ReturnType create(Mesh* mesh = nullptr, Options* options = nullptr,
                    CELL_LOC loc = CELL_CENTRE) {
    return Factory::create(getType(options), mesh, options, loc);
  }

  static void ensureRegistered();
};

template <class DerivedType>
class RegisterLaplaceXZ {
public:
  RegisterLaplaceXZ(const std::string& name) {
    LaplaceXZFactory::getInstance().add(
      name, [](Mesh* mesh, Options* options, CELL_LOC loc) -> std::unique_ptr<LaplaceXZ> {
        return std::make_unique<DerivedType>(mesh, options, loc);
      });
  }
};

class LaplaceXZ {
public:
  LaplaceXZ(Mesh* m = nullptr, Options* UNUSED(options) = nullptr,
            const CELL_LOC loc = CELL_CENTRE)
      : localmesh(m == nullptr ? bout::globals::mesh : m), location(loc) {}
  virtual ~LaplaceXZ() = default;

  virtual void setCoefs(const Field2D &A, const Field2D &B) = 0;
  virtual void setCoefs(const Field3D &A, const Field3D &B) { setCoefs(DC(A), DC(B)); }

  virtual Field3D solve(const Field3D &b, const Field3D &x0) = 0;

  static std::unique_ptr<LaplaceXZ> create(Mesh* m = nullptr, Options* opt = nullptr,
                                           CELL_LOC loc = CELL_CENTRE) {
    return LaplaceXZFactory::getInstance().create(m, opt, loc);
  }

protected:
  static const int INVERT_DC_GRAD  = 1;
  static const int INVERT_AC_GRAD  = 2;  // Use zero neumann (NOTE: AC is a misnomer)
  static const int INVERT_SET      = 16; // Set boundary to x0 value
  static const int INVERT_RHS      = 32; // Set boundary to b value
  Mesh* localmesh;   ///< The mesh this operates on, provides metrics and communication
  CELL_LOC location;
private:

};

#endif // __LAPLACEXZ_H__
