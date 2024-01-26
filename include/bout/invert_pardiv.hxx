/************************************************************************
 * Inversion of parallel divergence
 * Intended for use in SNB heat flux calculation
 *
 * Inverts a matrix of the form
 *
 * A + Div_par( B * Grad_par )
 *
 **************************************************************************
 * Copyright 2010-2022 BOUT++ contributors
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
 ************************************************************************/

#ifndef INV_PARDIV_H
#define INV_PARDIV_H

#include "field2d.hxx"
#include "field3d.hxx"
#include "options.hxx"
#include "unused.hxx"
#include "bout/generic_factory.hxx"

// Pardivergence implementations
constexpr auto PARDIVCYCLIC = "cyclic";

class InvertParDiv;

class InvertParDivFactory
    : public Factory<InvertParDiv, InvertParDivFactory, Options*, CELL_LOC, Mesh*> {
public:
  static constexpr auto type_name = "InvertParDiv";
  static constexpr auto section_name = "pardiv";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = PARDIVCYCLIC;

  ReturnType create(Options* options = nullptr, CELL_LOC location = CELL_CENTRE,
                    Mesh* mesh = nullptr) const {
    return Factory::create(getType(options), options, location, mesh);
  }
  ReturnType create(const std::string& type, Options* options) const {
    return Factory::create(type, options, CELL_CENTRE, nullptr);
  }
  static void ensureRegistered();
};

template <class DerivedType>
using RegisterInvertParDiv = InvertParDivFactory::RegisterInFactory<DerivedType>;

using RegisterUnavailableInvertParDiv = InvertParDivFactory::RegisterUnavailableInFactory;

/// Base class for parallel inversion solvers
/*!
 *
 * Inverts a matrix of the form
 *
 * A + Div_par( B * Grad_par )
 *
 * Example
 * -------
 *
 * auto inv = InvertParDiv::create();
 * inv->setCoefA(1.0);
 * inv->setCoefB(-0.1);
 *
 * Field3D result = inv->solve(rhs);
 */
class InvertParDiv {
public:
  /*!
   * Constructor. Note that this is a base class,
   * with pure virtual members, so can't be created directly.
   * To create an InvertParDiv object call the create() static function.
   */
  InvertParDiv(Options* UNUSED(opt), CELL_LOC location_in, Mesh* mesh_in = nullptr)
      : location(location_in),
        localmesh(mesh_in == nullptr ? bout::globals::mesh : mesh_in) {}
  virtual ~InvertParDiv() = default;

  /*!
   * Create an instance of InvertParDiv
   */
  static std::unique_ptr<InvertParDiv> create(Options* opt_in = nullptr,
                                              CELL_LOC location_in = CELL_CENTRE,
                                              Mesh* mesh_in = nullptr) {
    return InvertParDivFactory::getInstance().create(opt_in, location_in, mesh_in);
  }

  /*!
   * Solve the system of equations
   * Warning: Default implementation very inefficient. This converts
   * the Field2D to a Field3D then calls solve() on the 3D variable
   */
  virtual Field2D solve(const Field2D& f);

  /*!
   * Solve the system of equations
   *
   * This method must be implemented
   */
  virtual Field3D solve(const Field3D& f) = 0;

  /*!
   * Solve, given an initial guess for the solution
   * This can help if using an iterative scheme
   */
  virtual Field3D solve(const Field2D& f, const Field2D& UNUSED(start)) {
    return solve(f);
  }
  virtual Field3D solve(const Field3D& f, const Field3D& UNUSED(start)) {
    return solve(f);
  }

  /*!
   * Set the constant coefficient A
   */
  virtual void setCoefA(const Field2D& f) = 0;
  virtual void setCoefA(const Field3D& f) { setCoefA(DC(f)); }
  virtual void setCoefA(BoutReal f) {
    auto A = Field2D(f, localmesh);
    A.setLocation(location);
    setCoefA(A);
  }

  /*!
   * Set the Div_par(Grad_par()) coefficient B
   */
  virtual void setCoefB(const Field2D& f) = 0;
  virtual void setCoefB(const Field3D& f) { setCoefB(DC(f)); }
  virtual void setCoefB(BoutReal f) {
    auto B = Field2D(f, localmesh);
    B.setLocation(location);
    setCoefB(B);
  }

protected:
  CELL_LOC location;
  Mesh* localmesh; ///< Mesh object for this solver

private:
};

#endif // INV_PARDIV_H
