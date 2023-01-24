/************************************************************************
 * Inversion of parallel derivatives
 *
 * Inverts a matrix of the form
 *
 * (A + Div_par( B * Grad_par ) x = r
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

#include "impls/cyclic/pardiv_cyclic.hxx"
#include <bout/invert_pardiv.hxx>

Field2D InvertParDiv::solve(const Field2D& f) {
  Field3D var(f);

  var = solve(var);
  return DC(var);
}

// DO NOT REMOVE: ensures linker keeps all symbols in this TU
void InvertParDivFactory::ensureRegistered() {}
constexpr decltype(InvertParDivFactory::type_name) InvertParDivFactory::type_name;
constexpr decltype(InvertParDivFactory::section_name) InvertParDivFactory::section_name;
constexpr decltype(InvertParDivFactory::option_name) InvertParDivFactory::option_name;
constexpr decltype(InvertParDivFactory::default_type) InvertParDivFactory::default_type;
