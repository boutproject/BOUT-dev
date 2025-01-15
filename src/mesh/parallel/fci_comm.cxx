/**************************************************************************
 * Communication for Flux-coordinate Independent interpolation
 *
 **************************************************************************
 * Copyright 2025 BOUT++ contributors
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

#include "fci_comm.hxx"

#include <vector>

const BoutReal& GlobalField3DAccessInstance::operator[](IndG3D ind) const {
  auto it = gfa.mapping.find(ind.ind);
  ASSERT2(it != gfa.mapping.end());
  return data[it->second];
}
