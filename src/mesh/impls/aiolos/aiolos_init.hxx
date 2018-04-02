/********************************************************
 * BOUT++ Library - Write fluid simulations in curviilinear geometry
 * Copyright (C) 2016, 2017, 2018 David Schwörer
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
 *****************************************************************/

#pragma once

#include "aiolosmesh.hxx"

// This file is auto-generated - do not edit!
extern DIFF_METHOD default_stencil[8][3];
enum AIOLOS_DIFF_TYPE {
  AIOLOS_First = 0,
  AIOLOS_FirstStag = 1,
  AIOLOS_Second = 2,
  AIOLOS_SecondStag = 3,
  AIOLOS_Upwind = 4,
  AIOLOS_UpwindStag = 5,
  AIOLOS_Flux = 6,
  AIOLOS_FluxStag = 7,
};
