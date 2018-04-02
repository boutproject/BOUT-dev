
#include "aiolosmesh.hxx"
#include <interpolation.hxx>

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

static void interp_to_CtoL_Field3D_x(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x - 2) * Ny + y) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] -
            static_cast<const BoutReal>(1. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z];
      }
    }
  }
  int x = 1;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x + 0) * Ny + y) * Nz + z] =
          +static_cast<const BoutReal>(5. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(1. / 16.) * in_ptr[((x + 2) * Ny + y) * Nz + z];
      result_ptr[((x - 1) * Ny + y) * Nz + z] =
          +static_cast<const BoutReal>(35. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(21. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 2) * Ny + y) * Nz + z];
    }
  }
  x = Nx - 1;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x + 0) * Ny + y) * Nz + z] =
          +static_cast<const BoutReal>(1. / 16.) * in_ptr[((x - 3) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x - 2) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z];
    }
  }
}
static void interp_to_CtoL_Field3D_y(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y - 2) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] -
            static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z];
      }
    }
  }
  int y = 1;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x)*Ny + y + 0) * Nz + z] =
          +static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z] +
          static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y + 2) * Nz + z];
      result_ptr[((x)*Ny + y - 1) * Nz + z] =
          +static_cast<const BoutReal>(35. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] -
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] +
          static_cast<const BoutReal>(21. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 2) * Nz + z];
    }
  }
  y = Ny - 1;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x)*Ny + y + 0) * Nz + z] =
          +static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y - 3) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y - 2) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] +
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z];
    }
  }
}
static void interp_to_CtoL_Field3D_z(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1];
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1];
        }
        for (int z = 2; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1];
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];
        }
      }
    }
  }
}
static void interp_to_LtoC_Field3D_x(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z] -
            static_cast<const BoutReal>(1. / 16.) * in_ptr[((x + 2) * Ny + y) * Nz + z];
      }
    }
  }
  int x = 0;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x + 0) * Ny + y) * Nz + z] =
          +static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 2) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(1. / 16.) * in_ptr[((x + 3) * Ny + y) * Nz + z];
    }
  }
  x = Nx - 2;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x + 0) * Ny + y) * Nz + z] =
          +static_cast<const BoutReal>(1. / 16.) * in_ptr[((x - 2) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z];
      result_ptr[((x + 1) * Ny + y) * Nz + z] =
          -static_cast<const BoutReal>(5. / 16.) * in_ptr[((x - 2) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(21. / 16.) * in_ptr[((x - 1) * Ny + y) * Nz + z] -
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x + 0) * Ny + y) * Nz + z] +
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x + 1) * Ny + y) * Nz + z];
    }
  }
}
static void interp_to_LtoC_Field3D_y(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] +
            static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z] -
            static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y + 2) * Nz + z];
      }
    }
  }
  int y = 0;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x)*Ny + y + 0) * Nz + z] =
          +static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 2) * Nz + z] +
          static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y + 3) * Nz + z];
    }
  }
  y = Ny - 2;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result_ptr[((x)*Ny + y + 0) * Nz + z] =
          +static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y - 2) * Nz + z] -
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] +
          static_cast<const BoutReal>(15. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] +
          static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z];
      result_ptr[((x)*Ny + y + 1) * Nz + z] =
          -static_cast<const BoutReal>(5. / 16.) * in_ptr[((x)*Ny + y - 2) * Nz + z] +
          static_cast<const BoutReal>(21. / 16.) * in_ptr[((x)*Ny + y - 1) * Nz + z] -
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x)*Ny + y + 0) * Nz + z] +
          static_cast<const BoutReal>(35. / 16.) * in_ptr[((x)*Ny + y + 1) * Nz + z];
    }
  }
}
static void interp_to_LtoC_Field3D_z(BoutReal *__restrict__ result_ptr,
                                     const BoutReal *__restrict__ in_ptr,
                                     Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 2];
        }
        for (int z = 1; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 2];
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz];
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) * in_ptr[((x)*Ny + y) * Nz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) * in_ptr[((x)*Ny + y) * Nz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz];
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)];
        }
      }
    }
  }
}
const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc) const {
  Mesh *localmesh = f.getMesh();
  Field3D result(localmesh);
  result.allocate();
  Indices i0{0, 0, 0};
  if (f.getLocation() != CELL_CENTRE) {
    // we first need to go back to centre before we can go anywhere else
    switch (f.getLocation()) {
    case CELL_XLOW:
      interp_to_LtoC_Field3D_x(&result[i0], &f[i0], localmesh);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc);
      break;
    case CELL_YLOW:
      interp_to_LtoC_Field3D_y(&result[i0], &f[i0], localmesh);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc);
      break;
    case CELL_ZLOW:
      interp_to_LtoC_Field3D_z(&result[i0], &f[i0], localmesh);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc);
      break;
    default:
      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                          strLocation(loc));
    }
  }
  // we are in CELL_CENTRE and need to go somewhere else ...
  switch (loc) {
  case CELL_XLOW:
    interp_to_CtoL_Field3D_x(&result[i0], &f[i0], localmesh);
    result.setLocation(CELL_XLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  case CELL_YLOW:
    interp_to_CtoL_Field3D_y(&result[i0], &f[i0], localmesh);
    result.setLocation(CELL_YLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  case CELL_ZLOW:
    interp_to_CtoL_Field3D_z(&result[i0], &f[i0], localmesh);
    result.setLocation(CELL_ZLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  default:
    throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                        strLocation(loc));
  }
}
