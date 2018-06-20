
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

void AiolosMesh::interp_to_CtoL_Field3D_x(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
      }
    }
  }
  int x = 1;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
      result_ptr[((x - 1) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
    }
  }
  x = LocalNx - 1;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x - 3) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_CtoL_Field3D_y(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
      }
    }
  }
  int y = 1;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
      result_ptr[((x)*LocalNy + y - 1) * LocalNz + z] =
          +static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
    }
  }
  y = LocalNy - 1;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y - 3) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_CtoL_Field3D_z(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        for (int z = 2; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 2 + 2 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 1 + 1 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];
        }
      }
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_x(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
      }
    }
  }
  int x = 0;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x + 3) * LocalNy + y) * LocalNz + z];
    }
  }
  x = LocalNx - 2;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
      result_ptr[((x + 1) * LocalNy + y) * LocalNz + z] =
          -static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_y(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
      }
    }
  }
  int y = 0;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y + 3) * LocalNz + z];
    }
  }
  y = LocalNy - 2;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
      result_ptr[((x)*LocalNy + y + 1) * LocalNz + z] =
          -static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_z(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2];
        }
        for (int z = 1; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2];
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 1 + 1 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)];
        }
      }
    }
  }
}
const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc,
                                       REGION region) const {
  Field3D result((AiolosMesh *)this);
  result.allocate();
  Indices i0{0, 0, 0};
  if (f.getLocation() != CELL_CENTRE) {
    // we first need to go back to centre before we can go anywhere else
    switch (f.getLocation()) {
    case CELL_XLOW:
      interp_to_LtoC_Field3D_x(&result[i0], &f[i0]);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc, region);
      break;
    case CELL_YLOW:
      interp_to_LtoC_Field3D_y(&result[i0], &f[i0]);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc, region);
      break;
    case CELL_ZLOW:
      interp_to_LtoC_Field3D_z(&result[i0], &f[i0]);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc, region);
      break;
    default:
      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                          strLocation(loc));
    }
  }
  // we are in CELL_CENTRE and need to go somewhere else ...
  switch (loc) {
  case CELL_XLOW:
    interp_to_CtoL_Field3D_x(&result[i0], &f[i0]);
    result.setLocation(CELL_XLOW);
    // return or interpolate again
    return interp_to(result, loc, region);
    break;
  case CELL_YLOW:
    interp_to_CtoL_Field3D_y(&result[i0], &f[i0]);
    result.setLocation(CELL_YLOW);
    // return or interpolate again
    return interp_to(result, loc, region);
    break;
  case CELL_ZLOW:
    interp_to_CtoL_Field3D_z(&result[i0], &f[i0]);
    result.setLocation(CELL_ZLOW);
    // return or interpolate again
    return interp_to(result, loc, region);
    break;
  default:
    throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                        strLocation(loc));
  }
}
