
#include "aiolosmesh.hxx"
#include <interpolation.hxx>
#include <output.hxx>

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

const BoutReal WENO_SMALL = 1.0e-8;

void AiolosMesh::DDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            0.5 * (in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            0.5 * (in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 *
              (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          0.5 * (in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          0.5 * (in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (8. * in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             8. * in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
             in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (8. * in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             8. * in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
             in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz +
                           +((z - 1 + 1 * LocalNz) % LocalNz)] +
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 2 + 2 * LocalNz) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (8. * in_ptr[(x + 1) * LocalNy + y] - 8. * in_ptr[(x - 1) * LocalNy + y] +
           in_ptr[(x - 2) * LocalNy + y] - in_ptr[(x + 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (8. * in_ptr[(x)*LocalNy + y + 1] - 8. * in_ptr[(x)*LocalNy + y - 1] +
           in_ptr[(x)*LocalNy + y - 2] - in_ptr[(x)*LocalNy + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_CWENO2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal isl, isr, isc;
        BoutReal al, ar, ac, sa;
        BoutReal dl, dr, dc;
        dc = 0.5 * (in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                    in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
        dl = in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x - 1) * LocalNy + y) * LocalNz + z];
        dr = in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        isl = SQ(dl);
        isr = SQ(dr);
        isc = (13. / 3.) * SQ(in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                              2. * in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) +
              0.25 * SQ(in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                        in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
        al = 0.25 / SQ(WENO_SMALL + isl);
        ar = 0.25 / SQ(WENO_SMALL + isr);
        ac = 0.5 / SQ(WENO_SMALL + isc);
        sa = al + ar + ac;

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (al * dl + ar * dr + ac * dc) / sa;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_CWENO2_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_CWENO2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_CWENO2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal isl, isr, isc;
        BoutReal al, ar, ac, sa;
        BoutReal dl, dr, dc;
        dc = 0.5 * (in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                    in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
        dl = in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y - 1) * LocalNz + z];
        dr = in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        isl = SQ(dl);
        isr = SQ(dr);
        isc = (13. / 3.) * SQ(in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                              2. * in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) +
              0.25 * SQ(in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                        in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
        al = 0.25 / SQ(WENO_SMALL + isl);
        ar = 0.25 / SQ(WENO_SMALL + isr);
        ac = 0.5 / SQ(WENO_SMALL + isc);
        sa = al + ar + ac;

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (al * dl + ar * dr + ac * dc) / sa;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_CWENO2_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_CWENO2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_CWENO2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
          dl = in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz];
          dr = in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) +
                0.25 * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                          in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (al * dl + ar * dr + ac * dc) / sa;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          dl = in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1];
          dr = in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) +
                0.25 * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                          in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (al * dl + ar * dr + ac * dc) / sa;
        }
        {
          int z = LocalNz - 1;
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          dl = in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1];
          dr = in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                                2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) +
                0.25 * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                          in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (al * dl + ar * dr + ac * dc) / sa;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 *
               (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)]);
          dl = in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)];
          dr = in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) *
                    SQ(in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                       2. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                       in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)]) +
                0.25 * SQ(in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                          in_ptr[((x)*LocalNy + y) * LocalNz +
                                 +((z - 1 + 1 * LocalNz) % LocalNz)]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (al * dl + ar * dr + ac * dc) / sa;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_CWENO2_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_CWENO2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_CWENO2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal isl, isr, isc;
      BoutReal al, ar, ac, sa;
      BoutReal dl, dr, dc;
      dc = 0.5 * (in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y]);
      dl = in_ptr[(x + 0) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y];
      dr = in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x + 0) * LocalNy + y];
      isl = SQ(dl);
      isr = SQ(dr);
      isc = (13. / 3.) *
                SQ(in_ptr[(x + 1) * LocalNy + y] - 2. * in_ptr[(x + 0) * LocalNy + y] +
                   in_ptr[(x - 1) * LocalNy + y]) +
            0.25 * SQ(in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y]);
      al = 0.25 / SQ(WENO_SMALL + isl);
      ar = 0.25 / SQ(WENO_SMALL + isr);
      ac = 0.5 / SQ(WENO_SMALL + isc);
      sa = al + ar + ac;

      result_ptr[(x + 0) * LocalNy + y] = (al * dl + ar * dr + ac * dc) / sa;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_CWENO2_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_CWENO2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_CWENO2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_CWENO2_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_CWENO2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal isl, isr, isc;
      BoutReal al, ar, ac, sa;
      BoutReal dl, dr, dc;
      dc = 0.5 * (in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y - 1]);
      dl = in_ptr[(x)*LocalNy + y + 0] - in_ptr[(x)*LocalNy + y - 1];
      dr = in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y + 0];
      isl = SQ(dl);
      isr = SQ(dr);
      isc =
          (13. / 3.) * SQ(in_ptr[(x)*LocalNy + y + 1] - 2. * in_ptr[(x)*LocalNy + y + 0] +
                          in_ptr[(x)*LocalNy + y - 1]) +
          0.25 * SQ(in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y - 1]);
      al = 0.25 / SQ(WENO_SMALL + isl);
      ar = 0.25 / SQ(WENO_SMALL + isr);
      ac = 0.5 / SQ(WENO_SMALL + isc);
      sa = al + ar + ac;

      result_ptr[(x)*LocalNy + y + 0] = (al * dl + ar * dr + ac * dc) / sa;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_CWENO2_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_CWENO2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_CWENO2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_CWENO2_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_S2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (8. * in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                            8. * in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                            in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
                            in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
                           12.;
        result_ += SIGN(in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) *
                   (in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] -
                    4. * in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
                    6. * in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                    4. * in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                    in_ptr[((x - 2) * LocalNy + y) * LocalNz + z]) /
                   12.;

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_S2_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_S2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_S2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                            8. * in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                            in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
                            in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
                           12.;
        result_ += SIGN(in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) *
                   (in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] -
                    4. * in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
                    6. * in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                    4. * in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                    in_ptr[((x)*LocalNy + y - 2) * LocalNz + z]) /
                   12.;

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_S2_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_S2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_S2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                              8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                              in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                              in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                     (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                      6. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) /
                     12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                              8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                              in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                              in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                     (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                      6. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) /
                     12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {
          BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                              8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                              in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                              in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                     (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                      6. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
                     12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 2;
          BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                              8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                              in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                              in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                     (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                      6. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
                     12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_ = (8. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                              8. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                              in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                              in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                     (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
                      6. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      4. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
                     12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_ =
              (8. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               8. * in_ptr[((x)*LocalNy + y) * LocalNz +
                           +((z - 1 + 1 * LocalNz) % LocalNz)] +
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 2 + 2 * LocalNz) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
          result_ +=
              SIGN(in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) *
              (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] -
               4. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
               6. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
               4. * in_ptr[((x)*LocalNy + y) * LocalNz +
                           +((z - 1 + 1 * LocalNz) % LocalNz)] +
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 2 + 2 * LocalNz) % LocalNz)]) /
              12.;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_S2_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_S2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_S2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_ =
          (8. * in_ptr[(x + 1) * LocalNy + y] - 8. * in_ptr[(x - 1) * LocalNy + y] +
           in_ptr[(x - 2) * LocalNy + y] - in_ptr[(x + 2) * LocalNy + y]) /
          12.;
      result_ += SIGN(in_ptr[(x + 0) * LocalNy + y]) *
                 (in_ptr[(x + 2) * LocalNy + y] - 4. * in_ptr[(x + 1) * LocalNy + y] +
                  6. * in_ptr[(x + 0) * LocalNy + y] -
                  4. * in_ptr[(x - 1) * LocalNy + y] + in_ptr[(x - 2) * LocalNy + y]) /
                 12.;

      result_ptr[(x + 0) * LocalNy + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_S2_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_S2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_S2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_S2_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_S2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      BoutReal result_ =
          (8. * in_ptr[(x)*LocalNy + y + 1] - 8. * in_ptr[(x)*LocalNy + y - 1] +
           in_ptr[(x)*LocalNy + y - 2] - in_ptr[(x)*LocalNy + y + 2]) /
          12.;
      result_ += SIGN(in_ptr[(x)*LocalNy + y + 0]) *
                 (in_ptr[(x)*LocalNy + y + 2] - 4. * in_ptr[(x)*LocalNy + y + 1] +
                  6. * in_ptr[(x)*LocalNy + y + 0] - 4. * in_ptr[(x)*LocalNy + y - 1] +
                  in_ptr[(x)*LocalNy + y - 2]) /
                 12.;

      result_ptr[(x)*LocalNy + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_S2_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_S2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_S2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_S2_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
            in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
            2. * in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
            in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
            2. * in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
              in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
              2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
              in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
              2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
              in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
              2. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)] -
              2. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] = in_ptr[(x + 1) * LocalNy + y] +
                                          in_ptr[(x - 1) * LocalNy + y] -
                                          2. * in_ptr[(x + 0) * LocalNy + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] = in_ptr[(x)*LocalNy + y + 1] +
                                        in_ptr[(x)*LocalNy + y - 1] -
                                        2. * in_ptr[(x)*LocalNy + y + 0];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (-in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] +
             16. * in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             30. * in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
             16. * in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x - 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C4_x_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_x_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (-in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] +
             16. * in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             30. * in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
             16. * in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y - 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C4_y_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_y_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (-in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               30. * in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
               16. * in_ptr[((x)*LocalNy + y) * LocalNz +
                            +((z - 1 + 1 * LocalNz) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 2 + 2 * LocalNz) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C4_z_norm(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_z_Field3D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (-in_ptr[(x + 2) * LocalNy + y] + 16. * in_ptr[(x + 1) * LocalNy + y] -
           30. * in_ptr[(x + 0) * LocalNy + y] + 16. * in_ptr[(x - 1) * LocalNy + y] -
           in_ptr[(x - 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C4_x_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C4_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C4_x_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (-in_ptr[(x)*LocalNy + y + 2] + 16. * in_ptr[(x)*LocalNy + y + 1] -
           30. * in_ptr[(x)*LocalNy + y + 0] + 16. * in_ptr[(x)*LocalNy + y - 1] -
           in_ptr[(x)*LocalNy + y - 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C4_y_norm(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C4_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C4_y_Field2D_norm(result_ptr, in_ptr);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] * 0.5 *
            (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] * 0.5 *
            (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * 0.5 *
              (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * 0.5 *
              (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * 0.5 *
              (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] * 0.5 *
              (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 1 + 1 * LocalNz) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          v_in_ptr[(x + 0) * LocalNy + y] * 0.5 *
          (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          v_in_ptr[(x)*LocalNy + y + 0] * 0.5 *
          (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
            (8. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             8. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
             f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
             f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
            (8. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             8. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
             f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
             f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
              (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               8. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] +
               f_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 2 + 2 * LocalNz) % LocalNz)] -
               f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          v_in_ptr[(x + 0) * LocalNy + y] *
          (8. * f_in_ptr[(x + 1) * LocalNy + y] - 8. * f_in_ptr[(x - 1) * LocalNy + y] +
           f_in_ptr[(x - 2) * LocalNy + y] - f_in_ptr[(x + 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          v_in_ptr[(x)*LocalNy + y + 0] *
          (8. * f_in_ptr[(x)*LocalNy + y + 1] - 8. * f_in_ptr[(x)*LocalNy + y - 1] +
           f_in_ptr[(x)*LocalNy + y - 2] - f_in_ptr[(x)*LocalNy + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                       f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z])
                : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                       f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                       f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z])
                : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                       f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                  +((z - 1 + 1 * LocalNz) % LocalNz)])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          v_in_ptr[(x + 0) * LocalNy + y] >= 0.0
              ? v_in_ptr[(x + 0) * LocalNy + y] *
                    (f_in_ptr[(x + 0) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y])
              : v_in_ptr[(x + 0) * LocalNy + y] *
                    (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x + 0) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          v_in_ptr[(x)*LocalNy + y + 0] >= 0.0
              ? v_in_ptr[(x)*LocalNy + y + 0] *
                    (f_in_ptr[(x)*LocalNy + y + 0] - f_in_ptr[(x)*LocalNy + y - 1])
              : v_in_ptr[(x)*LocalNy + y + 0] *
                    (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y + 0]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                       2.0 * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                       0.5 * f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z])
                : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (-0.5 * f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] +
                       2.0 * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                       1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                       2.0 * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                       0.5 * f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z])
                : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (-0.5 * f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] +
                       2.0 * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                       1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z + 0) % LocalNz)] -
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z - 1 + 1 * LocalNz) % LocalNz)] +
                         0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z - 2 + 2 * LocalNz) % LocalNz)])
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (-0.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                         +((z + 2) % LocalNz)] +
                         2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z + 1) % LocalNz)] -
                         1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z + 0) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          v_in_ptr[(x + 0) * LocalNy + y] >= 0.0
              ? v_in_ptr[(x + 0) * LocalNy + y] * (1.5 * f_in_ptr[(x + 0) * LocalNy + y] -
                                                   2.0 * f_in_ptr[(x - 1) * LocalNy + y] +
                                                   0.5 * f_in_ptr[(x - 2) * LocalNy + y])
              : v_in_ptr[(x + 0) * LocalNy + y] *
                    (-0.5 * f_in_ptr[(x + 2) * LocalNy + y] +
                     2.0 * f_in_ptr[(x + 1) * LocalNy + y] -
                     1.5 * f_in_ptr[(x + 0) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          v_in_ptr[(x)*LocalNy + y + 0] >= 0.0
              ? v_in_ptr[(x)*LocalNy + y + 0] * (1.5 * f_in_ptr[(x)*LocalNy + y + 0] -
                                                 2.0 * f_in_ptr[(x)*LocalNy + y - 1] +
                                                 0.5 * f_in_ptr[(x)*LocalNy + y - 2])
              : v_in_ptr[(x)*LocalNy + y + 0] * (-0.5 * f_in_ptr[(x)*LocalNy + y + 2] +
                                                 2.0 * f_in_ptr[(x)*LocalNy + y + 1] -
                                                 1.5 * f_in_ptr[(x)*LocalNy + y + 0]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U3_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (4. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                       12. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                       2. * f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
                       6. * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) /
                      12.
                : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                      (-4. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                       12. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                       2. * f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] -
                       6. * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) /
                      12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U3_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U3_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0.0
                ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (4. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                       12. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                       2. * f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
                       6. * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) /
                      12.
                : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                      (-4. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                       12. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                       2. * f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] -
                       6. * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) /
                      12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U3_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U3_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
                        12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0.0
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (4. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                       +((z + 1) % LocalNz)] -
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z - 1 + 1 * LocalNz) % LocalNz)] +
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                       +((z - 2 + 2 * LocalNz) % LocalNz)] +
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                       +((z + 0) % LocalNz)]) /
                        12.
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        (-4. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z - 1 + 1 * LocalNz) % LocalNz)] +
                         12. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                        +((z + 1) % LocalNz)] -
                         2. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                       +((z + 2) % LocalNz)] -
                         6. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                       +((z + 0) % LocalNz)]) /
                        12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U3_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U3_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          v_in_ptr[(x + 0) * LocalNy + y] >= 0.0
              ? v_in_ptr[(x + 0) * LocalNy + y] * (4. * f_in_ptr[(x + 1) * LocalNy + y] -
                                                   12. * f_in_ptr[(x - 1) * LocalNy + y] +
                                                   2. * f_in_ptr[(x - 2) * LocalNy + y] +
                                                   6. * f_in_ptr[(x + 0) * LocalNy + y]) /
                    12.
              : v_in_ptr[(x + 0) * LocalNy + y] * (-4. * f_in_ptr[(x - 1) * LocalNy + y] +
                                                   12. * f_in_ptr[(x + 1) * LocalNy + y] -
                                                   2. * f_in_ptr[(x + 2) * LocalNy + y] -
                                                   6. * f_in_ptr[(x + 0) * LocalNy + y]) /
                    12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U3_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U3_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U3_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U3_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          v_in_ptr[(x)*LocalNy + y + 0] >= 0.0
              ? v_in_ptr[(x)*LocalNy + y + 0] * (4. * f_in_ptr[(x)*LocalNy + y + 1] -
                                                 12. * f_in_ptr[(x)*LocalNy + y - 1] +
                                                 2. * f_in_ptr[(x)*LocalNy + y - 2] +
                                                 6. * f_in_ptr[(x)*LocalNy + y + 0]) /
                    12.
              : v_in_ptr[(x)*LocalNy + y + 0] * (-4. * f_in_ptr[(x)*LocalNy + y - 1] +
                                                 12. * f_in_ptr[(x)*LocalNy + y + 1] -
                                                 2. * f_in_ptr[(x)*LocalNy + y + 2] -
                                                 6. * f_in_ptr[(x)*LocalNy + y + 0]) /
                    12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U3_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U3_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U3_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_WENO3_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal deriv, w, r;
        if (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] > 0.0) {
          r = (WENO_SMALL + SQ(f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                               2.0 * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                               f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                               2.0 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                               f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                         f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) -
                  0.5 * w * (-f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
                             3. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
                             3. * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                             f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z]);
        } else {
          r = (WENO_SMALL + SQ(f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] -
                               2.0 * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
                               f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                               2.0 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                               f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                         f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) -
                  0.5 * w * (-f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                             3. * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                             3. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
                             f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]);
        }

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] * deriv;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_WENO3_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_WENO3_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal deriv, w, r;
        if (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] > 0.0) {
          r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                               2.0 * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                               f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                               2.0 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                               f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                         f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) -
                  0.5 * w * (-f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
                             3. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
                             3. * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                             f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z]);
        } else {
          r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] -
                               2.0 * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
                               f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                               2.0 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                               f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                         f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) -
                  0.5 * w * (-f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                             3. * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                             3. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
                             f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]);
        }

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] * deriv;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_WENO3_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_WENO3_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0.0) {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])) /
                (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv =
                0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) -
                0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
                           3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
                           3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0])) /
                (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * deriv;
        }
        {
          int z = 1;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0.0) {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * deriv;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * deriv;
        }
        {
          int z = LocalNz - 2;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * deriv;
        }
        {
          int z = LocalNz - 1;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                               3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz]);
          } else {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                                 2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                                 f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv =
                0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
                0.5 * w * (-f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                           3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                           3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
                           f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] * deriv;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] > 0.0) {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                   +((z - 1 + 1 * LocalNz) % LocalNz)] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 2 + 2 * LocalNz) % LocalNz)])) /
                (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv =
                0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)]) -
                0.5 * w *
                    (-f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 2 + 2 * LocalNz) % LocalNz)] +
                     3. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                   +((z - 1 + 1 * LocalNz) % LocalNz)] -
                     3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)]);
          } else {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)])) /
                (WENO_SMALL +
                 SQ(f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                    2.0 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                    f_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv =
                0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)]) -
                0.5 * w *
                    (-f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 1 + 1 * LocalNz) % LocalNz)] +
                     3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                     3. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] * deriv;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_WENO3_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_WENO3_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal deriv, w, r;
      if (v_in_ptr[(x + 0) * LocalNy + y] > 0.0) {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x + 0) * LocalNy + y] - 2.0 * f_in_ptr[(x - 1) * LocalNy + y] +
                f_in_ptr[(x - 2) * LocalNy + y])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x + 1) * LocalNy + y] - 2.0 * f_in_ptr[(x + 0) * LocalNy + y] +
                f_in_ptr[(x - 1) * LocalNy + y]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv =
            0.5 * (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]) -
            0.5 * w *
                (-f_in_ptr[(x - 2) * LocalNy + y] + 3. * f_in_ptr[(x - 1) * LocalNy + y] -
                 3. * f_in_ptr[(x + 0) * LocalNy + y] + f_in_ptr[(x + 1) * LocalNy + y]);
      } else {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x + 2) * LocalNy + y] - 2.0 * f_in_ptr[(x + 1) * LocalNy + y] +
                f_in_ptr[(x + 0) * LocalNy + y])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x + 1) * LocalNy + y] - 2.0 * f_in_ptr[(x + 0) * LocalNy + y] +
                f_in_ptr[(x - 1) * LocalNy + y]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv =
            0.5 * (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]) -
            0.5 * w *
                (-f_in_ptr[(x - 1) * LocalNy + y] + 3. * f_in_ptr[(x + 0) * LocalNy + y] -
                 3. * f_in_ptr[(x + 1) * LocalNy + y] + f_in_ptr[(x + 2) * LocalNy + y]);
      }

      result_ptr[(x + 0) * LocalNy + y] = v_in_ptr[(x + 0) * LocalNy + y] * deriv;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_WENO3_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_WENO3_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_WENO3_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_WENO3_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      BoutReal deriv, w, r;
      if (v_in_ptr[(x)*LocalNy + y + 0] > 0.0) {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x)*LocalNy + y + 0] - 2.0 * f_in_ptr[(x)*LocalNy + y - 1] +
                f_in_ptr[(x)*LocalNy + y - 2])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x)*LocalNy + y + 1] - 2.0 * f_in_ptr[(x)*LocalNy + y + 0] +
                f_in_ptr[(x)*LocalNy + y - 1]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]) -
                0.5 * w *
                    (-f_in_ptr[(x)*LocalNy + y - 2] + 3. * f_in_ptr[(x)*LocalNy + y - 1] -
                     3. * f_in_ptr[(x)*LocalNy + y + 0] + f_in_ptr[(x)*LocalNy + y + 1]);
      } else {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x)*LocalNy + y + 2] - 2.0 * f_in_ptr[(x)*LocalNy + y + 1] +
                f_in_ptr[(x)*LocalNy + y + 0])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x)*LocalNy + y + 1] - 2.0 * f_in_ptr[(x)*LocalNy + y + 0] +
                f_in_ptr[(x)*LocalNy + y - 1]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]) -
                0.5 * w *
                    (-f_in_ptr[(x)*LocalNy + y - 1] + 3. * f_in_ptr[(x)*LocalNy + y + 0] -
                     3. * f_in_ptr[(x)*LocalNy + y + 1] + f_in_ptr[(x)*LocalNy + y + 2]);
      }

      result_ptr[(x)*LocalNy + y + 0] = v_in_ptr[(x)*LocalNy + y + 0] * deriv;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_WENO3_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_WENO3_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_WENO3_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal vs = 0.5 * (v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                             v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);
        BoutReal result_ = (vs >= 0.0)
                               ? vs * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]
                               : vs * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        vs = 0.5 * (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                    v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z]);
        result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]
                               : vs * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal vs = 0.5 * (v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                             v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);
        BoutReal result_ = (vs >= 0.0)
                               ? vs * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]
                               : vs * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        vs = 0.5 * (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                    v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z]);
        result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]
                               : vs * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          BoutReal result_ =
              (vs >= 0.0) ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]
                          : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                                 : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          BoutReal result_ = (vs >= 0.0)
                                 ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                                 : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          BoutReal result_ = (vs >= 0.0)
                                 ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz]);
          result_ -= (vs >= 0.0)
                         ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal vs =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]);
          BoutReal result_ =
              (vs >= 0.0)
                  ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                  +((z - 1 + 1 * LocalNz) % LocalNz)]
                  : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          vs = 0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)]);
          result_ -=
              (vs >= 0.0)
                  ? vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]
                  : vs * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal vs =
          0.5 * (v_in_ptr[(x - 1) * LocalNy + y] + v_in_ptr[(x + 0) * LocalNy + y]);
      BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[(x - 1) * LocalNy + y]
                                     : vs * f_in_ptr[(x + 0) * LocalNy + y];
      vs = 0.5 * (v_in_ptr[(x + 0) * LocalNy + y] + v_in_ptr[(x + 1) * LocalNy + y]);
      result_ -= (vs >= 0.0) ? vs * f_in_ptr[(x + 0) * LocalNy + y]
                             : vs * f_in_ptr[(x + 1) * LocalNy + y];

      result_ptr[(x + 0) * LocalNy + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal vs = 0.5 * (v_in_ptr[(x)*LocalNy + y - 1] + v_in_ptr[(x)*LocalNy + y + 0]);
      BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[(x)*LocalNy + y - 1]
                                     : vs * f_in_ptr[(x)*LocalNy + y + 0];
      vs = 0.5 * (v_in_ptr[(x)*LocalNy + y + 0] + v_in_ptr[(x)*LocalNy + y + 1]);
      result_ -= (vs >= 0.0) ? vs * f_in_ptr[(x)*LocalNy + y + 0]
                             : vs * f_in_ptr[(x)*LocalNy + y + 1];

      result_ptr[(x)*LocalNy + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            0.5 * (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                       f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                       f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C2_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            0.5 * (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                       f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                       f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C2_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                     v_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)] *
                         f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                  +((z - 1 + 1 * LocalNz) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C2_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          0.5 * (v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y] -
                 v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_C2_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          0.5 * (v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 1] -
                 v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_C2_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (8. * v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                 f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             8. * v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                 f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
             v_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] *
                 f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
             v_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] *
                 f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C4_x_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (8. * v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                 f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             8. * v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                 f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
             v_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] *
                 f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
             v_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] *
                 f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C4_y_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (8. * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               8. * v_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz +
                            +((z - 1 + 1 * LocalNz) % LocalNz)] +
               v_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 2 + 2 * LocalNz) % LocalNz)] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz +
                            +((z - 2 + 2 * LocalNz) % LocalNz)] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] *
                   f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_C4_z_norm(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (8. * v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y] -
           8. * v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y] +
           v_in_ptr[(x - 2) * LocalNy + y] * f_in_ptr[(x - 2) * LocalNy + y] -
           v_in_ptr[(x + 2) * LocalNy + y] * f_in_ptr[(x + 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_C4_x_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::FDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C4_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (8. * v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 1] -
           8. * v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y - 1] +
           v_in_ptr[(x)*LocalNy + y - 2] * f_in_ptr[(x)*LocalNy + y - 2] -
           v_in_ptr[(x)*LocalNy + y + 2] * f_in_ptr[(x)*LocalNy + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_C4_y_norm(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::FDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C4_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx + 0; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
            in_ptr[((x - 1) * LocalNy + y) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_x_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_x_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy + 0; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
            in_ptr[((x)*LocalNy + y - 1) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_y_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_y_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz];
        }
        for (int z = 1; z < LocalNz - 0; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              in_ptr[((x)*LocalNy + y) * LocalNz + z - 1];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_z_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_z_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx + 0; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          in_ptr[(x + 0) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_stag_x_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_x_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy + 0; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          in_ptr[(x)*LocalNy + y + 0] - in_ptr[(x)*LocalNy + y - 1];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_stag_y_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_y_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 0; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
            in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_x_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_x_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 0; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
            in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_y_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_y_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
              in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
              in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C2_stag_z_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_z_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 0; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x + 0) * LocalNy + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_stag_x_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_x_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 0; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y + 0];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C2_stag_y_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_y_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                    in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) -
             (in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_x_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_x_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                    in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) -
             (in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_y_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_y_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])) /
              24.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz])) /
              24.;
        }
        for (int z = 2; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])) /
              24.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 2])) /
              24.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                      in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                in_ptr[((x)*LocalNy + y) * LocalNz +
                       +((z - 2 + 2 * LocalNz) % LocalNz)])) /
              24.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_z_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_z_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (27. * (in_ptr[(x + 0) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y]) -
           (in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x - 2) * LocalNy + y])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_stag_x_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_x_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (27. * (in_ptr[(x)*LocalNy + y + 0] - in_ptr[(x)*LocalNy + y - 1]) -
           (in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y - 2])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_stag_y_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_y_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                    in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) -
             (in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] -
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_x_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_x_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                    in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) -
             (in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] -
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_y_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_y_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz])) /
              24.;
        }
        for (int z = 1; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1])) /
              24.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1])) /
              24.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] -
                in_ptr[((x)*LocalNy + y) * LocalNz + z - 1])) /
              24.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (27. * (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                      in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) -
               (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] -
                in_ptr[((x)*LocalNy + y) * LocalNz +
                       +((z - 1 + 1 * LocalNz) % LocalNz)])) /
              24.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::DDX_C4_stag_z_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_z_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (27. * (in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x + 0) * LocalNy + y]) -
           (in_ptr[(x + 2) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_stag_x_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_x_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::DDX_C4_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (27. * (in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y + 0]) -
           (in_ptr[(x)*LocalNy + y + 2] - in_ptr[(x)*LocalNy + y - 1])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::DDX_C4_stag_y_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method DDX_C4_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_y_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
             in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_x_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_x_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
             in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_y_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_y_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) /
              2.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) /
              2.;
        }
        for (int z = 2; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) /
              2.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) /
              2.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 2 + 2 * LocalNz) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)]) /
              2.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_z_CtoL(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_z_Field3D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (in_ptr[(x + 1) * LocalNy + y] + in_ptr[(x - 2) * LocalNy + y] -
           in_ptr[(x + 0) * LocalNy + y] - in_ptr[(x - 1) * LocalNy + y]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_stag_x_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_x_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (in_ptr[(x)*LocalNy + y + 1] + in_ptr[(x)*LocalNy + y - 2] -
           in_ptr[(x)*LocalNy + y + 0] - in_ptr[(x)*LocalNy + y - 1]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_stag_y_CtoL(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_y_Field2D_CtoL(result_ptr, in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] +
             in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
             in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_x_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_x_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] +
             in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
             in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_y_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_y_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
              2.;
        }
        for (int z = 1; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
              2.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
              2.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz] +
               in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
               in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) /
              2.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)] +
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
               in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) /
              2.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::D2DX2_C2_stag_z_LtoC(const Field3D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_z_Field3D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (in_ptr[(x + 2) * LocalNy + y] + in_ptr[(x - 1) * LocalNy + y] -
           in_ptr[(x + 1) * LocalNy + y] - in_ptr[(x + 0) * LocalNy + y]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_stag_x_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_x_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::D2DX2_C2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                              const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (in_ptr[(x)*LocalNy + y + 2] + in_ptr[(x)*LocalNy + y - 1] -
           in_ptr[(x)*LocalNy + y + 1] - in_ptr[(x)*LocalNy + y + 0]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::D2DX2_C2_stag_y_LtoC(const Field2D &in) const {
  ASSERT1(this == in.getMesh());
  output_debug.write("Using method D2DX2_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_y_Field2D_LtoC(result_ptr, in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_x_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]
                               : v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        result_ -= (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]
                       : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                   (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                    v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_y_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]
                               : v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        result_ -= (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]
                       : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                   (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                    v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_z_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 1 + 1 * LocalNz) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                 +((z - 1 + 1 * LocalNz) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          result_ -=
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 1 + 1 * LocalNz) % LocalNz)]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_x_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x - 1) * LocalNy + y] >= 0)
              ? v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y]
              : v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y];
      result_ -= (v_in_ptr[(x + 0) * LocalNy + y] >= 0)
                     ? v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y]
                     : v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y];
      result_ *= -1;
      result_ -= f_in_ptr[(x + 0) * LocalNy + y] *
                 (v_in_ptr[(x + 0) * LocalNy + y] - v_in_ptr[(x - 1) * LocalNy + y]);

      result_ptr[(x + 0) * LocalNy + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_y_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x)*LocalNy + y - 1] >= 0)
              ? v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y - 1]
              : v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y + 0];
      result_ -= (v_in_ptr[(x)*LocalNy + y + 0] >= 0)
                     ? v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 0]
                     : v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 1];
      result_ *= -1;
      result_ -= f_in_ptr[(x)*LocalNy + y + 0] *
                 (v_in_ptr[(x)*LocalNy + y + 0] - v_in_ptr[(x)*LocalNy + y - 1]);

      result_ptr[(x)*LocalNy + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_x_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]
                               : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        result_ -= (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]
                       : v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                   (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                    v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_y_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]
                               : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        result_ -= (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]
                       : v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                   (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                    v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_z_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                 +((z - 1 + 1 * LocalNz) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          result_ -=
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                     (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                      v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]);

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U1_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_x_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x + 0) * LocalNy + y] >= 0)
              ? v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y]
              : v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y];
      result_ -= (v_in_ptr[(x + 1) * LocalNy + y] >= 0)
                     ? v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y]
                     : v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y];
      result_ *= -1;
      result_ -= f_in_ptr[(x + 0) * LocalNy + y] *
                 (v_in_ptr[(x + 1) * LocalNy + y] - v_in_ptr[(x + 0) * LocalNy + y]);

      result_ptr[(x + 0) * LocalNy + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U1_stag_y_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x)*LocalNy + y + 0] >= 0)
              ? v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y - 1]
              : v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 0];
      result_ -= (v_in_ptr[(x)*LocalNy + y + 1] >= 0)
                     ? v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 0]
                     : v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 1];
      result_ *= -1;
      result_ -= f_in_ptr[(x)*LocalNy + y + 0] *
                 (v_in_ptr[(x)*LocalNy + y + 1] - v_in_ptr[(x)*LocalNy + y + 0]);

      result_ptr[(x)*LocalNy + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U1_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_x_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] > 0 &&
            v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
                     .5 * v_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z]) *
                    (.5 * f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
                     2. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                     1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);
        } else if (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] < 0 &&
                   v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                     .5 * v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z]) *
                    (-1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                     2. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                     .5 * f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                           v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) *
                    (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
        }

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_y_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] > 0 &&
            v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
                     .5 * v_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z]) *
                    (.5 * f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
                     2. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                     1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);
        } else if (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] < 0 &&
                   v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                     .5 * v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z]) *
                    (-1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                     2. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                     .5 * f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                           v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) *
                    (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
        }

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_z_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 2;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z - 1 + 1 * LocalNz) % LocalNz)] >
                  0) {
            result_ =
                (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)] -
                 .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 2 + 2 * LocalNz) % LocalNz)]) *
                (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 2 + 2 * LocalNz) % LocalNz)] -
                 2. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 1 + 1 * LocalNz) % LocalNz)] +
                 1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)] < 0) {
            result_ =
                (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                 .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)]) *
                (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                 2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                 .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]);
          } else {
            result_ = .25 *
                      (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                       v_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_x_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x + 0) * LocalNy + y] > 0 && v_in_ptr[(x - 1) * LocalNy + y] > 0) {
        result_ =
            (1.5 * v_in_ptr[(x - 1) * LocalNy + y] -
             .5 * v_in_ptr[(x - 2) * LocalNy + y]) *
            (.5 * f_in_ptr[(x - 2) * LocalNy + y] - 2. * f_in_ptr[(x - 1) * LocalNy + y] +
             1.5 * f_in_ptr[(x + 0) * LocalNy + y]);
      } else if (v_in_ptr[(x + 0) * LocalNy + y] < 0 &&
                 v_in_ptr[(x - 1) * LocalNy + y] < 0) {
        result_ =
            (1.5 * v_in_ptr[(x + 0) * LocalNy + y] -
             .5 * v_in_ptr[(x + 1) * LocalNy + y]) *
            (-1.5 * f_in_ptr[(x + 0) * LocalNy + y] +
             2. * f_in_ptr[(x + 1) * LocalNy + y] - .5 * f_in_ptr[(x + 2) * LocalNy + y]);
      } else {
        result_ = .25 *
                  (v_in_ptr[(x + 0) * LocalNy + y] + v_in_ptr[(x - 1) * LocalNy + y]) *
                  (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]);
      }

      result_ptr[(x + 0) * LocalNy + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_y_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x)*LocalNy + y + 0] > 0 && v_in_ptr[(x)*LocalNy + y - 1] > 0) {
        result_ =
            (1.5 * v_in_ptr[(x)*LocalNy + y - 1] - .5 * v_in_ptr[(x)*LocalNy + y - 2]) *
            (.5 * f_in_ptr[(x)*LocalNy + y - 2] - 2. * f_in_ptr[(x)*LocalNy + y - 1] +
             1.5 * f_in_ptr[(x)*LocalNy + y + 0]);
      } else if (v_in_ptr[(x)*LocalNy + y + 0] < 0 && v_in_ptr[(x)*LocalNy + y - 1] < 0) {
        result_ =
            (1.5 * v_in_ptr[(x)*LocalNy + y + 0] - .5 * v_in_ptr[(x)*LocalNy + y + 1]) *
            (-1.5 * f_in_ptr[(x)*LocalNy + y + 0] + 2. * f_in_ptr[(x)*LocalNy + y + 1] -
             .5 * f_in_ptr[(x)*LocalNy + y + 2]);
      } else {
        result_ = .25 * (v_in_ptr[(x)*LocalNy + y + 0] + v_in_ptr[(x)*LocalNy + y - 1]) *
                  (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]);
      }

      result_ptr[(x)*LocalNy + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_x_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] > 0 &&
            v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
                     .5 * v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) *
                    (.5 * f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
                     2. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                     1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]);
        } else if (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] < 0 &&
                   v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                     .5 * v_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) *
                    (-1.5 * f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                     2. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                     .5 * f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
                           v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) *
                    (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
        }

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_y_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] > 0 &&
            v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
                     .5 * v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) *
                    (.5 * f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
                     2. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                     1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]);
        } else if (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] < 0 &&
                   v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                     .5 * v_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) *
                    (-1.5 * f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                     2. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                     .5 * f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
                           v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) *
                    (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
        }

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_z_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 2;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
                      (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                       1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) *
                      (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                       2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
                             v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] > 0 &&
              v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] > 0) {
            result_ =
                (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
                 .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 1 + 1 * LocalNz) % LocalNz)]) *
                (.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 2 + 2 * LocalNz) % LocalNz)] -
                 2. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                               +((z - 1 + 1 * LocalNz) % LocalNz)] +
                 1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]);
          } else if (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] < 0 &&
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] < 0) {
            result_ =
                (1.5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                 .5 * v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) *
                (-1.5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                 2. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                 .5 * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]);
          } else {
            result_ = .25 *
                      (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
                       v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) *
                      (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                       f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                +((z - 1 + 1 * LocalNz) % LocalNz)]);
          }

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_U2_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_x_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x + 1) * LocalNy + y] > 0 && v_in_ptr[(x + 0) * LocalNy + y] > 0) {
        result_ =
            (1.5 * v_in_ptr[(x + 0) * LocalNy + y] -
             .5 * v_in_ptr[(x - 1) * LocalNy + y]) *
            (.5 * f_in_ptr[(x - 2) * LocalNy + y] - 2. * f_in_ptr[(x - 1) * LocalNy + y] +
             1.5 * f_in_ptr[(x + 0) * LocalNy + y]);
      } else if (v_in_ptr[(x + 1) * LocalNy + y] < 0 &&
                 v_in_ptr[(x + 0) * LocalNy + y] < 0) {
        result_ =
            (1.5 * v_in_ptr[(x + 1) * LocalNy + y] -
             .5 * v_in_ptr[(x + 2) * LocalNy + y]) *
            (-1.5 * f_in_ptr[(x + 0) * LocalNy + y] +
             2. * f_in_ptr[(x + 1) * LocalNy + y] - .5 * f_in_ptr[(x + 2) * LocalNy + y]);
      } else {
        result_ = .25 *
                  (v_in_ptr[(x + 1) * LocalNy + y] + v_in_ptr[(x + 0) * LocalNy + y]) *
                  (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]);
      }

      result_ptr[(x + 0) * LocalNy + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_U2_stag_y_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x)*LocalNy + y + 1] > 0 && v_in_ptr[(x)*LocalNy + y + 0] > 0) {
        result_ =
            (1.5 * v_in_ptr[(x)*LocalNy + y + 0] - .5 * v_in_ptr[(x)*LocalNy + y - 1]) *
            (.5 * f_in_ptr[(x)*LocalNy + y - 2] - 2. * f_in_ptr[(x)*LocalNy + y - 1] +
             1.5 * f_in_ptr[(x)*LocalNy + y + 0]);
      } else if (v_in_ptr[(x)*LocalNy + y + 1] < 0 && v_in_ptr[(x)*LocalNy + y + 0] < 0) {
        result_ =
            (1.5 * v_in_ptr[(x)*LocalNy + y + 1] - .5 * v_in_ptr[(x)*LocalNy + y + 2]) *
            (-1.5 * f_in_ptr[(x)*LocalNy + y + 0] + 2. * f_in_ptr[(x)*LocalNy + y + 1] -
             .5 * f_in_ptr[(x)*LocalNy + y + 2]);
      } else {
        result_ = .25 * (v_in_ptr[(x)*LocalNy + y + 1] + v_in_ptr[(x)*LocalNy + y + 0]) *
                  (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]);
      }

      result_ptr[(x)*LocalNy + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_U2_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_x_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            0.5 * (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                   v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]) *
            0.5 * (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_y_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            0.5 * (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                   v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]) *
            0.5 * (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_z_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_x_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          0.5 * (v_in_ptr[(x + 0) * LocalNy + y] + v_in_ptr[(x - 1) * LocalNy + y]) *
          0.5 * (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_y_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          0.5 * (v_in_ptr[(x)*LocalNy + y + 0] + v_in_ptr[(x)*LocalNy + y - 1]) * 0.5 *
          (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_x_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            0.5 * (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
                   v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) *
            0.5 * (f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_y_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            0.5 * (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
                   v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) *
            0.5 * (f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_z_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]);
        }
        for (int z = 1; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              0.5 * (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) *
              0.5 * (f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C2_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_x_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          0.5 * (v_in_ptr[(x + 1) * LocalNy + y] + v_in_ptr[(x + 0) * LocalNy + y]) *
          0.5 * (f_in_ptr[(x + 1) * LocalNy + y] - f_in_ptr[(x - 1) * LocalNy + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C2_stag_y_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          0.5 * (v_in_ptr[(x)*LocalNy + y + 1] + v_in_ptr[(x)*LocalNy + y + 0]) * 0.5 *
          (f_in_ptr[(x)*LocalNy + y + 1] - f_in_ptr[(x)*LocalNy + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C2_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_x_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (9. * (v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                   v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]) -
             v_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
             v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z]) /
            16. * (8. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   8. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                   f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
                   f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_y_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (9. * (v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                   v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]) -
             v_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
             v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z]) /
            16. * (8. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   8. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                   f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
                   f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_z_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 1 + 1 * LocalNz) % LocalNz)] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 2 + 2 * LocalNz) % LocalNz)] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                   +((z - 1 + 1 * LocalNz) % LocalNz)] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 2 + 2 * LocalNz) % LocalNz)] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_x_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (9. * (v_in_ptr[(x - 1) * LocalNy + y] + v_in_ptr[(x + 0) * LocalNy + y]) -
           v_in_ptr[(x - 2) * LocalNy + y] - v_in_ptr[(x + 1) * LocalNy + y]) /
          16. *
          (8. * f_in_ptr[(x + 1) * LocalNy + y] - 8. * f_in_ptr[(x - 1) * LocalNy + y] +
           f_in_ptr[(x - 2) * LocalNy + y] - f_in_ptr[(x + 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_y_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (9. * (v_in_ptr[(x)*LocalNy + y - 1] + v_in_ptr[(x)*LocalNy + y + 0]) -
           v_in_ptr[(x)*LocalNy + y - 2] - v_in_ptr[(x)*LocalNy + y + 1]) /
          16. * (8. * f_in_ptr[(x)*LocalNy + y + 1] - 8. * f_in_ptr[(x)*LocalNy + y - 1] +
                 f_in_ptr[(x)*LocalNy + y - 2] - f_in_ptr[(x)*LocalNy + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_x_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            (9. * (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
                   v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z]) -
             v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
             v_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            16. * (8. * f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
                   8. * f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
                   f_in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
                   f_in_ptr[((x + 2) * LocalNy + y) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_y_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            (9. * (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
                   v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z]) -
             v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
             v_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            16. * (8. * f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
                   8. * f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
                   f_in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
                   f_in_ptr[((x)*LocalNy + y + 2) * LocalNz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_z_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        for (int z = 2; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2]) /
              12.;
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              (9. * (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
                     v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)]) -
               v_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 1 + 1 * LocalNz) % LocalNz)] -
               v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              16. * (8. * f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
                     8. * f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                   +((z - 1 + 1 * LocalNz) % LocalNz)] +
                     f_in_ptr[((x)*LocalNy + y) * LocalNz +
                              +((z - 2 + 2 * LocalNz) % LocalNz)] -
                     f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::VDDX_C4_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_x_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {

      result_ptr[(x + 0) * LocalNy + y] =
          (9. * (v_in_ptr[(x + 0) * LocalNy + y] + v_in_ptr[(x + 1) * LocalNy + y]) -
           v_in_ptr[(x - 1) * LocalNy + y] - v_in_ptr[(x + 2) * LocalNy + y]) /
          16. *
          (8. * f_in_ptr[(x + 1) * LocalNy + y] - 8. * f_in_ptr[(x - 1) * LocalNy + y] +
           f_in_ptr[(x - 2) * LocalNy + y] - f_in_ptr[(x + 2) * LocalNy + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::VDDX_C4_stag_y_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 2; ++y) {

      result_ptr[(x)*LocalNy + y + 0] =
          (9. * (v_in_ptr[(x)*LocalNy + y + 0] + v_in_ptr[(x)*LocalNy + y + 1]) -
           v_in_ptr[(x)*LocalNy + y - 1] - v_in_ptr[(x)*LocalNy + y + 2]) /
          16. * (8. * f_in_ptr[(x)*LocalNy + y + 1] - 8. * f_in_ptr[(x)*LocalNy + y - 1] +
                 f_in_ptr[(x)*LocalNy + y - 2] - f_in_ptr[(x)*LocalNy + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::VDDX_C4_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_x_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]
                               : v_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        result_ -= (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]
                       : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_y_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]
                               : v_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        result_ -= (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]
                       : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_z_Field3D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz +
                        +((z - 1 + 1 * LocalNz) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                 +((z - 1 + 1 * LocalNz) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz +
                             +((z - 1 + 1 * LocalNz) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          result_ -=
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_x_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x - 1) * LocalNy + y] >= 0)
              ? v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y]
              : v_in_ptr[(x - 1) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y];
      result_ -= (v_in_ptr[(x + 0) * LocalNy + y] >= 0)
                     ? v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y]
                     : v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y];

      result_ptr[(x + 0) * LocalNy + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_y_Field2D_CtoL(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x)*LocalNy + y - 1] >= 0)
              ? v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y - 1]
              : v_in_ptr[(x)*LocalNy + y - 1] * f_in_ptr[(x)*LocalNy + y + 0];
      result_ -= (v_in_ptr[(x)*LocalNy + y + 0] >= 0)
                     ? v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 0]
                     : v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 1];

      result_ptr[(x)*LocalNy + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_x_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x - 1) * LocalNy + y) * LocalNz + z]
                               : v_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] *
                                     f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
        result_ -= (v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 0) * LocalNy + y) * LocalNz + z]
                       : v_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] *
                             f_in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];

        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_y_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] >= 0)
                               ? v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y - 1) * LocalNz + z]
                               : v_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] *
                                     f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
        result_ -= (v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] >= 0)
                       ? v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 0) * LocalNz + z]
                       : v_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] *
                             f_in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];

        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_z_Field3D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        for (int z = 1; z < LocalNz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
        {
          int z = LocalNz - 1;
          BoutReal result_ = (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z - 1]
                                 : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] *
                                       f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0];
          result_ -= (v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] >= 0)
                         ? v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 0]
                         : v_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] *
                               f_in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz +
                                 +((z - 1 + 1 * LocalNz) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)];
          result_ -=
              (v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] >= 0)
                  ? v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)]
                  : v_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] *
                        f_in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D AiolosMesh::FDDX_U1_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_x_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x + 0) * LocalNy + y] >= 0)
              ? v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x - 1) * LocalNy + y]
              : v_in_ptr[(x + 0) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y];
      result_ -= (v_in_ptr[(x + 1) * LocalNy + y] >= 0)
                     ? v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 0) * LocalNy + y]
                     : v_in_ptr[(x + 1) * LocalNy + y] * f_in_ptr[(x + 1) * LocalNy + y];

      result_ptr[(x + 0) * LocalNy + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

void AiolosMesh::FDDX_U1_stag_y_Field2D_LtoC(
    BoutReal *__restrict__ result_ptr, const BoutReal *__restrict__ v_in_ptr,
    const BoutReal *__restrict__ f_in_ptr) const {
#if CHECK > 0
  if (ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 1; ++y) {
      BoutReal result_ =
          (v_in_ptr[(x)*LocalNy + y + 0] >= 0)
              ? v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y - 1]
              : v_in_ptr[(x)*LocalNy + y + 0] * f_in_ptr[(x)*LocalNy + y + 0];
      result_ -= (v_in_ptr[(x)*LocalNy + y + 1] >= 0)
                     ? v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 0]
                     : v_in_ptr[(x)*LocalNy + y + 1] * f_in_ptr[(x)*LocalNy + y + 1];

      result_ptr[(x)*LocalNy + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D AiolosMesh::FDDX_U1_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) const {
  ASSERT1(this == v_in.getMesh());
  ASSERT1(this == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result((AiolosMesh *)this);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}
