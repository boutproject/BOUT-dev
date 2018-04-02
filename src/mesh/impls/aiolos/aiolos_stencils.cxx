
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

static void DDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            0.5 *
            (in_ptr[((x + 1) * Ny + y) * Nz + z] - in_ptr[((x - 1) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            0.5 * (in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 *
              (in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                     in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          0.5 * (in_ptr[(x + 1) * Ny + y] - in_ptr[(x - 1) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          0.5 * (in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            (8. * in_ptr[((x + 1) * Ny + y) * Nz + z] -
             8. * in_ptr[((x - 1) * Ny + y) * Nz + z] +
             in_ptr[((x - 2) * Ny + y) * Nz + z] - in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            (8. * in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * in_ptr[((x)*Ny + y - 1) * Nz + z] + in_ptr[((x)*Ny + y - 2) * Nz + z] -
             in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2] - in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2] -
               in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2] -
               in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               8. * in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
               in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (8. * in_ptr[(x + 1) * Ny + y] - 8. * in_ptr[(x - 1) * Ny + y] +
           in_ptr[(x - 2) * Ny + y] - in_ptr[(x + 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (8. * in_ptr[(x)*Ny + y + 1] - 8. * in_ptr[(x)*Ny + y - 1] +
           in_ptr[(x)*Ny + y - 2] - in_ptr[(x)*Ny + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_CWENO2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal isl, isr, isc;
        BoutReal al, ar, ac, sa;
        BoutReal dl, dr, dc;
        dc = 0.5 *
             (in_ptr[((x + 1) * Ny + y) * Nz + z] - in_ptr[((x - 1) * Ny + y) * Nz + z]);
        dl = in_ptr[((x + 0) * Ny + y) * Nz + z] - in_ptr[((x - 1) * Ny + y) * Nz + z];
        dr = in_ptr[((x + 1) * Ny + y) * Nz + z] - in_ptr[((x + 0) * Ny + y) * Nz + z];
        isl = SQ(dl);
        isr = SQ(dr);
        isc = (13. / 3.) * SQ(in_ptr[((x + 1) * Ny + y) * Nz + z] -
                              2. * in_ptr[((x + 0) * Ny + y) * Nz + z] +
                              in_ptr[((x - 1) * Ny + y) * Nz + z]) +
              0.25 * SQ(in_ptr[((x + 1) * Ny + y) * Nz + z] -
                        in_ptr[((x - 1) * Ny + y) * Nz + z]);
        al = 0.25 / SQ(WENO_SMALL + isl);
        ar = 0.25 / SQ(WENO_SMALL + isr);
        ac = 0.5 / SQ(WENO_SMALL + isc);
        sa = al + ar + ac;

        result_ptr[((x + 0) * Ny + y) * Nz + z] = (al * dl + ar * dr + ac * dc) / sa;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_CWENO2_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_CWENO2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_CWENO2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal isl, isr, isc;
        BoutReal al, ar, ac, sa;
        BoutReal dl, dr, dc;
        dc =
            0.5 * (in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z]);
        dl = in_ptr[((x)*Ny + y + 0) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z];
        dr = in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y + 0) * Nz + z];
        isl = SQ(dl);
        isr = SQ(dr);
        isc = (13. / 3.) * SQ(in_ptr[((x)*Ny + y + 1) * Nz + z] -
                              2. * in_ptr[((x)*Ny + y + 0) * Nz + z] +
                              in_ptr[((x)*Ny + y - 1) * Nz + z]) +
              0.25 * SQ(in_ptr[((x)*Ny + y + 1) * Nz + z] -
                        in_ptr[((x)*Ny + y - 1) * Nz + z]);
        al = 0.25 / SQ(WENO_SMALL + isl);
        ar = 0.25 / SQ(WENO_SMALL + isr);
        ac = 0.5 / SQ(WENO_SMALL + isc);
        sa = al + ar + ac;

        result_ptr[((x)*Ny + y + 0) * Nz + z] = (al * dl + ar * dr + ac * dc) / sa;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_CWENO2_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_CWENO2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_CWENO2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
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
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
          dl = in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz];
          dr = in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                2. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) +
                0.25 * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1] -
                          in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = (al * dl + ar * dr + ac * dc) / sa;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 *
               (in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z - 1]);
          dl = in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1];
          dr = in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                2. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                in_ptr[((x)*Ny + y) * Nz + z - 1]) +
                0.25 * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1] -
                          in_ptr[((x)*Ny + y) * Nz + z - 1]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = (al * dl + ar * dr + ac * dc) / sa;
        }
        {
          int z = Nz - 1;
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1]);
          dl = in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1];
          dr = in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] - in_ptr[((x)*Ny + y) * Nz + z + 0];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                                2. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                in_ptr[((x)*Ny + y) * Nz + z - 1]) +
                0.25 * SQ(in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                          in_ptr[((x)*Ny + y) * Nz + z - 1]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = (al * dl + ar * dr + ac * dc) / sa;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal isl, isr, isc;
          BoutReal al, ar, ac, sa;
          BoutReal dl, dr, dc;
          dc = 0.5 * (in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                      in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
          dl = in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)];
          dr = in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          isl = SQ(dl);
          isr = SQ(dr);
          isc = (13. / 3.) * SQ(in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                                2. * in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                                in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) +
                0.25 * SQ(in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                          in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
          al = 0.25 / SQ(WENO_SMALL + isl);
          ar = 0.25 / SQ(WENO_SMALL + isr);
          ac = 0.5 / SQ(WENO_SMALL + isc);
          sa = al + ar + ac;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = (al * dl + ar * dr + ac * dc) / sa;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_CWENO2_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_CWENO2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::DDX_CWENO2_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_CWENO2_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_CWENO2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal isl, isr, isc;
      BoutReal al, ar, ac, sa;
      BoutReal dl, dr, dc;
      dc = 0.5 * (in_ptr[(x + 1) * Ny + y] - in_ptr[(x - 1) * Ny + y]);
      dl = in_ptr[(x + 0) * Ny + y] - in_ptr[(x - 1) * Ny + y];
      dr = in_ptr[(x + 1) * Ny + y] - in_ptr[(x + 0) * Ny + y];
      isl = SQ(dl);
      isr = SQ(dr);
      isc = (13. / 3.) * SQ(in_ptr[(x + 1) * Ny + y] - 2. * in_ptr[(x + 0) * Ny + y] +
                            in_ptr[(x - 1) * Ny + y]) +
            0.25 * SQ(in_ptr[(x + 1) * Ny + y] - in_ptr[(x - 1) * Ny + y]);
      al = 0.25 / SQ(WENO_SMALL + isl);
      ar = 0.25 / SQ(WENO_SMALL + isr);
      ac = 0.5 / SQ(WENO_SMALL + isc);
      sa = al + ar + ac;

      result_ptr[(x + 0) * Ny + y] = (al * dl + ar * dr + ac * dc) / sa;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_CWENO2_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_CWENO2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_CWENO2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_CWENO2_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_CWENO2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal isl, isr, isc;
      BoutReal al, ar, ac, sa;
      BoutReal dl, dr, dc;
      dc = 0.5 * (in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y - 1]);
      dl = in_ptr[(x)*Ny + y + 0] - in_ptr[(x)*Ny + y - 1];
      dr = in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y + 0];
      isl = SQ(dl);
      isr = SQ(dr);
      isc = (13. / 3.) * SQ(in_ptr[(x)*Ny + y + 1] - 2. * in_ptr[(x)*Ny + y + 0] +
                            in_ptr[(x)*Ny + y - 1]) +
            0.25 * SQ(in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y - 1]);
      al = 0.25 / SQ(WENO_SMALL + isl);
      ar = 0.25 / SQ(WENO_SMALL + isr);
      ac = 0.5 / SQ(WENO_SMALL + isc);
      sa = al + ar + ac;

      result_ptr[(x)*Ny + y + 0] = (al * dl + ar * dr + ac * dc) / sa;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_CWENO2_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_CWENO2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::DDX_CWENO2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_CWENO2_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_S2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ =
            (8. * in_ptr[((x + 1) * Ny + y) * Nz + z] -
             8. * in_ptr[((x - 1) * Ny + y) * Nz + z] +
             in_ptr[((x - 2) * Ny + y) * Nz + z] - in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
        result_ += SIGN(in_ptr[((x + 0) * Ny + y) * Nz + z]) *
                   (in_ptr[((x + 2) * Ny + y) * Nz + z] -
                    4. * in_ptr[((x + 1) * Ny + y) * Nz + z] +
                    6. * in_ptr[((x + 0) * Ny + y) * Nz + z] -
                    4. * in_ptr[((x - 1) * Ny + y) * Nz + z] +
                    in_ptr[((x - 2) * Ny + y) * Nz + z]) /
                   12.;

        result_ptr[((x + 0) * Ny + y) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_S2_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_S2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_S2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ =
            (8. * in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * in_ptr[((x)*Ny + y - 1) * Nz + z] + in_ptr[((x)*Ny + y - 2) * Nz + z] -
             in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
        result_ +=
            SIGN(in_ptr[((x)*Ny + y + 0) * Nz + z]) *
            (in_ptr[((x)*Ny + y + 2) * Nz + z] - 4. * in_ptr[((x)*Ny + y + 1) * Nz + z] +
             6. * in_ptr[((x)*Ny + y + 0) * Nz + z] -
             4. * in_ptr[((x)*Ny + y - 1) * Nz + z] + in_ptr[((x)*Ny + y - 2) * Nz + z]) /
            12.;

        result_ptr[((x)*Ny + y + 0) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_S2_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_S2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_S2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_ = (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
                              8. * in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                              in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                              in_ptr[((x)*Ny + y) * Nz + z + 2]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                     (in_ptr[((x)*Ny + y) * Nz + z + 2] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z + 1] +
                      6. * in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                      in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_ = (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
                              8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                              in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                              in_ptr[((x)*Ny + y) * Nz + z + 2]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                     (in_ptr[((x)*Ny + y) * Nz + z + 2] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z + 1] +
                      6. * in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                      in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        for (int z = 2; z < Nz - 2; ++z) {
          BoutReal result_ =
              (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2] - in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                     (in_ptr[((x)*Ny + y) * Nz + z + 2] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z + 1] +
                      6. * in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                      in_ptr[((x)*Ny + y) * Nz + z - 2]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 2;
          BoutReal result_ = (8. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
                              8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                              in_ptr[((x)*Ny + y) * Nz + z - 2] -
                              in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                     (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z + 1] +
                      6. * in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                      in_ptr[((x)*Ny + y) * Nz + z - 2]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_ = (8. * in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                              8. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                              in_ptr[((x)*Ny + y) * Nz + z - 2] -
                              in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                     (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                      6. * in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      4. * in_ptr[((x)*Ny + y) * Nz + z - 1] +
                      in_ptr[((x)*Ny + y) * Nz + z - 2]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_ = (8. * in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                              8. * in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                              in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
                              in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
                             12.;
          result_ += SIGN(in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) *
                     (in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] -
                      4. * in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
                      6. * in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                      4. * in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                      in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)]) /
                     12.;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_S2_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_S2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::DDX_S2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_S2_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_S2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_ = (8. * in_ptr[(x + 1) * Ny + y] - 8. * in_ptr[(x - 1) * Ny + y] +
                          in_ptr[(x - 2) * Ny + y] - in_ptr[(x + 2) * Ny + y]) /
                         12.;
      result_ += SIGN(in_ptr[(x + 0) * Ny + y]) *
                 (in_ptr[(x + 2) * Ny + y] - 4. * in_ptr[(x + 1) * Ny + y] +
                  6. * in_ptr[(x + 0) * Ny + y] - 4. * in_ptr[(x - 1) * Ny + y] +
                  in_ptr[(x - 2) * Ny + y]) /
                 12.;

      result_ptr[(x + 0) * Ny + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_S2_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_S2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_S2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_S2_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void DDX_S2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                  const BoutReal *__restrict__ in_ptr, Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {
      BoutReal result_ = (8. * in_ptr[(x)*Ny + y + 1] - 8. * in_ptr[(x)*Ny + y - 1] +
                          in_ptr[(x)*Ny + y - 2] - in_ptr[(x)*Ny + y + 2]) /
                         12.;
      result_ += SIGN(in_ptr[(x)*Ny + y + 0]) *
                 (in_ptr[(x)*Ny + y + 2] - 4. * in_ptr[(x)*Ny + y + 1] +
                  6. * in_ptr[(x)*Ny + y + 0] - 4. * in_ptr[(x)*Ny + y - 1] +
                  in_ptr[(x)*Ny + y - 2]) /
                 12.;

      result_ptr[(x)*Ny + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_S2_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_S2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::DDX_S2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_S2_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            in_ptr[((x + 1) * Ny + y) * Nz + z] + in_ptr[((x - 1) * Ny + y) * Nz + z] -
            2. * in_ptr[((x + 0) * Ny + y) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] = in_ptr[((x)*Ny + y + 1) * Nz + z] +
                                                in_ptr[((x)*Ny + y - 1) * Nz + z] -
                                                2. * in_ptr[((x)*Ny + y + 0) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
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

          result_ptr[((x)*Ny + y) * Nz + z + 0] = in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                                  in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
                                                  2. * in_ptr[((x)*Ny + y) * Nz + z + 0];
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] = in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                                  in_ptr[((x)*Ny + y) * Nz + z - 1] -
                                                  2. * in_ptr[((x)*Ny + y) * Nz + z + 0];
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] = in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                                                  in_ptr[((x)*Ny + y) * Nz + z - 1] -
                                                  2. * in_ptr[((x)*Ny + y) * Nz + z + 0];
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
              in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
              2. * in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] = in_ptr[(x + 1) * Ny + y] + in_ptr[(x - 1) * Ny + y] -
                                     2. * in_ptr[(x + 0) * Ny + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          in_ptr[(x)*Ny + y + 1] + in_ptr[(x)*Ny + y - 1] - 2. * in_ptr[(x)*Ny + y + 0];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            (-in_ptr[((x + 2) * Ny + y) * Nz + z] +
             16. * in_ptr[((x + 1) * Ny + y) * Nz + z] -
             30. * in_ptr[((x + 0) * Ny + y) * Nz + z] +
             16. * in_ptr[((x - 1) * Ny + y) * Nz + z] -
             in_ptr[((x - 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C4_x_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_x_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] = (-in_ptr[((x)*Ny + y + 2) * Nz + z] +
                                                 16. * in_ptr[((x)*Ny + y + 1) * Nz + z] -
                                                 30. * in_ptr[((x)*Ny + y + 0) * Nz + z] +
                                                 16. * in_ptr[((x)*Ny + y - 1) * Nz + z] -
                                                 in_ptr[((x)*Ny + y - 2) * Nz + z]) /
                                                12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C4_y_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_y_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
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
              (-in_ptr[((x)*Ny + y) * Nz + z + 2] +
               16. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               30. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
               16. * in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (-in_ptr[((x)*Ny + y) * Nz + z + 2] +
               16. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               30. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
               16. * in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (-in_ptr[((x)*Ny + y) * Nz + z + 2] +
               16. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               30. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
               16. * in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z - 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (-in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
               16. * in_ptr[((x)*Ny + y) * Nz + z + 1] -
               30. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
               16. * in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z - 2]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (-in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
               16. * in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               30. * in_ptr[((x)*Ny + y) * Nz + z + 0] +
               16. * in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z - 2]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (-in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] +
               16. * in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               30. * in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
               16. * in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C4_z_norm(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C4_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C4_z_Field3D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (-in_ptr[(x + 2) * Ny + y] + 16. * in_ptr[(x + 1) * Ny + y] -
           30. * in_ptr[(x + 0) * Ny + y] + 16. * in_ptr[(x - 1) * Ny + y] -
           in_ptr[(x - 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C4_x_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C4_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C4_x_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void D2DX2_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                    const BoutReal *__restrict__ in_ptr,
                                    Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (-in_ptr[(x)*Ny + y + 2] + 16. * in_ptr[(x)*Ny + y + 1] -
           30. * in_ptr[(x)*Ny + y + 0] + 16. * in_ptr[(x)*Ny + y - 1] -
           in_ptr[(x)*Ny + y - 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C4_y_norm(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C4_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C4_y_Field2D_norm(result_ptr, in_ptr, localmesh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] = v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                                                  0.5 *
                                                  (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                                                   f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] * 0.5 *
            (f_in_ptr[((x)*Ny + y + 1) * Nz + z] - f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * 0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * 0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + z + 1] - f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * 0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] * 0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          v_in_ptr[(x + 0) * Ny + y] * 0.5 *
          (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] = v_in_ptr[(x)*Ny + y + 0] * 0.5 *
                                   (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
            (8. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
             8. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
             f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
             f_in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
            (8. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
             f_in_ptr[((x)*Ny + y - 2) * Nz + z] - f_in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
               f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
              (8. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               8. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
               f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
               f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          v_in_ptr[(x + 0) * Ny + y] *
          (8. * f_in_ptr[(x + 1) * Ny + y] - 8. * f_in_ptr[(x - 1) * Ny + y] +
           f_in_ptr[(x - 2) * Ny + y] - f_in_ptr[(x + 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          v_in_ptr[(x)*Ny + y + 0] *
          (8. * f_in_ptr[(x)*Ny + y + 1] - 8. * f_in_ptr[(x)*Ny + y - 1] +
           f_in_ptr[(x)*Ny + y - 2] - f_in_ptr[(x)*Ny + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U1_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0.0
                ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (f_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                       f_in_ptr[((x - 1) * Ny + y) * Nz + z])
                : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                       f_in_ptr[((x + 0) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U1_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0.0
                ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (f_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                       f_in_ptr[((x)*Ny + y - 1) * Nz + z])
                : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                       f_in_ptr[((x)*Ny + y + 0) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U1_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         f_in_ptr[((x)*Ny + y) * Nz + z - 1])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         f_in_ptr[((x)*Ny + y) * Nz + z - 1])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                         f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                         f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)])
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                         f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U1_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          v_in_ptr[(x + 0) * Ny + y] >= 0.0
              ? v_in_ptr[(x + 0) * Ny + y] *
                    (f_in_ptr[(x + 0) * Ny + y] - f_in_ptr[(x - 1) * Ny + y])
              : v_in_ptr[(x + 0) * Ny + y] *
                    (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x + 0) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U1_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          v_in_ptr[(x)*Ny + y + 0] >= 0.0
              ? v_in_ptr[(x)*Ny + y + 0] *
                    (f_in_ptr[(x)*Ny + y + 0] - f_in_ptr[(x)*Ny + y - 1])
              : v_in_ptr[(x)*Ny + y + 0] *
                    (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y + 0]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0.0
                ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                       2.0 * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                       0.5 * f_in_ptr[((x - 2) * Ny + y) * Nz + z])
                : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (-0.5 * f_in_ptr[((x + 2) * Ny + y) * Nz + z] +
                       2.0 * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                       1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0.0
                ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                       2.0 * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                       0.5 * f_in_ptr[((x)*Ny + y - 2) * Nz + z])
                : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (-0.5 * f_in_ptr[((x)*Ny + y + 2) * Nz + z] +
                       2.0 * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                       1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2])
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                         0.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)])
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (-0.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] +
                         2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                         1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          v_in_ptr[(x + 0) * Ny + y] >= 0.0
              ? v_in_ptr[(x + 0) * Ny + y] *
                    (1.5 * f_in_ptr[(x + 0) * Ny + y] - 2.0 * f_in_ptr[(x - 1) * Ny + y] +
                     0.5 * f_in_ptr[(x - 2) * Ny + y])
              : v_in_ptr[(x + 0) * Ny + y] *
                    (-0.5 * f_in_ptr[(x + 2) * Ny + y] +
                     2.0 * f_in_ptr[(x + 1) * Ny + y] - 1.5 * f_in_ptr[(x + 0) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          v_in_ptr[(x)*Ny + y + 0] >= 0.0
              ? v_in_ptr[(x)*Ny + y + 0] *
                    (1.5 * f_in_ptr[(x)*Ny + y + 0] - 2.0 * f_in_ptr[(x)*Ny + y - 1] +
                     0.5 * f_in_ptr[(x)*Ny + y - 2])
              : v_in_ptr[(x)*Ny + y + 0] *
                    (-0.5 * f_in_ptr[(x)*Ny + y + 2] + 2.0 * f_in_ptr[(x)*Ny + y + 1] -
                     1.5 * f_in_ptr[(x)*Ny + y + 0]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U3_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0.0
                ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (4. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                       12. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                       2. * f_in_ptr[((x - 2) * Ny + y) * Nz + z] +
                       6. * f_in_ptr[((x + 0) * Ny + y) * Nz + z]) /
                      12.
                : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                      (-4. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                       12. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                       2. * f_in_ptr[((x + 2) * Ny + y) * Nz + z] -
                       6. * f_in_ptr[((x + 0) * Ny + y) * Nz + z]) /
                      12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U3_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U3_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0.0
                ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (4. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                       12. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                       2. * f_in_ptr[((x)*Ny + y - 2) * Nz + z] +
                       6. * f_in_ptr[((x)*Ny + y + 0) * Nz + z]) /
                      12.
                : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                      (-4. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                       12. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                       2. * f_in_ptr[((x)*Ny + y + 2) * Nz + z] -
                       6. * f_in_ptr[((x)*Ny + y + 0) * Nz + z]) /
                      12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U3_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U3_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + z + 0]) /
                        12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0.0
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (4. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                         12. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                         2. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] +
                         6. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) /
                        12.
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        (-4. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                         12. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                         2. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] -
                         6. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) /
                        12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U3_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U3_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U3_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U3_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          v_in_ptr[(x + 0) * Ny + y] >= 0.0
              ? v_in_ptr[(x + 0) * Ny + y] *
                    (4. * f_in_ptr[(x + 1) * Ny + y] - 12. * f_in_ptr[(x - 1) * Ny + y] +
                     2. * f_in_ptr[(x - 2) * Ny + y] + 6. * f_in_ptr[(x + 0) * Ny + y]) /
                    12.
              : v_in_ptr[(x + 0) * Ny + y] *
                    (-4. * f_in_ptr[(x - 1) * Ny + y] + 12. * f_in_ptr[(x + 1) * Ny + y] -
                     2. * f_in_ptr[(x + 2) * Ny + y] - 6. * f_in_ptr[(x + 0) * Ny + y]) /
                    12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U3_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U3_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U3_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_U3_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          v_in_ptr[(x)*Ny + y + 0] >= 0.0
              ? v_in_ptr[(x)*Ny + y + 0] *
                    (4. * f_in_ptr[(x)*Ny + y + 1] - 12. * f_in_ptr[(x)*Ny + y - 1] +
                     2. * f_in_ptr[(x)*Ny + y - 2] + 6. * f_in_ptr[(x)*Ny + y + 0]) /
                    12.
              : v_in_ptr[(x)*Ny + y + 0] *
                    (-4. * f_in_ptr[(x)*Ny + y - 1] + 12. * f_in_ptr[(x)*Ny + y + 1] -
                     2. * f_in_ptr[(x)*Ny + y + 2] - 6. * f_in_ptr[(x)*Ny + y + 0]) /
                    12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U3_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U3_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U3_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U3_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_WENO3_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ v_in_ptr,
                                      const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal deriv, w, r;
        if (v_in_ptr[((x + 0) * Ny + y) * Nz + z] > 0.0) {
          r = (WENO_SMALL + SQ(f_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                               2.0 * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                               f_in_ptr[((x - 2) * Ny + y) * Nz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                               2.0 * f_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                               f_in_ptr[((x - 1) * Ny + y) * Nz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                         f_in_ptr[((x - 1) * Ny + y) * Nz + z]) -
                  0.5 * w * (-f_in_ptr[((x - 2) * Ny + y) * Nz + z] +
                             3. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] -
                             3. * f_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                             f_in_ptr[((x + 1) * Ny + y) * Nz + z]);
        } else {
          r = (WENO_SMALL + SQ(f_in_ptr[((x + 2) * Ny + y) * Nz + z] -
                               2.0 * f_in_ptr[((x + 1) * Ny + y) * Nz + z] +
                               f_in_ptr[((x + 0) * Ny + y) * Nz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                               2.0 * f_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                               f_in_ptr[((x - 1) * Ny + y) * Nz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                         f_in_ptr[((x - 1) * Ny + y) * Nz + z]) -
                  0.5 * w * (-f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                             3. * f_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                             3. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] +
                             f_in_ptr[((x + 2) * Ny + y) * Nz + z]);
        }

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] * deriv;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_WENO3_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_WENO3_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ v_in_ptr,
                                      const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal deriv, w, r;
        if (v_in_ptr[((x)*Ny + y + 0) * Nz + z] > 0.0) {
          r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                               2.0 * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                               f_in_ptr[((x)*Ny + y - 2) * Nz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                               2.0 * f_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                               f_in_ptr[((x)*Ny + y - 1) * Nz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                         f_in_ptr[((x)*Ny + y - 1) * Nz + z]) -
                  0.5 * w * (-f_in_ptr[((x)*Ny + y - 2) * Nz + z] +
                             3. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] -
                             3. * f_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                             f_in_ptr[((x)*Ny + y + 1) * Nz + z]);
        } else {
          r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y + 2) * Nz + z] -
                               2.0 * f_in_ptr[((x)*Ny + y + 1) * Nz + z] +
                               f_in_ptr[((x)*Ny + y + 0) * Nz + z])) /
              (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                               2.0 * f_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                               f_in_ptr[((x)*Ny + y - 1) * Nz + z]));
          w = 1.0 / (1.0 + 2.0 * r * r);
          deriv = 0.5 * (f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                         f_in_ptr[((x)*Ny + y - 1) * Nz + z]) -
                  0.5 * w * (-f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                             3. * f_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                             3. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] +
                             f_in_ptr[((x)*Ny + y + 2) * Nz + z]);
        }

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] * deriv;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_WENO3_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_WENO3_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ v_in_ptr,
                                      const BoutReal *__restrict__ f_in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * deriv;
        }
        {
          int z = 1;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * deriv;
        }
        for (int z = 2; z < Nz - 2; ++z) {
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 2] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * deriv;
        }
        {
          int z = Nz - 2;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * deriv;
        }
        {
          int z = Nz - 1;
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0.0) {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 2])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 2] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z + 0])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                                 f_in_ptr[((x)*Ny + y) * Nz + z - 1]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                           f_in_ptr[((x)*Ny + y) * Nz + z - 1]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                               f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] * deriv;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal deriv, w, r;
          if (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] > 0.0) {
            r = (WENO_SMALL +
                 SQ(f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                    2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                    f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                                 f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv =
                0.5 * (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) -
                0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] +
                           3. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
                           3. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                           f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)]);
          } else {
            r = (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
                                 f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)])) /
                (WENO_SMALL + SQ(f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                                 2.0 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                                 f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]));
            w = 1.0 / (1.0 + 2.0 * r * r);
            deriv = 0.5 * (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                           f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) -
                    0.5 * w * (-f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                               3. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                               3. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] * deriv;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_WENO3_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_WENO3_z_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_WENO3_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_WENO3_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ v_in_ptr,
                                      const BoutReal *__restrict__ f_in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal deriv, w, r;
      if (v_in_ptr[(x + 0) * Ny + y] > 0.0) {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x + 0) * Ny + y] - 2.0 * f_in_ptr[(x - 1) * Ny + y] +
                f_in_ptr[(x - 2) * Ny + y])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x + 1) * Ny + y] - 2.0 * f_in_ptr[(x + 0) * Ny + y] +
                f_in_ptr[(x - 1) * Ny + y]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]) -
                0.5 * w * (-f_in_ptr[(x - 2) * Ny + y] + 3. * f_in_ptr[(x - 1) * Ny + y] -
                           3. * f_in_ptr[(x + 0) * Ny + y] + f_in_ptr[(x + 1) * Ny + y]);
      } else {
        r = (WENO_SMALL +
             SQ(f_in_ptr[(x + 2) * Ny + y] - 2.0 * f_in_ptr[(x + 1) * Ny + y] +
                f_in_ptr[(x + 0) * Ny + y])) /
            (WENO_SMALL +
             SQ(f_in_ptr[(x + 1) * Ny + y] - 2.0 * f_in_ptr[(x + 0) * Ny + y] +
                f_in_ptr[(x - 1) * Ny + y]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]) -
                0.5 * w * (-f_in_ptr[(x - 1) * Ny + y] + 3. * f_in_ptr[(x + 0) * Ny + y] -
                           3. * f_in_ptr[(x + 1) * Ny + y] + f_in_ptr[(x + 2) * Ny + y]);
      }

      result_ptr[(x + 0) * Ny + y] = v_in_ptr[(x + 0) * Ny + y] * deriv;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_WENO3_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_WENO3_x_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_WENO3_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void VDDX_WENO3_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                      const BoutReal *__restrict__ v_in_ptr,
                                      const BoutReal *__restrict__ f_in_ptr,
                                      Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {
      BoutReal deriv, w, r;
      if (v_in_ptr[(x)*Ny + y + 0] > 0.0) {
        r = (WENO_SMALL + SQ(f_in_ptr[(x)*Ny + y + 0] - 2.0 * f_in_ptr[(x)*Ny + y - 1] +
                             f_in_ptr[(x)*Ny + y - 2])) /
            (WENO_SMALL + SQ(f_in_ptr[(x)*Ny + y + 1] - 2.0 * f_in_ptr[(x)*Ny + y + 0] +
                             f_in_ptr[(x)*Ny + y - 1]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]) -
                0.5 * w * (-f_in_ptr[(x)*Ny + y - 2] + 3. * f_in_ptr[(x)*Ny + y - 1] -
                           3. * f_in_ptr[(x)*Ny + y + 0] + f_in_ptr[(x)*Ny + y + 1]);
      } else {
        r = (WENO_SMALL + SQ(f_in_ptr[(x)*Ny + y + 2] - 2.0 * f_in_ptr[(x)*Ny + y + 1] +
                             f_in_ptr[(x)*Ny + y + 0])) /
            (WENO_SMALL + SQ(f_in_ptr[(x)*Ny + y + 1] - 2.0 * f_in_ptr[(x)*Ny + y + 0] +
                             f_in_ptr[(x)*Ny + y - 1]));
        w = 1.0 / (1.0 + 2.0 * r * r);
        deriv = 0.5 * (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]) -
                0.5 * w * (-f_in_ptr[(x)*Ny + y - 1] + 3. * f_in_ptr[(x)*Ny + y + 0] -
                           3. * f_in_ptr[(x)*Ny + y + 1] + f_in_ptr[(x)*Ny + y + 2]);
      }

      result_ptr[(x)*Ny + y + 0] = v_in_ptr[(x)*Ny + y + 0] * deriv;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_WENO3_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_WENO3_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_WENO3_y_norm - Not enough guards cells "
                        "to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_WENO3_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_U1_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal vs = 0.5 * (v_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                             v_in_ptr[((x + 0) * Ny + y) * Nz + z]);
        BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[((x - 1) * Ny + y) * Nz + z]
                                       : vs * f_in_ptr[((x + 0) * Ny + y) * Nz + z];
        vs = 0.5 * (v_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                    v_in_ptr[((x + 1) * Ny + y) * Nz + z]);
        result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x + 0) * Ny + y) * Nz + z]
                               : vs * f_in_ptr[((x + 1) * Ny + y) * Nz + z];

        result_ptr[((x + 0) * Ny + y) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_U1_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal vs = 0.5 * (v_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                             v_in_ptr[((x)*Ny + y + 0) * Nz + z]);
        BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y - 1) * Nz + z]
                                       : vs * f_in_ptr[((x)*Ny + y + 0) * Nz + z];
        vs = 0.5 *
             (v_in_ptr[((x)*Ny + y + 0) * Nz + z] + v_in_ptr[((x)*Ny + y + 1) * Nz + z]);
        result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y + 0) * Nz + z]
                               : vs * f_in_ptr[((x)*Ny + y + 1) * Nz + z];

        result_ptr[((x)*Ny + y + 0) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_U1_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                               v_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]
                                         : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                      v_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                                 : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               v_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                         : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                      v_in_ptr[((x)*Ny + y) * Nz + z + 1]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                                 : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        {
          int z = Nz - 1;
          BoutReal vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                               v_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          BoutReal result_ = (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                         : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                      v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                                 : vs * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                               v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);
          BoutReal result_ =
              (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]
                          : vs * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          vs = 0.5 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                      v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)]);
          result_ -= (vs >= 0.0) ? vs * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]
                                 : vs * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_U1_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal vs = 0.5 * (v_in_ptr[(x - 1) * Ny + y] + v_in_ptr[(x + 0) * Ny + y]);
      BoutReal result_ =
          (vs >= 0.0) ? vs * f_in_ptr[(x - 1) * Ny + y] : vs * f_in_ptr[(x + 0) * Ny + y];
      vs = 0.5 * (v_in_ptr[(x + 0) * Ny + y] + v_in_ptr[(x + 1) * Ny + y]);
      result_ -=
          (vs >= 0.0) ? vs * f_in_ptr[(x + 0) * Ny + y] : vs * f_in_ptr[(x + 1) * Ny + y];

      result_ptr[(x + 0) * Ny + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_U1_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal vs = 0.5 * (v_in_ptr[(x)*Ny + y - 1] + v_in_ptr[(x)*Ny + y + 0]);
      BoutReal result_ =
          (vs >= 0.0) ? vs * f_in_ptr[(x)*Ny + y - 1] : vs * f_in_ptr[(x)*Ny + y + 0];
      vs = 0.5 * (v_in_ptr[(x)*Ny + y + 0] + v_in_ptr[(x)*Ny + y + 1]);
      result_ -=
          (vs >= 0.0) ? vs * f_in_ptr[(x)*Ny + y + 0] : vs * f_in_ptr[(x)*Ny + y + 1];

      result_ptr[(x)*Ny + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C2_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            0.5 * (v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                       f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                   v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                       f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C2_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C2_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            0.5 *
            (v_in_ptr[((x)*Ny + y + 1) * Nz + z] * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
             v_in_ptr[((x)*Ny + y - 1) * Nz + z] * f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C2_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C2_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              0.5 *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 1] * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 1] * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] * f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] * f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                         f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                         f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C2_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_C2_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C2_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C2_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          0.5 * (v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 1) * Ny + y] -
                 v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x - 1) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_C2_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_C2_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C2_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C2_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          0.5 * (v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 1] -
                 v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_C2_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C2_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_C2_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C2_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C4_x_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            (8. * v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                 f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
             8. * v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                 f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
             v_in_ptr[((x - 2) * Ny + y) * Nz + z] *
                 f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
             v_in_ptr[((x + 2) * Ny + y) * Nz + z] *
                 f_in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C4_x_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_x_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_x_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C4_y_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            (8. * v_in_ptr[((x)*Ny + y + 1) * Nz + z] *
                 f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * v_in_ptr[((x)*Ny + y - 1) * Nz + z] *
                 f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
             v_in_ptr[((x)*Ny + y - 2) * Nz + z] * f_in_ptr[((x)*Ny + y - 2) * Nz + z] -
             v_in_ptr[((x)*Ny + y + 2) * Nz + z] * f_in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C4_y_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_y_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_y_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C4_z_Field3D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
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
              (8. * v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
               v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                   f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] *
                   f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (8. * v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                   f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               8. * v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                   f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
               v_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] *
                   f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
               v_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] *
                   f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_C4_z_norm(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_z_Field3D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::FDDX_C4_z_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_C4_z_Field3D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C4_x_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (8. * v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 1) * Ny + y] -
           8. * v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x - 1) * Ny + y] +
           v_in_ptr[(x - 2) * Ny + y] * f_in_ptr[(x - 2) * Ny + y] -
           v_in_ptr[(x + 2) * Ny + y] * f_in_ptr[(x + 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_C4_x_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_x_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::FDDX_C4_x_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C4_x_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void FDDX_C4_y_Field2D_norm(BoutReal *__restrict__ result_ptr,
                                   const BoutReal *__restrict__ v_in_ptr,
                                   const BoutReal *__restrict__ f_in_ptr,
                                   Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (8. * v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 1] -
           8. * v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y - 1] +
           v_in_ptr[(x)*Ny + y - 2] * f_in_ptr[(x)*Ny + y - 2] -
           v_in_ptr[(x)*Ny + y + 2] * f_in_ptr[(x)*Ny + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_C4_y_norm(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_C4_y_Field2D_norm!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::FDDX_C4_y_norm - Not enough guards cells to "
                        "take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_C4_y_Field2D_norm(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

static void DDX_C2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx + 0; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            in_ptr[((x + 0) * Ny + y) * Nz + z] - in_ptr[((x - 1) * Ny + y) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_x_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_x_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void DDX_C2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny + 0; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            in_ptr[((x)*Ny + y + 0) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_y_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_y_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void DDX_C2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
              in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz];
        }
        for (int z = 1; z < Nz - 0; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1];
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
              in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_z_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_z_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void DDX_C2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx + 0; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] = in_ptr[(x + 0) * Ny + y] - in_ptr[(x - 1) * Ny + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_stag_x_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_x_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void DDX_C2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny + 0; ++y) {

      result_ptr[(x)*Ny + y + 0] = in_ptr[(x)*Ny + y + 0] - in_ptr[(x)*Ny + y - 1];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_stag_y_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_y_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void DDX_C2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 0; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            in_ptr[((x + 1) * Ny + y) * Nz + z] - in_ptr[((x + 0) * Ny + y) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_x_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_x_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y + 0) * Nz + z];
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_y_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_y_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z + 0];
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] - in_ptr[((x)*Ny + y) * Nz + z + 0];
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
              in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C2_stag_z_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 2) {
    throw BoutException("Field3D AiolosMesh::DDX_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C2_stag_z_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 0; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] = in_ptr[(x + 1) * Ny + y] - in_ptr[(x + 0) * Ny + y];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_stag_x_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_x_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] = in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y + 0];
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C2_stag_y_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 2) {
    throw BoutException("Field2D AiolosMesh::DDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C2_stag_y_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C4_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x + 0) * Ny + y) * Nz + z] -
                    in_ptr[((x - 1) * Ny + y) * Nz + z]) -
             (in_ptr[((x + 1) * Ny + y) * Nz + z] -
              in_ptr[((x - 2) * Ny + y) * Nz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_x_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_x_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void DDX_C4_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (27. *
                 (in_ptr[((x)*Ny + y + 0) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z]) -
             (in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y - 2) * Nz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_y_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_y_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void DDX_C4_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])) /
              24.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz])) /
              24.;
        }
        for (int z = 2; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z - 2])) /
              24.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      in_ptr[((x)*Ny + y) * Nz + z - 1]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                in_ptr[((x)*Ny + y) * Nz + z - 2])) /
              24.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                      in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) -
               (in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)])) /
              24.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_z_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_z_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void DDX_C4_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (27. * (in_ptr[(x + 0) * Ny + y] - in_ptr[(x - 1) * Ny + y]) -
           (in_ptr[(x + 1) * Ny + y] - in_ptr[(x - 2) * Ny + y])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_stag_x_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_x_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void DDX_C4_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (27. * (in_ptr[(x)*Ny + y + 0] - in_ptr[(x)*Ny + y - 1]) -
           (in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y - 2])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_stag_y_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_y_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void DDX_C4_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (27. * (in_ptr[((x + 1) * Ny + y) * Nz + z] -
                    in_ptr[((x + 0) * Ny + y) * Nz + z]) -
             (in_ptr[((x + 2) * Ny + y) * Nz + z] -
              in_ptr[((x - 1) * Ny + y) * Nz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_x_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_x_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C4_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (27. *
                 (in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y + 0) * Nz + z]) -
             (in_ptr[((x)*Ny + y + 2) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z])) /
            24.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_y_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_y_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C4_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                      in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 2] -
                in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz])) /
              24.;
        }
        for (int z = 1; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                      in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 2] - in_ptr[((x)*Ny + y) * Nz + z - 1])) /
              24.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 1] -
                      in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                in_ptr[((x)*Ny + y) * Nz + z - 1])) /
              24.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                      in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] -
                in_ptr[((x)*Ny + y) * Nz + z - 1])) /
              24.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (27. * (in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                      in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) -
               (in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] -
                in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)])) /
              24.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D DDX_C4_stag_z_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::DDX_C4_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  DDX_C4_stag_z_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C4_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (27. * (in_ptr[(x + 1) * Ny + y] - in_ptr[(x + 0) * Ny + y]) -
           (in_ptr[(x + 2) * Ny + y] - in_ptr[(x - 1) * Ny + y])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_stag_x_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_x_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void DDX_C4_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                       const BoutReal *__restrict__ in_ptr,
                                       Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (27. * (in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y + 0]) -
           (in_ptr[(x)*Ny + y + 2] - in_ptr[(x)*Ny + y - 1])) /
          24.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D DDX_C4_stag_y_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method DDX_C4_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::DDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  DDX_C4_stag_y_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x + 1) * Ny + y) * Nz + z] + in_ptr[((x - 2) * Ny + y) * Nz + z] -
             in_ptr[((x + 0) * Ny + y) * Nz + z] - in_ptr[((x - 1) * Ny + y) * Nz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_x_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_x_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x)*Ny + y + 1) * Nz + z] + in_ptr[((x)*Ny + y - 2) * Nz + z] -
             in_ptr[((x)*Ny + y + 0) * Nz + z] - in_ptr[((x)*Ny + y - 1) * Nz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_y_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_y_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
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
              (in_ptr[((x)*Ny + y) * Nz + z + 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 0] -
               in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) /
              2.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 1] +
               in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1]) /
              2.;
        }
        for (int z = 2; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 1] + in_ptr[((x)*Ny + y) * Nz + z - 2] -
               in_ptr[((x)*Ny + y) * Nz + z + 0] - in_ptr[((x)*Ny + y) * Nz + z - 1]) /
              2.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
               in_ptr[((x)*Ny + y) * Nz + z - 2] - in_ptr[((x)*Ny + y) * Nz + z + 0] -
               in_ptr[((x)*Ny + y) * Nz + z - 1]) /
              2.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
               in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) /
              2.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_z_CtoL(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_z_Field3D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr,
                                         Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (in_ptr[(x + 1) * Ny + y] + in_ptr[(x - 2) * Ny + y] -
           in_ptr[(x + 0) * Ny + y] - in_ptr[(x - 1) * Ny + y]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_stag_x_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_x_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr,
                                         Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] = (in_ptr[(x)*Ny + y + 1] + in_ptr[(x)*Ny + y - 2] -
                                    in_ptr[(x)*Ny + y + 0] - in_ptr[(x)*Ny + y - 1]) /
                                   2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_stag_y_CtoL(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_y_Field2D_CtoL(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x + 2) * Ny + y) * Nz + z] + in_ptr[((x - 1) * Ny + y) * Nz + z] -
             in_ptr[((x + 1) * Ny + y) * Nz + z] - in_ptr[((x + 0) * Ny + y) * Nz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_x_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_x_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
            (in_ptr[((x)*Ny + y + 2) * Nz + z] + in_ptr[((x)*Ny + y - 1) * Nz + z] -
             in_ptr[((x)*Ny + y + 1) * Nz + z] - in_ptr[((x)*Ny + y + 0) * Nz + z]) /
            2.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_y_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_y_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
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
              (in_ptr[((x)*Ny + y) * Nz + z + 2] +
               in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z + 0]) /
              2.;
        }
        for (int z = 1; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 2] + in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z + 1] - in_ptr[((x)*Ny + y) * Nz + z + 0]) /
              2.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
               in_ptr[((x)*Ny + y) * Nz + z - 1] - in_ptr[((x)*Ny + y) * Nz + z + 1] -
               in_ptr[((x)*Ny + y) * Nz + z + 0]) /
              2.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz] +
               in_ptr[((x)*Ny + y) * Nz + z - 1] -
               in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
               in_ptr[((x)*Ny + y) * Nz + z + 0]) /
              2.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)] +
               in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
               in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) /
              2.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D D2DX2_C2_stag_z_LtoC(const Field3D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 4) {
    throw BoutException("Field3D AiolosMesh::D2DX2_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  D2DX2_C2_stag_z_Field3D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr,
                                         Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (in_ptr[(x + 2) * Ny + y] + in_ptr[(x - 1) * Ny + y] -
           in_ptr[(x + 1) * Ny + y] - in_ptr[(x + 0) * Ny + y]) /
          2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_stag_x_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_x_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void D2DX2_C2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr,
                                         Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] = (in_ptr[(x)*Ny + y + 2] + in_ptr[(x)*Ny + y - 1] -
                                    in_ptr[(x)*Ny + y + 1] - in_ptr[(x)*Ny + y + 0]) /
                                   2.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D D2DX2_C2_stag_y_LtoC(const Field2D &in) {
  Mesh *localmesh = in.getMesh();
  output_debug.write("Using method D2DX2_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 4) {
    throw BoutException("Field2D AiolosMesh::D2DX2_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  D2DX2_C2_stag_y_Field2D_LtoC(result_ptr, in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x - 1) * Ny + y) * Nz + z] >= 0)
                               ? v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]
                               : v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x + 0) * Ny + y) * Nz + z];
        result_ -= (v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0)
                       ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 0) * Ny + y) * Nz + z]
                       : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 1) * Ny + y) * Nz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                   (v_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                    v_in_ptr[((x - 1) * Ny + y) * Nz + z]);

        result_ptr[((x + 0) * Ny + y) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*Ny + y - 1) * Nz + z] >= 0)
                               ? v_in_ptr[((x)*Ny + y - 1) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y - 1) * Nz + z]
                               : v_in_ptr[((x)*Ny + y - 1) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y + 0) * Nz + z];
        result_ -= (v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0)
                       ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 0) * Nz + z]
                       : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 1) * Nz + z];
        result_ *= -1;
        result_ -=
            f_in_ptr[((x)*Ny + y + 0) * Nz + z] *
            (v_in_ptr[((x)*Ny + y + 0) * Nz + z] - v_in_ptr[((x)*Ny + y - 1) * Nz + z]);

        result_ptr[((x)*Ny + y + 0) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                     (v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                      v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];
          result_ *= -1;
          result_ -=
              f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 0] - v_in_ptr[((x)*Ny + y) * Nz + z - 1]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];
          result_ *= -1;
          result_ -=
              f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 0] - v_in_ptr[((x)*Ny + y) * Nz + z - 1]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] >= 0)
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]
                         : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                     (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                      v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_ = (v_in_ptr[(x - 1) * Ny + y] >= 0)
                             ? v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x - 1) * Ny + y]
                             : v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x + 0) * Ny + y];
      result_ -= (v_in_ptr[(x + 0) * Ny + y] >= 0)
                     ? v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 0) * Ny + y]
                     : v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 1) * Ny + y];
      result_ *= -1;
      result_ -= f_in_ptr[(x + 0) * Ny + y] *
                 (v_in_ptr[(x + 0) * Ny + y] - v_in_ptr[(x - 1) * Ny + y]);

      result_ptr[(x + 0) * Ny + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal result_ = (v_in_ptr[(x)*Ny + y - 1] >= 0)
                             ? v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y - 1]
                             : v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y + 0];
      result_ -= (v_in_ptr[(x)*Ny + y + 0] >= 0)
                     ? v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 0]
                     : v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 1];
      result_ *= -1;
      result_ -= f_in_ptr[(x)*Ny + y + 0] *
                 (v_in_ptr[(x)*Ny + y + 0] - v_in_ptr[(x)*Ny + y - 1]);

      result_ptr[(x)*Ny + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0)
                               ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]
                               : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x + 0) * Ny + y) * Nz + z];
        result_ -= (v_in_ptr[((x + 1) * Ny + y) * Nz + z] >= 0)
                       ? v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 0) * Ny + y) * Nz + z]
                       : v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 1) * Ny + y) * Nz + z];
        result_ *= -1;
        result_ -= f_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                   (v_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                    v_in_ptr[((x + 0) * Ny + y) * Nz + z]);

        result_ptr[((x + 0) * Ny + y) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0)
                               ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y - 1) * Nz + z]
                               : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y + 0) * Nz + z];
        result_ -= (v_in_ptr[((x)*Ny + y + 1) * Nz + z] >= 0)
                       ? v_in_ptr[((x)*Ny + y + 1) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 0) * Nz + z]
                       : v_in_ptr[((x)*Ny + y + 1) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 1) * Nz + z];
        result_ *= -1;
        result_ -=
            f_in_ptr[((x)*Ny + y + 0) * Nz + z] *
            (v_in_ptr[((x)*Ny + y + 1) * Nz + z] - v_in_ptr[((x)*Ny + y + 0) * Nz + z]);

        result_ptr[((x)*Ny + y + 0) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];
          result_ *= -1;
          result_ -=
              f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 1] - v_in_ptr[((x)*Ny + y) * Nz + z + 0]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];
          result_ *= -1;
          result_ -=
              f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
              (v_in_ptr[((x)*Ny + y) * Nz + z + 1] - v_in_ptr[((x)*Ny + y) * Nz + z + 0]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                     (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                      v_in_ptr[((x)*Ny + y) * Nz + z + 0]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0)
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]
                         : v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];
          result_ *= -1;
          result_ -= f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                     (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                      v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U1_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_U1_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U1_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_ = (v_in_ptr[(x + 0) * Ny + y] >= 0)
                             ? v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x - 1) * Ny + y]
                             : v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 0) * Ny + y];
      result_ -= (v_in_ptr[(x + 1) * Ny + y] >= 0)
                     ? v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 0) * Ny + y]
                     : v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 1) * Ny + y];
      result_ *= -1;
      result_ -= f_in_ptr[(x + 0) * Ny + y] *
                 (v_in_ptr[(x + 1) * Ny + y] - v_in_ptr[(x + 0) * Ny + y]);

      result_ptr[(x + 0) * Ny + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U1_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal result_ = (v_in_ptr[(x)*Ny + y + 0] >= 0)
                             ? v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y - 1]
                             : v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 0];
      result_ -= (v_in_ptr[(x)*Ny + y + 1] >= 0)
                     ? v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 0]
                     : v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 1];
      result_ *= -1;
      result_ -= f_in_ptr[(x)*Ny + y + 0] *
                 (v_in_ptr[(x)*Ny + y + 1] - v_in_ptr[(x)*Ny + y + 0]);

      result_ptr[(x)*Ny + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U1_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U1_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U1_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x + 0) * Ny + y) * Nz + z] > 0 &&
            v_in_ptr[((x - 1) * Ny + y) * Nz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x - 1) * Ny + y) * Nz + z] -
                     .5 * v_in_ptr[((x - 2) * Ny + y) * Nz + z]) *
                    (.5 * f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
                     2. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                     1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z]);
        } else if (v_in_ptr[((x + 0) * Ny + y) * Nz + z] < 0 &&
                   v_in_ptr[((x - 1) * Ny + y) * Nz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                     .5 * v_in_ptr[((x + 1) * Ny + y) * Nz + z]) *
                    (-1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                     2. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                     .5 * f_in_ptr[((x + 2) * Ny + y) * Nz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                           v_in_ptr[((x - 1) * Ny + y) * Nz + z]) *
                    (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
        }

        result_ptr[((x + 0) * Ny + y) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x)*Ny + y + 0) * Nz + z] > 0 &&
            v_in_ptr[((x)*Ny + y - 1) * Nz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x)*Ny + y - 1) * Nz + z] -
                     .5 * v_in_ptr[((x)*Ny + y - 2) * Nz + z]) *
                    (.5 * f_in_ptr[((x)*Ny + y - 2) * Nz + z] -
                     2. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                     1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z]);
        } else if (v_in_ptr[((x)*Ny + y + 0) * Nz + z] < 0 &&
                   v_in_ptr[((x)*Ny + y - 1) * Nz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                     .5 * v_in_ptr[((x)*Ny + y + 1) * Nz + z]) *
                    (-1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                     2. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                     .5 * f_in_ptr[((x)*Ny + y + 2) * Nz + z]);
        } else {
          result_ =
              .25 * (v_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                     v_in_ptr[((x)*Ny + y - 1) * Nz + z]) *
              (f_in_ptr[((x)*Ny + y + 1) * Nz + z] - f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
        }

        result_ptr[((x)*Ny + y + 0) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                             v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                             v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        for (int z = 2; z < Nz - 2; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 2]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                             v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 2;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 2]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                             v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z - 1] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 2]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                             v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                             v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x + 0) * Ny + y] > 0 && v_in_ptr[(x - 1) * Ny + y] > 0) {
        result_ = (1.5 * v_in_ptr[(x - 1) * Ny + y] - .5 * v_in_ptr[(x - 2) * Ny + y]) *
                  (.5 * f_in_ptr[(x - 2) * Ny + y] - 2. * f_in_ptr[(x - 1) * Ny + y] +
                   1.5 * f_in_ptr[(x + 0) * Ny + y]);
      } else if (v_in_ptr[(x + 0) * Ny + y] < 0 && v_in_ptr[(x - 1) * Ny + y] < 0) {
        result_ = (1.5 * v_in_ptr[(x + 0) * Ny + y] - .5 * v_in_ptr[(x + 1) * Ny + y]) *
                  (-1.5 * f_in_ptr[(x + 0) * Ny + y] + 2. * f_in_ptr[(x + 1) * Ny + y] -
                   .5 * f_in_ptr[(x + 2) * Ny + y]);
      } else {
        result_ = .25 * (v_in_ptr[(x + 0) * Ny + y] + v_in_ptr[(x - 1) * Ny + y]) *
                  (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]);
      }

      result_ptr[(x + 0) * Ny + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x)*Ny + y + 0] > 0 && v_in_ptr[(x)*Ny + y - 1] > 0) {
        result_ = (1.5 * v_in_ptr[(x)*Ny + y - 1] - .5 * v_in_ptr[(x)*Ny + y - 2]) *
                  (.5 * f_in_ptr[(x)*Ny + y - 2] - 2. * f_in_ptr[(x)*Ny + y - 1] +
                   1.5 * f_in_ptr[(x)*Ny + y + 0]);
      } else if (v_in_ptr[(x)*Ny + y + 0] < 0 && v_in_ptr[(x)*Ny + y - 1] < 0) {
        result_ = (1.5 * v_in_ptr[(x)*Ny + y + 0] - .5 * v_in_ptr[(x)*Ny + y + 1]) *
                  (-1.5 * f_in_ptr[(x)*Ny + y + 0] + 2. * f_in_ptr[(x)*Ny + y + 1] -
                   .5 * f_in_ptr[(x)*Ny + y + 2]);
      } else {
        result_ = .25 * (v_in_ptr[(x)*Ny + y + 0] + v_in_ptr[(x)*Ny + y - 1]) *
                  (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]);
      }

      result_ptr[(x)*Ny + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x + 1) * Ny + y) * Nz + z] > 0 &&
            v_in_ptr[((x + 0) * Ny + y) * Nz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x + 0) * Ny + y) * Nz + z] -
                     .5 * v_in_ptr[((x - 1) * Ny + y) * Nz + z]) *
                    (.5 * f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
                     2. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                     1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z]);
        } else if (v_in_ptr[((x + 1) * Ny + y) * Nz + z] < 0 &&
                   v_in_ptr[((x + 0) * Ny + y) * Nz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                     .5 * v_in_ptr[((x + 2) * Ny + y) * Nz + z]) *
                    (-1.5 * f_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                     2. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                     .5 * f_in_ptr[((x + 2) * Ny + y) * Nz + z]);
        } else {
          result_ = .25 * (v_in_ptr[((x + 1) * Ny + y) * Nz + z] +
                           v_in_ptr[((x + 0) * Ny + y) * Nz + z]) *
                    (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
        }

        result_ptr[((x + 0) * Ny + y) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_;
        if (v_in_ptr[((x)*Ny + y + 1) * Nz + z] > 0 &&
            v_in_ptr[((x)*Ny + y + 0) * Nz + z] > 0) {
          result_ = (1.5 * v_in_ptr[((x)*Ny + y + 0) * Nz + z] -
                     .5 * v_in_ptr[((x)*Ny + y - 1) * Nz + z]) *
                    (.5 * f_in_ptr[((x)*Ny + y - 2) * Nz + z] -
                     2. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                     1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z]);
        } else if (v_in_ptr[((x)*Ny + y + 1) * Nz + z] < 0 &&
                   v_in_ptr[((x)*Ny + y + 0) * Nz + z] < 0) {
          result_ = (1.5 * v_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                     .5 * v_in_ptr[((x)*Ny + y + 2) * Nz + z]) *
                    (-1.5 * f_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                     2. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
                     .5 * f_in_ptr[((x)*Ny + y + 2) * Nz + z]);
        } else {
          result_ =
              .25 * (v_in_ptr[((x)*Ny + y + 1) * Nz + z] +
                     v_in_ptr[((x)*Ny + y + 0) * Nz + z]) *
              (f_in_ptr[((x)*Ny + y + 1) * Nz + z] - f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
        }

        result_ptr[((x)*Ny + y + 0) * Nz + z] = result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                             v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = 1;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                             v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        for (int z = 2; z < Nz - 2; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 2]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                             v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 2;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 1] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                             v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + z + 0] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 0] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                             v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_;
          if (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] > 0 &&
              v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] > 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) *
                      (.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
                       2. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                       1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]);
          } else if (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] < 0 &&
                     v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] < 0) {
            result_ = (1.5 * v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       .5 * v_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) *
                      (-1.5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                       2. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       .5 * f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]);
          } else {
            result_ = .25 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
                             v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) *
                      (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                       f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
          }

          result_ptr[((x)*Ny + y) * Nz + z + 0] = result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_U2_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_U2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_U2_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x + 1) * Ny + y] > 0 && v_in_ptr[(x + 0) * Ny + y] > 0) {
        result_ = (1.5 * v_in_ptr[(x + 0) * Ny + y] - .5 * v_in_ptr[(x - 1) * Ny + y]) *
                  (.5 * f_in_ptr[(x - 2) * Ny + y] - 2. * f_in_ptr[(x - 1) * Ny + y] +
                   1.5 * f_in_ptr[(x + 0) * Ny + y]);
      } else if (v_in_ptr[(x + 1) * Ny + y] < 0 && v_in_ptr[(x + 0) * Ny + y] < 0) {
        result_ = (1.5 * v_in_ptr[(x + 1) * Ny + y] - .5 * v_in_ptr[(x + 2) * Ny + y]) *
                  (-1.5 * f_in_ptr[(x + 0) * Ny + y] + 2. * f_in_ptr[(x + 1) * Ny + y] -
                   .5 * f_in_ptr[(x + 2) * Ny + y]);
      } else {
        result_ = .25 * (v_in_ptr[(x + 1) * Ny + y] + v_in_ptr[(x + 0) * Ny + y]) *
                  (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]);
      }

      result_ptr[(x + 0) * Ny + y] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_U2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {
      BoutReal result_;
      if (v_in_ptr[(x)*Ny + y + 1] > 0 && v_in_ptr[(x)*Ny + y + 0] > 0) {
        result_ = (1.5 * v_in_ptr[(x)*Ny + y + 0] - .5 * v_in_ptr[(x)*Ny + y - 1]) *
                  (.5 * f_in_ptr[(x)*Ny + y - 2] - 2. * f_in_ptr[(x)*Ny + y - 1] +
                   1.5 * f_in_ptr[(x)*Ny + y + 0]);
      } else if (v_in_ptr[(x)*Ny + y + 1] < 0 && v_in_ptr[(x)*Ny + y + 0] < 0) {
        result_ = (1.5 * v_in_ptr[(x)*Ny + y + 1] - .5 * v_in_ptr[(x)*Ny + y + 2]) *
                  (-1.5 * f_in_ptr[(x)*Ny + y + 0] + 2. * f_in_ptr[(x)*Ny + y + 1] -
                   .5 * f_in_ptr[(x)*Ny + y + 2]);
      } else {
        result_ = .25 * (v_in_ptr[(x)*Ny + y + 1] + v_in_ptr[(x)*Ny + y + 0]) *
                  (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]);
      }

      result_ptr[(x)*Ny + y + 0] = result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_U2_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_U2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_U2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_U2_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            0.5 * (v_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                   v_in_ptr[((x - 1) * Ny + y) * Nz + z]) *
            0.5 * (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                   f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            0.5 *
            (v_in_ptr[((x)*Ny + y + 0) * Nz + z] + v_in_ptr[((x)*Ny + y - 1) * Nz + z]) *
            0.5 *
            (f_in_ptr[((x)*Ny + y + 1) * Nz + z] - f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
              0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + z + 1] - f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z - 1]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                     v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          0.5 * (v_in_ptr[(x + 0) * Ny + y] + v_in_ptr[(x - 1) * Ny + y]) * 0.5 *
          (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          0.5 * (v_in_ptr[(x)*Ny + y + 0] + v_in_ptr[(x)*Ny + y - 1]) * 0.5 *
          (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            0.5 * (v_in_ptr[((x + 1) * Ny + y) * Nz + z] +
                   v_in_ptr[((x + 0) * Ny + y) * Nz + z]) *
            0.5 * (f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                   f_in_ptr[((x - 1) * Ny + y) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            0.5 *
            (v_in_ptr[((x)*Ny + y + 1) * Nz + z] + v_in_ptr[((x)*Ny + y + 0) * Nz + z]) *
            0.5 *
            (f_in_ptr[((x)*Ny + y + 1) * Nz + z] - f_in_ptr[((x)*Ny + y - 1) * Nz + z]);
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]);
        }
        for (int z = 1; z < Nz - 1; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
              0.5 *
              (f_in_ptr[((x)*Ny + y) * Nz + z + 1] - f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z - 1]);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              0.5 * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] +
                     v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) *
              0.5 * (f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]);
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C2_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::VDDX_C2_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C2_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          0.5 * (v_in_ptr[(x + 1) * Ny + y] + v_in_ptr[(x + 0) * Ny + y]) * 0.5 *
          (f_in_ptr[(x + 1) * Ny + y] - f_in_ptr[(x - 1) * Ny + y]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C2_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          0.5 * (v_in_ptr[(x)*Ny + y + 1] + v_in_ptr[(x)*Ny + y + 0]) * 0.5 *
          (f_in_ptr[(x)*Ny + y + 1] - f_in_ptr[(x)*Ny + y - 1]);
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C2_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C2_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::VDDX_C2_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C2_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            (9. * (v_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                   v_in_ptr[((x + 0) * Ny + y) * Nz + z]) -
             v_in_ptr[((x - 2) * Ny + y) * Nz + z] -
             v_in_ptr[((x + 1) * Ny + y) * Nz + z]) /
            16. * (8. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                   8. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                   f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
                   f_in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            (9. * (v_in_ptr[((x)*Ny + y - 1) * Nz + z] +
                   v_in_ptr[((x)*Ny + y + 0) * Nz + z]) -
             v_in_ptr[((x)*Ny + y - 2) * Nz + z] - v_in_ptr[((x)*Ny + y + 1) * Nz + z]) /
            16. *
            (8. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
             f_in_ptr[((x)*Ny + y - 2) * Nz + z] - f_in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 1]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 0]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 2] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                     v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]) -
               v_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
               v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                     f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
                     f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (9. * (v_in_ptr[(x - 1) * Ny + y] + v_in_ptr[(x + 0) * Ny + y]) -
           v_in_ptr[(x - 2) * Ny + y] - v_in_ptr[(x + 1) * Ny + y]) /
          16. * (8. * f_in_ptr[(x + 1) * Ny + y] - 8. * f_in_ptr[(x - 1) * Ny + y] +
                 f_in_ptr[(x - 2) * Ny + y] - f_in_ptr[(x + 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (9. * (v_in_ptr[(x)*Ny + y - 1] + v_in_ptr[(x)*Ny + y + 0]) -
           v_in_ptr[(x)*Ny + y - 2] - v_in_ptr[(x)*Ny + y + 1]) /
          16. * (8. * f_in_ptr[(x)*Ny + y + 1] - 8. * f_in_ptr[(x)*Ny + y - 1] +
                 f_in_ptr[(x)*Ny + y - 2] - f_in_ptr[(x)*Ny + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x + 0) * Ny + y) * Nz + z] =
            (9. * (v_in_ptr[((x + 0) * Ny + y) * Nz + z] +
                   v_in_ptr[((x + 1) * Ny + y) * Nz + z]) -
             v_in_ptr[((x - 1) * Ny + y) * Nz + z] -
             v_in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            16. * (8. * f_in_ptr[((x + 1) * Ny + y) * Nz + z] -
                   8. * f_in_ptr[((x - 1) * Ny + y) * Nz + z] +
                   f_in_ptr[((x - 2) * Ny + y) * Nz + z] -
                   f_in_ptr[((x + 2) * Ny + y) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
    for (int y = 2; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {

        result_ptr[((x)*Ny + y + 0) * Nz + z] =
            (9. * (v_in_ptr[((x)*Ny + y + 0) * Nz + z] +
                   v_in_ptr[((x)*Ny + y + 1) * Nz + z]) -
             v_in_ptr[((x)*Ny + y - 1) * Nz + z] - v_in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            16. *
            (8. * f_in_ptr[((x)*Ny + y + 1) * Nz + z] -
             8. * f_in_ptr[((x)*Ny + y - 1) * Nz + z] +
             f_in_ptr[((x)*Ny + y - 2) * Nz + z] - f_in_ptr[((x)*Ny + y + 2) * Nz + z]) /
            12.;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
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
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 1]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 1]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2 + Nz] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        for (int z = 2; z < Nz - 2; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 1]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2]) /
              12.;
        }
        {
          int z = Nz - 2;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 1]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
        {
          int z = Nz - 1;

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + z + 0] +
                     v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz]) -
               v_in_ptr[((x)*Ny + y) * Nz + z - 1] -
               v_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + z - 1] +
                     f_in_ptr[((x)*Ny + y) * Nz + z - 2] -
                     f_in_ptr[((x)*Ny + y) * Nz + z + 2 - Nz]) /
              12.;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result_ptr[((x)*Ny + y) * Nz + z + 0] =
              (9. * (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] +
                     v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)]) -
               v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] -
               v_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              16. * (8. * f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] -
                     8. * f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] +
                     f_in_ptr[((x)*Ny + y) * Nz + +((z - 2 + 2 * Nz) % Nz)] -
                     f_in_ptr[((x)*Ny + y) * Nz + +((z + 2) % Nz)]) /
              12.;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D VDDX_C4_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 5) {
    throw BoutException("Field3D AiolosMesh::VDDX_C4_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  VDDX_C4_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {

      result_ptr[(x + 0) * Ny + y] =
          (9. * (v_in_ptr[(x + 0) * Ny + y] + v_in_ptr[(x + 1) * Ny + y]) -
           v_in_ptr[(x - 1) * Ny + y] - v_in_ptr[(x + 2) * Ny + y]) /
          16. * (8. * f_in_ptr[(x + 1) * Ny + y] - 8. * f_in_ptr[(x - 1) * Ny + y] +
                 f_in_ptr[(x - 2) * Ny + y] - f_in_ptr[(x + 2) * Ny + y]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void VDDX_C4_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 2; ++y) {

      result_ptr[(x)*Ny + y + 0] =
          (9. * (v_in_ptr[(x)*Ny + y + 0] + v_in_ptr[(x)*Ny + y + 1]) -
           v_in_ptr[(x)*Ny + y - 1] - v_in_ptr[(x)*Ny + y + 2]) /
          16. * (8. * f_in_ptr[(x)*Ny + y + 1] - 8. * f_in_ptr[(x)*Ny + y - 1] +
                 f_in_ptr[(x)*Ny + y - 2] - f_in_ptr[(x)*Ny + y + 2]) /
          12.;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D VDDX_C4_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method VDDX_C4_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 5) {
    throw BoutException("Field2D AiolosMesh::VDDX_C4_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  VDDX_C4_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_x_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x - 1) * Ny + y) * Nz + z] >= 0)
                               ? v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]
                               : v_in_ptr[((x - 1) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x + 0) * Ny + y) * Nz + z];
        result_ -= (v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0)
                       ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 0) * Ny + y) * Nz + z]
                       : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 1) * Ny + y) * Nz + z];

        result_ptr[((x + 0) * Ny + y) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_x_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_x_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_y_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*Ny + y - 1) * Nz + z] >= 0)
                               ? v_in_ptr[((x)*Ny + y - 1) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y - 1) * Nz + z]
                               : v_in_ptr[((x)*Ny + y - 1) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y + 0) * Nz + z];
        result_ -= (v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0)
                       ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 0) * Nz + z]
                       : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 1) * Nz + z];

        result_ptr[((x)*Ny + y + 0) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_y_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_y_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_z_Field3D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z - 1] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z - 1] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] >= 0)
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]
                         : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_z_CtoL(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_z_Field3D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_z_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_z_Field3D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_x_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_ = (v_in_ptr[(x - 1) * Ny + y] >= 0)
                             ? v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x - 1) * Ny + y]
                             : v_in_ptr[(x - 1) * Ny + y] * f_in_ptr[(x + 0) * Ny + y];
      result_ -= (v_in_ptr[(x + 0) * Ny + y] >= 0)
                     ? v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 0) * Ny + y]
                     : v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 1) * Ny + y];

      result_ptr[(x + 0) * Ny + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_stag_x_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_x_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_x_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_y_Field2D_CtoL(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal result_ = (v_in_ptr[(x)*Ny + y - 1] >= 0)
                             ? v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y - 1]
                             : v_in_ptr[(x)*Ny + y - 1] * f_in_ptr[(x)*Ny + y + 0];
      result_ -= (v_in_ptr[(x)*Ny + y + 0] >= 0)
                     ? v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 0]
                     : v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 1];

      result_ptr[(x)*Ny + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_stag_y_CtoL(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field2D_CtoL!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_y_CtoL - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_y_Field2D_CtoL(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_x_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x + 0) * Ny + y) * Nz + z] >= 0)
                               ? v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x - 1) * Ny + y) * Nz + z]
                               : v_in_ptr[((x + 0) * Ny + y) * Nz + z] *
                                     f_in_ptr[((x + 0) * Ny + y) * Nz + z];
        result_ -= (v_in_ptr[((x + 1) * Ny + y) * Nz + z] >= 0)
                       ? v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 0) * Ny + y) * Nz + z]
                       : v_in_ptr[((x + 1) * Ny + y) * Nz + z] *
                             f_in_ptr[((x + 1) * Ny + y) * Nz + z];

        result_ptr[((x + 0) * Ny + y) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_x_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_x_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_y_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        BoutReal result_ = (v_in_ptr[((x)*Ny + y + 0) * Nz + z] >= 0)
                               ? v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y - 1) * Nz + z]
                               : v_in_ptr[((x)*Ny + y + 0) * Nz + z] *
                                     f_in_ptr[((x)*Ny + y + 0) * Nz + z];
        result_ -= (v_in_ptr[((x)*Ny + y + 1) * Nz + z] >= 0)
                       ? v_in_ptr[((x)*Ny + y + 1) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 0) * Nz + z]
                       : v_in_ptr[((x)*Ny + y + 1) * Nz + z] *
                             f_in_ptr[((x)*Ny + y + 1) * Nz + z];

        result_ptr[((x)*Ny + y + 0) * Nz + z] = -result_;
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_y_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_y_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_z_Field3D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
  const int Nz = localmesh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1 + Nz]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        for (int z = 1; z < Nz - 1; ++z) {
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
        {
          int z = Nz - 1;
          BoutReal result_ = (v_in_ptr[((x)*Ny + y) * Nz + z + 0] >= 0)
                                 ? v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z - 1]
                                 : v_in_ptr[((x)*Ny + y) * Nz + z + 0] *
                                       f_in_ptr[((x)*Ny + y) * Nz + z + 0];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 0]
                         : v_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz] *
                               f_in_ptr[((x)*Ny + y) * Nz + z + 1 - Nz];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          BoutReal result_ =
              (v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] >= 0)
                  ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z - 1 + 1 * Nz) % Nz)]
                  : v_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)] *
                        f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)];
          result_ -= (v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] >= 0)
                         ? v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 0) % Nz)]
                         : v_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)] *
                               f_in_ptr[((x)*Ny + y) * Nz + +((z + 1) % Nz)];

          result_ptr[((x)*Ny + y) * Nz + z + 0] = -result_;
        }
      }
    }
  }
}

// This file is auto-generated - do not edit!
Field3D FDDX_U1_stag_z_LtoC(const Field3D &v_in, const Field3D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_z_Field3D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNz < 3) {
    throw BoutException("Field3D AiolosMesh::FDDX_U1_stag_z_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field3D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  FDDX_U1_stag_z_Field3D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_x_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->xstart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      BoutReal result_ = (v_in_ptr[(x + 0) * Ny + y] >= 0)
                             ? v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x - 1) * Ny + y]
                             : v_in_ptr[(x + 0) * Ny + y] * f_in_ptr[(x + 0) * Ny + y];
      result_ -= (v_in_ptr[(x + 1) * Ny + y] >= 0)
                     ? v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 0) * Ny + y]
                     : v_in_ptr[(x + 1) * Ny + y] * f_in_ptr[(x + 1) * Ny + y];

      result_ptr[(x + 0) * Ny + y] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_stag_x_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_x_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNx < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_x_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_x_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

static void FDDX_U1_stag_y_Field2D_LtoC(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ v_in_ptr,
                                        const BoutReal *__restrict__ f_in_ptr,
                                        Mesh *localmesh) {
  const int Nx = localmesh->LocalNx;
  const int Ny = localmesh->LocalNy;
#if CHECK > 0
  if (localmesh->ystart < 1) {
    throw BoutException(
        "Cannot compute derivative - need at least 1 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 1; ++y) {
      BoutReal result_ = (v_in_ptr[(x)*Ny + y + 0] >= 0)
                             ? v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y - 1]
                             : v_in_ptr[(x)*Ny + y + 0] * f_in_ptr[(x)*Ny + y + 0];
      result_ -= (v_in_ptr[(x)*Ny + y + 1] >= 0)
                     ? v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 0]
                     : v_in_ptr[(x)*Ny + y + 1] * f_in_ptr[(x)*Ny + y + 1];

      result_ptr[(x)*Ny + y + 0] = -result_;
    }
  }
}

// This file is auto-generated - do not edit!
Field2D FDDX_U1_stag_y_LtoC(const Field2D &v_in, const Field2D &f_in) {
  Mesh *localmesh = v_in.getMesh();
  ASSERT1(localmesh == f_in.getMesh());
  output_debug.write("Using method FDDX_U1_stag_y_Field2D_LtoC!\n");
#if CHECK > 0
  if (localmesh->LocalNy < 3) {
    throw BoutException("Field2D AiolosMesh::FDDX_U1_stag_y_LtoC - Not enough guards "
                        "cells to take derivative!");
  }
#endif
  Field2D result(localmesh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  FDDX_U1_stag_y_Field2D_LtoC(result_ptr, v_in_ptr, f_in_ptr, localmesh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}
