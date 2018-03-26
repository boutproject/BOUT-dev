static void interp_to_on_Field3D_x(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < Nx - 1; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result(x + 0, y, z) =
            -6.25000e-02 * in(x - 2, y, z) + 5.62500e-01 * in(x - 1, y, z) +
            5.62500e-01 * in(x + 0, y, z) - 6.25000e-02 * in(x + 1, y, z);
      }
    }
  }
  int x = 1;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result(x + 0, y, z) = +3.12500e-01 * in(x - 1, y, z) +
                            9.37500e-01 * in(x + 0, y, z) -
                            3.12500e-01 * in(x + 1, y, z) + 6.25000e-02 * in(x + 2, y, z);
      result(x - 1, y, z) = +2.18750e+00 * in(x - 1, y, z) -
                            2.18750e+00 * in(x + 0, y, z) +
                            1.31250e+00 * in(x + 1, y, z) - 3.12500e-01 * in(x + 2, y, z);
    }
  }
  x = Nx - 1;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result(x + 0, y, z) = +6.25000e-02 * in(x - 3, y, z) -
                            3.12500e-01 * in(x - 2, y, z) +
                            9.37500e-01 * in(x - 1, y, z) + 3.12500e-01 * in(x + 0, y, z);
    }
  }
}
static void interp_to_on_Field3D_y(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 2; y < Ny - 1; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result(x, y + 0, z) =
            -6.25000e-02 * in(x, y - 2, z) + 5.62500e-01 * in(x, y - 1, z) +
            5.62500e-01 * in(x, y + 0, z) - 6.25000e-02 * in(x, y + 1, z);
      }
    }
  }
  int y = 1;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result(x, y + 0, z) = +3.12500e-01 * in(x, y - 1, z) +
                            9.37500e-01 * in(x, y + 0, z) -
                            3.12500e-01 * in(x, y + 1, z) + 6.25000e-02 * in(x, y + 2, z);
      result(x, y - 1, z) = +2.18750e+00 * in(x, y - 1, z) -
                            2.18750e+00 * in(x, y + 0, z) +
                            1.31250e+00 * in(x, y + 1, z) - 3.12500e-01 * in(x, y + 2, z);
    }
  }
  y = Ny - 1;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result(x, y + 0, z) = +6.25000e-02 * in(x, y - 3, z) -
                            3.12500e-01 * in(x, y - 2, z) +
                            9.37500e-01 * in(x, y - 1, z) + 3.12500e-01 * in(x, y + 0, z);
    }
  }
}
static void interp_to_on_Field3D_z(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 2 + Nz) + 5.62500e-01 * in(x, y, z - 1 + Nz) +
              5.62500e-01 * in(x, y, z + 0) - 6.25000e-02 * in(x, y, z + 1);
        }
        {
          int z = 1;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 2 + Nz) + 5.62500e-01 * in(x, y, z - 1) +
              5.62500e-01 * in(x, y, z + 0) - 6.25000e-02 * in(x, y, z + 1);
        }
        for (int z = 2; z < Nz - 1; ++z) {

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 2) + 5.62500e-01 * in(x, y, z - 1) +
              5.62500e-01 * in(x, y, z + 0) - 6.25000e-02 * in(x, y, z + 1);
        }
        {
          int z = Nz - 1;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 2) + 5.62500e-01 * in(x, y, z - 1) +
              5.62500e-01 * in(x, y, z + 0) - 6.25000e-02 * in(x, y, z + 1 - Nz);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result(x, y, z + 0) = -6.25000e-02 * in(x, y, +((z - 2 + 2 * Nz) % Nz)) +
                                5.62500e-01 * in(x, y, +((z - 1 + 1 * Nz) % Nz)) +
                                5.62500e-01 * in(x, y, +((z + 0) % Nz)) -
                                6.25000e-02 * in(x, y, +((z + 1) % Nz));
        }
      }
    }
  }
}
static void interp_to_off_Field3D_x(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < Nx - 2; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result(x + 0, y, z) =
            -6.25000e-02 * in(x - 1, y, z) + 5.62500e-01 * in(x + 0, y, z) +
            5.62500e-01 * in(x + 1, y, z) - 6.25000e-02 * in(x + 2, y, z);
      }
    }
  }
  int x = 0;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result(x + 0, y, z) = +3.12500e-01 * in(x + 0, y, z) +
                            9.37500e-01 * in(x + 1, y, z) -
                            3.12500e-01 * in(x + 2, y, z) + 6.25000e-02 * in(x + 3, y, z);
    }
  }
  x = Nx - 2;
  for (int y = 0; y < Ny; ++y) {
    for (int z = 0; z < Nz; ++z) {
      result(x + 0, y, z) = +6.25000e-02 * in(x - 2, y, z) -
                            3.12500e-01 * in(x - 1, y, z) +
                            9.37500e-01 * in(x + 0, y, z) + 3.12500e-01 * in(x + 1, y, z);
      result(x + 1, y, z) = -3.12500e-01 * in(x - 2, y, z) +
                            1.31250e+00 * in(x - 1, y, z) -
                            2.18750e+00 * in(x + 0, y, z) + 2.18750e+00 * in(x + 1, y, z);
    }
  }
}
static void interp_to_off_Field3D_y(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < Nx; ++x) {
    for (int y = 1; y < Ny - 2; ++y) {
      for (int z = 0; z < Nz; ++z) {
        result(x, y + 0, z) =
            -6.25000e-02 * in(x, y - 1, z) + 5.62500e-01 * in(x, y + 0, z) +
            5.62500e-01 * in(x, y + 1, z) - 6.25000e-02 * in(x, y + 2, z);
      }
    }
  }
  int y = 0;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result(x, y + 0, z) = +3.12500e-01 * in(x, y + 0, z) +
                            9.37500e-01 * in(x, y + 1, z) -
                            3.12500e-01 * in(x, y + 2, z) + 6.25000e-02 * in(x, y + 3, z);
    }
  }
  y = Ny - 2;
  for (int x = 0; x < Nx; ++x) {
    for (int z = 0; z < Nz; ++z) {
      result(x, y + 0, z) = +6.25000e-02 * in(x, y - 2, z) -
                            3.12500e-01 * in(x, y - 1, z) +
                            9.37500e-01 * in(x, y + 0, z) + 3.12500e-01 * in(x, y + 1, z);
      result(x, y + 1, z) = -3.12500e-01 * in(x, y - 2, z) +
                            1.31250e+00 * in(x, y - 1, z) -
                            2.18750e+00 * in(x, y + 0, z) + 2.18750e+00 * in(x, y + 1, z);
    }
  }
}
static void interp_to_off_Field3D_z(Field3D &result, const Field3D &in, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
  if (Nz > 3) {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        {
          int z = 0;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 1 + Nz) + 5.62500e-01 * in(x, y, z + 0) +
              5.62500e-01 * in(x, y, z + 1) - 6.25000e-02 * in(x, y, z + 2);
        }
        for (int z = 1; z < Nz - 2; ++z) {

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 1) + 5.62500e-01 * in(x, y, z + 0) +
              5.62500e-01 * in(x, y, z + 1) - 6.25000e-02 * in(x, y, z + 2);
        }
        {
          int z = Nz - 2;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 1) + 5.62500e-01 * in(x, y, z + 0) +
              5.62500e-01 * in(x, y, z + 1) - 6.25000e-02 * in(x, y, z + 2 - Nz);
        }
        {
          int z = Nz - 1;

          result(x, y, z + 0) =
              -6.25000e-02 * in(x, y, z - 1) + 5.62500e-01 * in(x, y, z + 0) +
              5.62500e-01 * in(x, y, z + 1 - Nz) - 6.25000e-02 * in(x, y, z + 2 - Nz);
        }
      }
    }
  } else {
    for (int x = 0; x < Nx; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {

          result(x, y, z + 0) = -6.25000e-02 * in(x, y, +((z - 1 + 1 * Nz) % Nz)) +
                                5.62500e-01 * in(x, y, +((z + 0) % Nz)) +
                                5.62500e-01 * in(x, y, +((z + 1) % Nz)) -
                                6.25000e-02 * in(x, y, +((z + 2) % Nz));
        }
      }
    }
  }
}
const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc) const {
  Mesh *msh = f.getMesh();
  Field3D result(msh);
  result.allocate();
  if (f.getLocation() != CELL_CENTRE) {
    // we first need to go back to centre before we can go anywhere else
    switch (f.getLocation()) {
    case CELL_XLOW:
      interp_to_off_Field3D_x(result, f, msh);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc);
      break;
    case CELL_YLOW:
      interp_to_off_Field3D_y(result, f, msh);
      result.setLocation(CELL_CENTRE);
      // return or interpolate again
      return interp_to(result, loc);
      break;
    case CELL_ZLOW:
      interp_to_off_Field3D_z(result, f, msh);
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
    interp_to_on_Field3D_x(result, f, msh);
    result.setLocation(CELL_XLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  case CELL_YLOW:
    interp_to_on_Field3D_y(result, f, msh);
    result.setLocation(CELL_YLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  case CELL_ZLOW:
    interp_to_on_Field3D_z(result, f, msh);
    result.setLocation(CELL_ZLOW);
    // return or interpolate again
    return interp_to(result, loc);
    break;
  default:
    throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                        strLocation(loc));
  }
}
const BoutReal WENO_SMALL = 1.0e-8;

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexDDX_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_W2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexDDX_norm_DIFF_W2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_W2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_W2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_norm_DIFF_W2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexDDX_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_S2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexDDX_norm_DIFF_S2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_S2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_S2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_norm_DIFF_S2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexDDY_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_W2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexDDY_norm_DIFF_W2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_W2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_W2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_norm_DIFF_W2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexDDY_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_S2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexDDY_norm_DIFF_S2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_S2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_S2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_norm_DIFF_S2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_norm_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_norm_DIFF_W2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_norm_DIFF_W2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_norm_DIFF_W2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_norm_DIFF_W2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_norm_DIFF_W2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_norm_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_norm_DIFF_S2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_norm_DIFF_S2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_norm_DIFF_S2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_norm_DIFF_S2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_norm_DIFF_S2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexDDX_norm_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_norm_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_W2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexDDX_norm_DIFF_W2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_W2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_W2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_norm_DIFF_W2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexDDX_norm_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_norm_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_norm_DIFF_S2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexDDX_norm_DIFF_S2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_norm_DIFF_S2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_norm_DIFF_S2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_norm_DIFF_S2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexDDY_norm_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_norm_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_W2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexDDY_norm_DIFF_W2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_W2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_W2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_norm_DIFF_W2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexDDY_norm_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_norm_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_norm_DIFF_S2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexDDY_norm_DIFF_S2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_norm_DIFF_S2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_norm_DIFF_S2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_norm_DIFF_S2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexD2DX2_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DX2_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexD2DX2_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DX2_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexD2DY2_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DY2_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexD2DY2_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DY2_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DZ2_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexD2DZ2_norm_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DZ2_norm_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DZ2_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DZ2_norm_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DZ2_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexD2DZ2_norm_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DZ2_norm_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DZ2_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DZ2_norm_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexD2DX2_norm_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DX2_norm_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexD2DX2_norm_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DX2_norm_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexD2DY2_norm_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DY2_norm_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                            const BoutReal *__restrict__ in_ptr,
                                            Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexD2DY2_norm_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DY2_norm_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_norm_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_norm_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U3_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_norm_DIFF_U3(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U3!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U3 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_norm_DIFF_U3_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_norm_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_norm_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U3_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_norm_DIFF_U3(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U3!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U3 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_norm_DIFF_U3_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_norm_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_norm_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_norm_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_norm_DIFF_U2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_norm_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_norm_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_norm_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_norm_DIFF_U3_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_norm_DIFF_U3(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_norm_DIFF_U3!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_norm_DIFF_U3 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_norm_DIFF_U3_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_norm_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_norm_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_norm_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_norm_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_norm_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_norm_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_norm_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_U3_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_norm_DIFF_U3(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_U3!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_U3 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_norm_DIFF_U3_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_norm_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_norm_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_norm_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_norm_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_norm_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_norm_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_norm_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_norm_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_U3_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_norm_DIFF_U3(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_U3!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_U3 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_norm_DIFF_U3_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_norm_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_norm_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexFDDX_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDX_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexFDDX_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDX_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexFDDX_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDX_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexFDDY_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDY_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexFDDY_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDY_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexFDDY_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDY_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDZ_norm_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexFDDZ_norm_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDZ_norm_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDZ_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDZ_norm_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDZ_norm_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexFDDZ_norm_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDZ_norm_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDZ_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDZ_norm_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDZ_norm_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexFDDZ_norm_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDZ_norm_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDZ_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDZ_norm_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexFDDX_norm_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDX_norm_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexFDDX_norm_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDX_norm_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexFDDX_norm_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDX_norm_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexFDDY_norm_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDY_norm_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexFDDY_norm_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDY_norm_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_norm_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ v_in_ptr,
                                           const BoutReal *__restrict__ f_in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexFDDY_norm_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_norm_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_norm_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDY_norm_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(f_in.getLocation());
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexDDX_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexDDX_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_off_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexDDX_on_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_on_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexDDX_off_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_off_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDX_off_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexDDY_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexDDY_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_off_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexDDY_on_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_on_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexDDY_off_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_off_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDY_off_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_on_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_off_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 2) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_off_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_on_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_on_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_on_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDZ_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexDDZ_off_DIFF_C4(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDZ_off_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDZ_off_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexDDZ_off_DIFF_C4_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexDDX_on_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 2) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_on_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexDDX_off_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 2) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_off_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_off_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_on_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexDDX_on_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_on_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDX_off_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexDDX_off_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDX_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDX_off_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDX_off_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexDDY_on_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 2) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_on_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexDDY_off_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 2) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_off_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_off_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_on_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                        const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexDDY_on_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_on_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexDDY_off_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ in_ptr, Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexDDY_off_DIFF_C4(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexDDY_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexDDY_off_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexDDY_off_DIFF_C4_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexD2DX2_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_on_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DX2_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexD2DX2_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DX2_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexD2DY2_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_on_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DY2_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexD2DY2_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DY2_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DZ2_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexD2DZ2_on_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DZ2_on_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DZ2_on_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DZ2_on_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DZ2_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexD2DZ2_off_DIFF_C2(const Field3D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DZ2_off_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 4) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DZ2_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0, 0);
  indexD2DZ2_off_DIFF_C2_field3d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexD2DX2_on_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_on_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DX2_on_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DX2_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexD2DX2_off_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DX2_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DX2_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DX2_off_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexD2DY2_on_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_on_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DY2_on_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexD2DY2_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                           const BoutReal *__restrict__ in_ptr,
                                           Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexD2DY2_off_DIFF_C2(const Field2D &in) {
  Mesh *msh = in.getMesh();
  output_debug.write("Using method indexD2DY2_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 4) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexD2DY2_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(in);
  const BoutReal *__restrict__ in_ptr = &in(0, 0);
  indexD2DY2_off_DIFF_C2_field2d(result_ptr, in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_on_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_U2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_on_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_off_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_off_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_on_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_on_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexVDDX_off_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_off_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_on_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_on_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field3D indexVDDX_off_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDX_off_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_on_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_U2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_on_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_off_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_off_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_on_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_on_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexVDDY_off_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_off_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_on_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_on_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field3D indexVDDY_off_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDY_off_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_on_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_off_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_on_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_on_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_on_DIFF_U2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_on_DIFF_U2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_on_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_off_DIFF_U2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_off_DIFF_U2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_off_DIFF_U2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_off_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_off_DIFF_U2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_on_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_on_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_on_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_on_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_off_DIFF_C2_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_off_DIFF_C2(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_off_DIFF_C2!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_off_DIFF_C2_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_on_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_on_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_on_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_on_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDZ_off_DIFF_C4_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexVDDZ_off_DIFF_C4(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDZ_off_DIFF_C4!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 5) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDZ_off_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexVDDZ_off_DIFF_C4_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_on_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_on_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_off_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_off_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_on_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_U2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_on_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_off_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_off_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_on_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_on_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexVDDX_off_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_off_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_on_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_on_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_on_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDX_off_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 2) {
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
static Field2D indexVDDX_off_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDX_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNx < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDX_off_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDX_off_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_on_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_on_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_off_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_off_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_on_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_U2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_on_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_U2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_off_DIFF_U2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_U2!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_U2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_off_DIFF_U2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_on_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_C2 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_on_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_C2_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexVDDY_off_DIFF_C2(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_C2!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_C2 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_off_DIFF_C2_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_on_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_on_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_on_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_on_DIFF_C4 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_on_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexVDDY_off_DIFF_C4_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 2) {
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
static Field2D indexVDDY_off_DIFF_C4(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexVDDY_off_DIFF_C4!\n");
#if CHECK > 0
  if (msh->LocalNy < 5) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexVDDY_off_DIFF_C4 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexVDDY_off_DIFF_C4_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexFDDX_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDX_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field3D indexFDDX_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDX_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexFDDY_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDY_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field3D indexFDDY_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDY_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDZ_on_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexFDDZ_on_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDZ_on_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDZ_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDZ_on_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_ZLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDZ_off_DIFF_U1_field3d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
  const int Nz = msh->LocalNz;
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
static Field3D indexFDDZ_off_DIFF_U1(const Field3D &v_in, const Field3D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDZ_off_DIFF_U1!\n");
  if (msh->LocalNz == 1) {
    return Field3D(0., msh);
  }
#if CHECK > 0
  if (msh->LocalNz < 3) {
    if (msh->xstart == 0) {
      // return Field3D(0.,msh);
      Field3D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDZ_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field3D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0, 0);
  indexFDDZ_off_DIFF_U1_field3d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_on_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexFDDX_on_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDX_on_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_XLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDX_off_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->xstart < 1) {
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
static Field2D indexFDDX_off_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDX_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNx < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDX_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDX_off_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_on_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                         const BoutReal *__restrict__ v_in_ptr,
                                         const BoutReal *__restrict__ f_in_ptr,
                                         Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexFDDY_on_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_on_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_on_DIFF_U1 - Not enough guards cells to "
                          "take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDY_on_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_YLOW);
  checkData(result);
  return result;
}

// This file is auto-generated - do not edit!
static void indexFDDY_off_DIFF_U1_field2d(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ v_in_ptr,
                                          const BoutReal *__restrict__ f_in_ptr,
                                          Mesh *msh) {
  const int Nx = msh->LocalNx;
  const int Ny = msh->LocalNy;
#if CHECK > 0
  if (msh->ystart < 1) {
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
static Field2D indexFDDY_off_DIFF_U1(const Field2D &v_in, const Field2D &f_in) {
  Mesh *msh = v_in.getMesh();
  ASSERT1(msh == f_in.getMesh());
  output_debug.write("Using method indexFDDY_off_DIFF_U1!\n");
#if CHECK > 0
  if (msh->LocalNy < 3) {
    if (msh->xstart == 0) {
      // return Field2D(0.,msh);
      Field2D result{msh};
      result = 0;
      return result;
    } else {
      throw BoutException("AiolosMesh::indexFDDY_off_DIFF_U1 - Not enough guards cells "
                          "to take derivative!");
    }
  }
#endif
  Field2D result(msh);
  result.allocate();
  BoutReal *__restrict__ result_ptr = &result(0, 0);
  checkData(v_in);
  const BoutReal *__restrict__ v_in_ptr = &v_in(0, 0);
  checkData(f_in);
  const BoutReal *__restrict__ f_in_ptr = &f_in(0, 0);
  indexFDDY_off_DIFF_U1_field2d(result_ptr, v_in_ptr, f_in_ptr, msh);
  result.setLocation(CELL_CENTRE);
  checkData(result);
  return result;
}
