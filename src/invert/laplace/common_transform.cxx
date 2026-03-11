#include "common_transform.hxx"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/constants.hxx"
#include "bout/dcomplex.hxx"
#include "bout/fft.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/mesh.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/utils.hxx"

#include <iterator>

FFTTransform::FFTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye,
                           int zs, int ze, int inbndry, int outbndry,
                           bool inner_boundary_set_on_first_x,
                           bool outer_boundary_set_on_last_x, bool zero_DC)
    : zlength(getUniform(mesh.getCoordinatesConst()->zlength())), nmode(nmode), xs(xs),
      xe(xe), ys(ys), ye(ye), zs(zs), ze(ze), nx(xe - xs + 1), ny(ye - ys + 1),
      nz(ze - zs + 1), nxny(nx * ny), nsys(nmode * ny), inbndry(inbndry),
      outbndry(outbndry), inner_boundary_set_on_first_x(inner_boundary_set_on_first_x),
      outer_boundary_set_on_last_x(outer_boundary_set_on_last_x), zero_DC(zero_DC),
      localmesh(&mesh) {}

auto FFTTransform::forward(const Laplacian& laplacian, const Field3D& rhs,
                           const Field3D& x0, const Field2D& Acoef, const Field2D& C1coef,
                           const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nsys, nx);

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X and Y indices, including boundaries but not guard cells
    // (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      // Take FFT in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        rfft(&(x0(ix, iy, zs)), nz, std::begin(k1d));
      } else {
        rfft(&(rhs(ix, iy, zs)), nz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(((iy - ys) * nmode) + kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nsys; ind++) {
      const int iy = ys + (ind / nmode);
      const int kz = ind % nmode;

      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(ind, 0), &result.b(ind, 0), &result.c(ind, 0),
                             &result.bcmplx(ind, 0), iy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto FFTTransform::backward(const Field3D& rhs,
                            const Matrix<dcomplex>& xcmplx3D) const -> Field3D {
  Field3D x{emptyFrom(rhs)};

  // FFT back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx3D(((iy - ys) * nmode) + kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      irfft(std::begin(k1d), nz, &(x(ix, iy, zs)));
    }
  }

  return x;
}

auto FFTTransform::forward(const Laplacian& laplacian, const FieldPerp& rhs,
                           const FieldPerp& x0, const Field2D& Acoef,
                           const Field2D& C1coef, const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nmode, nx);
  const int jy = rhs.getIndex();

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X indices, including boundaries but not guard
    // cells (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ix = xs; ix <= xe; ix++) {
      // Take FFT in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        rfft(&(x0(ix, zs)), nz, std::begin(k1d));
      } else {
        rfft(&(rhs(ix, zs)), nz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int kz = 0; kz < nmode; kz++) {
      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(kz, 0), &result.b(kz, 0), &result.c(kz, 0),
                             &result.bcmplx(kz, 0), jy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto FFTTransform::backward(const FieldPerp& rhs,
                            const Matrix<dcomplex>& xcmplx) const -> FieldPerp {
  FieldPerp x{emptyFrom(rhs)};

  // FFT back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ix = xs; ix <= xe; ++ix) {
      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx(kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      irfft(std::begin(k1d), nz, &(x(ix, zs)));
    }
  }

  checkData(x);

  return x;
}

DSTTransform::DSTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye,
                           int zs, int ze, int inbndry, int outbndry,
                           bool inner_boundary_set_on_first_x,
                           bool outer_boundary_set_on_last_x, bool zero_DC)
    : zlength(getUniform(mesh.getCoordinatesConst()->zlength())), nmode(nmode), xs(xs),
      xe(xe), ys(ys), ye(ye), zs(zs), ze(ze), nx(xe - xs + 1), ny(ye - ys + 1),
      nz(ze - zs + 1), nxny(nx * ny), nsys(nmode * ny), inbndry(inbndry),
      outbndry(outbndry), inner_boundary_set_on_first_x(inner_boundary_set_on_first_x),
      outer_boundary_set_on_last_x(outer_boundary_set_on_last_x), zero_DC(zero_DC),
      localmesh(&mesh) {}

auto DSTTransform::forward(const Laplacian& laplacian, const Field3D& rhs,
                           const Field3D& x0, const Field2D& Acoef, const Field2D& C1coef,
                           const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nsys, nx);

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>(nz);

    // Loop over X and Y indices, including boundaries but not guard cells
    // (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      // Take DST in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        bout::fft::DST(&(x0(ix, iy, zs + 1)), nz - 2, std::begin(k1d));
      } else {
        bout::fft::DST(&(rhs(ix, iy, zs + 1)), nz - 2, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(((iy - ys) * nmode) + kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nsys; ind++) {
      const int iy = ys + (ind / nmode);
      const int kz = ind % nmode;

      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(ind, 0), &result.b(ind, 0), &result.c(ind, 0),
                             &result.bcmplx(ind, 0), iy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto DSTTransform::backward(const Field3D& rhs,
                            const Matrix<dcomplex>& xcmplx3D) const -> Field3D {
  Field3D x{emptyFrom(rhs)};

  // DST back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx3D(((iy - ys) * nmode) + kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      bout::fft::DST_rev(std::begin(k1d), nz - 2, &(x(ix, iy, zs + 1)));

      x(ix, iy, zs) = -x(ix, iy, zs + 2);
      x(ix, iy, ze) = -x(ix, iy, ze - 2);
    }
  }

  return x;
}

auto DSTTransform::forward(const Laplacian& laplacian, const FieldPerp& rhs,
                           const FieldPerp& x0, const Field2D& Acoef,
                           const Field2D& C1coef, const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nmode, nx);
  const int jy = rhs.getIndex();

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X indices, including boundaries but not guard
    // cells (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ix = xs; ix <= xe; ix++) {
      // Take DST in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        rfft(&(x0(ix, zs)), nz, std::begin(k1d));
      } else {
        rfft(&(rhs(ix, zs)), nz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int kz = 0; kz < nmode; kz++) {
      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(kz, 0), &result.b(kz, 0), &result.c(kz, 0),
                             &result.bcmplx(kz, 0), jy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto DSTTransform::backward(const FieldPerp& rhs,
                            const Matrix<dcomplex>& xcmplx) const -> FieldPerp {
  FieldPerp x{emptyFrom(rhs)};

  // DST back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ix = xs; ix < xe; ++ix) {
      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx(kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      irfft(std::begin(k1d), nz, &(x(ix, zs)));
    }
  }

  checkData(x);

  return x;
}
