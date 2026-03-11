#pragma once

#include "bout/bout_types.hxx"
#include "bout/dcomplex.hxx"
#include "bout/utils.hxx"

class Laplacian;
class Field2D;
class Field3D;
class FieldPerp;
class Mesh;

/// Helper class to do FFTs for `LaplaceCyclic`, `LaplacePCR`, `LaplacePCR_THOMAS`
class FFTTransform {
  /// Physical length of the Z domain
  BoutReal zlength;
  /// Number of unfiltered Fourier modes
  int nmode;

  int xs;
  int xe;
  int ys;
  int ye;
  int zs;
  int ze;

  int nx;
  int ny;
  /// Number of (interior) points in Z
  int nz;
  int nxny;
  /// Number of systems to solve = number of unfiltered Fourier modes times number of y
  /// points
  int nsys;
  /// Number of inner boundary cells
  int inbndry;
  /// Number of outer boundary cells
  int outbndry;

  bool inner_boundary_set_on_first_x;
  bool outer_boundary_set_on_last_x;

  bool zero_DC;

  const Mesh* localmesh;

public:
  FFTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye, int zs,
               int ze, int inbndry, int outbndry, bool inner_boundary_set_on_first_x,
               bool outer_boundary_set_on_last_x, bool zero_DC);

  struct Matrices {
    Matrix<dcomplex> a;
    Matrix<dcomplex> b;
    Matrix<dcomplex> c;
    Matrix<dcomplex> bcmplx;

    Matrices(int nsys, int nx)
        : a(nsys, nx), b(nsys, nx), c(nsys, nx), bcmplx(nsys, nx) {}
  };

  auto forward(const Laplacian& laplacian, const Field3D& rhs, const Field3D& x0,
               const Field2D& Acoef, const Field2D& C1coef, const Field2D& C2coef,
               const Field2D& Dcoef) const -> Matrices;
  auto backward(const Field3D& rhs, const Matrix<dcomplex>& xcmplx3D) const -> Field3D;

  auto forward(const Laplacian& laplacian, const FieldPerp& rhs, const FieldPerp& x0,
               const Field2D& Acoef, const Field2D& C1coef, const Field2D& C2coef,
               const Field2D& Dcoef) const -> Matrices;
  auto backward(const FieldPerp& rhs, const Matrix<dcomplex>& xcmplx) const -> FieldPerp;
};

/// Helper class to do discrete sin transforms (DSTs) for `LaplaceCyclic`,
/// `LaplacePCR`, `LaplacePCR_THOMAS`
class DSTTransform {
  /// Physical length of the Z domain
  BoutReal zlength;
  /// Number of unfiltered Fourier modes
  int nmode;

  int xs;
  int xe;
  int ys;
  int ye;
  int zs;
  int ze;

  int nx;
  int ny;
  /// Number of (interior) points in Z
  int nz;
  int nxny;
  /// Number of systems to solve = number of unfiltered Fourier modes times number of y
  /// points
  int nsys;
  /// Number of inner boundary cells
  int inbndry;
  /// Number of outer boundary cells
  int outbndry;

  bool inner_boundary_set_on_first_x;
  bool outer_boundary_set_on_last_x;

  bool zero_DC;

  const Mesh* localmesh;

public:
  DSTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye, int zs,
               int ze, int inbndry, int outbndry, bool inner_boundary_set_on_first_x,
               bool outer_boundary_set_on_last_x, bool zero_DC);

  struct Matrices {
    Matrix<dcomplex> a;
    Matrix<dcomplex> b;
    Matrix<dcomplex> c;
    Matrix<dcomplex> bcmplx;

    Matrices(int nsys, int nx)
        : a(nsys, nx), b(nsys, nx), c(nsys, nx), bcmplx(nsys, nx) {}
  };

  auto forward(const Laplacian& laplacian, const Field3D& rhs, const Field3D& x0,
               const Field2D& Acoef, const Field2D& C1coef, const Field2D& C2coef,
               const Field2D& Dcoef) const -> Matrices;
  auto backward(const Field3D& rhs, const Matrix<dcomplex>& xcmplx3D) const -> Field3D;

  auto forward(const Laplacian& laplacian, const FieldPerp& rhs, const FieldPerp& x0,
               const Field2D& Acoef, const Field2D& C1coef, const Field2D& C2coef,
               const Field2D& Dcoef) const -> Matrices;
  auto backward(const FieldPerp& rhs, const Matrix<dcomplex>& xcmplx) const -> FieldPerp;
};
