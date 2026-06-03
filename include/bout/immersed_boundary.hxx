#pragma once
#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <bout/array.hxx>
#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/region.hxx>
#include <bout/utils.hxx>

class Mesh;

/// Class for handling immersed boundary methods, meant for the
/// device wall.
class ImmersedBoundary {
public:
  ImmersedBoundary();
  void SetBoundary(Field3D& f);
  void ComputeBoundaryFluxes(const Field3D& a, const Field3D& f, Field3D& result) const;

  bool IsInside(const Ind3D& ind) const;
  bool IsGhost(const Ind3D& ind) const;
  bool IsCutCell(const Ind3D& ind) const;
  float xFaceFrac(const Ind3D& ind) const;
  float zFaceFrac(const Ind3D& ind) const;
  int xFaceGradOffset(const Ind3D& ind) const;
  int zFaceGradOffset(const Ind3D& ind) const;

  void FieldSetup(Field3D& f);
  void FloorField(Field3D& f, const float val = 0.0) const;
  bool CheckFieldSetUp(const std::string& name);

private:
  enum class BoundCond {DIRICHLET, NEUMANN, SIZE};
  using BC_Info = std::pair<BoundCond, BoutReal>;

  Mesh* localmesh_;

  /// Mask function that is 1 in the plasma, 0 in the wall
  Field3D bndry_mask;
  /// Mask function with ids into ghost cell data arrays.
  Field3D ghost_ids;
  /// Coordinate fields.
  Field3D R3;
  Field3D Z3;
  //Partial face info.
  Field3D x_face_frac;
  Field3D z_face_frac;
  Field3D fx_grad_offset;
  Field3D fz_grad_offset;

  /// Ghost cell data arrays.
  int num_ghosts = 0;
  //IB_TODO: Allow for reading arrays of ints not just BoutReals from mesh. Then dont need to lround data to nearest int w/ floating point error in indexing.
  Matrix<BoutReal> image_inds;
  Matrix<BoutReal> ghost_points;
  Matrix<BoutReal> image_points;
  Matrix<BoutReal> bndry_points;
  Matrix<BoutReal> normals;
  Array<BoutReal> norm_dist;

  /// Image cell weights/ghost flag arrays.
  int num_weights = 0;
  Matrix<BoutReal> weights;
  Matrix<BoutReal> is_plasma;
  Array<int> all_plasma;
  Array<int> ghost_count;

  //Boundary information.
  Field3D bound_ids;
  int num_bounds = 0;
  BoutReal s1 = 0;
  BoutReal s2 = 0;
  Array<BoutReal> bds;
  Array<BoutReal> bd_area;
  Matrix<BoutReal> mid_pts;
  Matrix<BoutReal> bnorms;
  Matrix<BoutReal> bbase_inds;
  Matrix<BoutReal> bweights;
  std::unordered_map<std::string, BC_Info> bc_map;

  BoutReal GetGhostValue(const BoutReal image_val, const int gid,
                    const BoutReal bc, const BoundCond bc_type) const;
  BoutReal GetImageValue(Field3D& f, const int gid, const BoutReal bc_val,
                    const BoundCond bc_type) const;

  const std::string bc_exception = "Invalid boundary condition specified for immersed boundary.";
  const std::string bc_key = "bndry_wall";
  std::pair<BoundCond, BoutReal> ReadBC(const std::string& bc_info) const;

  bool IsBadInterpForMPI(const int global_indx) const;
  void CheckInterpOkWithMPI(const int global_indx, const int cell_id, const int proc_idx,
                      const std::string& description) const;

  // Solve 4x4 A x = b (A is copied by value; partial pivoting)
  // IB_TODO: Clean up/double check this functionality.
  template <typename T>
  static std::array<T,4> solve4x4(std::array<std::array<T,4>,4> A,
                                  const std::array<T,4> b) {
    int p[4] = {0,1,2,3};
    // LU with partial pivoting (Doolittle)
    for (int k = 0; k < 4; ++k) {
      // pivot
      int imax = k;
      T amax = std::abs(A[p[k]][k]);
      for (int i = k+1; i < 4; ++i) {
        T v = std::abs(A[p[i]][k]);
        if (v > amax) { amax = v; imax = i; }
      }
      std::swap(p[k], p[imax]);
      T pivot = A[p[k]][k];
      if (std::abs(pivot) == T(0)) { throw BoutException("Singular 4x4 Vandermonde"); }
      // eliminate
      for (int i = k+1; i < 4; ++i) {
        T m = A[p[i]][k] / pivot;
        A[p[i]][k] = m;
        for (int j = k+1; j < 4; ++j)
          A[p[i]][j] -= m * A[p[k]][j];
      }
    }
    // forward subst: Ly = b(p)
    T y[4];
    for (int i = 0; i < 4; ++i) {
      T sum = b[p[i]];
      for (int j = 0; j < i; ++j) sum -= A[p[i]][j] * y[j];
      y[i] = sum;
    }
    // back subst: Ux = y
    std::array<T,4> x;
    for (int i = 3; i >= 0; --i) {
      T sum = y[i];
      for (int j = i+1; j < 4; ++j) sum -= A[p[i]][j] * x[j];
      x[i] = sum / A[p[i]][i];
    }
    return x;
  }
};

inline std::unique_ptr<ImmersedBoundary> immBndry; //C++17 global.

#endif // IMMERSED_BOUNDARY_H
