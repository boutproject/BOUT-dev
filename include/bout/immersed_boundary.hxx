#pragma once
#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
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

  //Mask function that is 1 in the plasma, 0 in the wall
  Field3D bndry_mask;
  //Mask function with ids into ghost cell data arrays.
  Field3D ghost_ids;
  //Coordinate fields.
  Field3D R3;
  Field3D Z3;
  //Partial face info.
  Field3D x_face_frac;
  Field3D z_face_frac;
  Field3D fx_grad_offset;
  Field3D fz_grad_offset;

  //Ghost cell data arrays.
  int num_ghosts = 0;
  //IB_TODO: Allow for reading arrays of ints not just BoutReals from mesh. Then dont need to lround data to nearest int w/ floating point error in indexing.
  Matrix<BoutReal> image_inds;
  Matrix<BoutReal> ghost_points;
  Matrix<BoutReal> image_points;
  Matrix<BoutReal> bndry_points;
  Matrix<BoutReal> normals;
  Array<BoutReal> norm_dist;

  //Image cell weights/ghost flag arrays.
  int num_weights = 0;
  Matrix<BoutReal> weights_dir;
  Matrix<BoutReal> weights_neu;
  Matrix<BoutReal> is_plasma;
  Matrix<BoutReal> ghost_use_bc;
  Array<BoutReal>  ghost_use_gs;

  //Boundary information.
  Field3D bound_ids;
  Field3D bound_counts;
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

  //IB_TODO: Useful function for converting BoutReal types to int/bool.
  //Needed because can only read Array/Matrix<BoutReal> from grid file now.
  //In python, the ints/bools are converted to float but represented exactly.
  template <typename T>
  inline T get_as(BoutReal x, BoutReal tol = 1e-12) const {
    static_assert(std::is_integral_v<T>,
      "Currently, get_as<T> only supports integral and bool types.");

    const long r = std::lround(x);

    // Catch values that are not really integer-like, e.g. 2.3, 0.49, etc.
    if (std::abs(x - static_cast<BoutReal>(r)) > tol) {
      throw BoutException("IB value not integer like...");
    }

    return static_cast<T>(r);
  }
  //Use this to get gridpoint inds associated with general point.
  inline int GetGridInd(BoutReal f) const {
    return static_cast<int>(f);
  }
};


inline std::unique_ptr<ImmersedBoundary> immBndry; //C++17 global.

#endif // IMMERSED_BOUNDARY_H
