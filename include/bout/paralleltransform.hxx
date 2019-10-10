/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include "bout_types.hxx"
#include "field3d.hxx"
#include "unused.hxx"

class Mesh;

/*!
 * Calculates the values of a field along the magnetic
 * field (y direction)
 *
 * This is used inside Mesh to represent a variety of possible
 * numerical schemes, including Shifted Metric and FCI
 *
 */
class ParallelTransform {
public:
  ParallelTransform(Mesh& mesh_in) : mesh(mesh_in) {}
  virtual ~ParallelTransform() = default;

  /// Given a 3D field, calculate and set the Y up down fields
  virtual void calcParallelSlices(Field3D &f) = 0;

  [[gnu::deprecated("Please use ParallelTransform::calcParallelSlices instead")]]
  void calcYupYdown(Field3D& f) {
    calcParallelSlices(f);
  }

  /// Calculate Yup and Ydown fields by integrating over mapped points
  /// This should be used for parallel divergence operators
  virtual void integrateParallelSlices(Field3D &f) {
    return calcParallelSlices(f);
  }

  [[gnu::deprecated("Please use ParallelTransform::integrateParallelSlices instead")]]
  void integrateYupYdown(Field3D& f) {
    integrateParallelSlices(f);
  }
  
  /// Convert a field into field-aligned coordinates
  /// so that the y index is along the magnetic field
  virtual const Field3D toFieldAligned(const Field3D &f, const std::string& region = "RGN_ALL") = 0;
  [[gnu::deprecated("Please use toFieldAligned(const Field3D& f, "
      "const std::string& region = \"RGN_ALL\") instead")]]
  const Field3D toFieldAligned(const Field3D &f, REGION region) {
    return toFieldAligned(f, toString(region));
  }
  virtual const FieldPerp toFieldAligned(const FieldPerp &f, const std::string& region = "RGN_ALL") = 0;
  [[gnu::deprecated("Please use toFieldAligned(const FieldPerp& f, "
      "const std::string& region = \"RGN_ALL\") instead")]]
  const FieldPerp toFieldAligned(const FieldPerp &f, REGION region) {
    return toFieldAligned(f, toString(region));
  }
  
  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual const Field3D fromFieldAligned(const Field3D &f, const std::string& region = "RGN_ALL") = 0;
  [[gnu::deprecated("Please use fromFieldAligned(const Field3D& f, "
      "const std::string& region = \"RGN_ALL\") instead")]]
  const Field3D fromFieldAligned(const Field3D &f, REGION region) {
    return fromFieldAligned(f, toString(region));
  }
  virtual const FieldPerp fromFieldAligned(const FieldPerp &f, const std::string& region = "RGN_ALL") = 0;
  [[gnu::deprecated("Please use fromFieldAligned(const FieldPerp& f, "
      "const std::string& region = \"RGN_ALL\") instead")]]
  const FieldPerp fromFieldAligned(const FieldPerp &f, REGION region) {
    return fromFieldAligned(f, toString(region));
  }

  virtual bool canToFromFieldAligned() = 0;

  /// Output variables used by a ParallelTransform instance to the dump files
  virtual void outputVars(Datafile& UNUSED(file)) {}

  /// If \p twist_shift_enabled is true, apply twist-shift to `Field3D f`
  virtual void applyTwistShift(Field3D& f, bool twist_shift_enabled) = 0;

  /// Shift part of a Field3D by a given angle in Z
  ///
  /// @param[inout] var  The variable to be modified in-place
  /// @param[in] jx      X index
  /// @param[in] jy      Y index
  /// @param[in] zangle  The Z angle to apply
  virtual void shiftZ(Field3D& f, int jx, int jy, BoutReal zangle);

  /// Apply a phase shift by a given angle \p zangle in Z to all points
  ///
  /// @param[inout] var  The variable to modify in-place
  /// @param[in] zangle  The angle to shift by in Z
  /// @param[in] rgn     The region to calculate the result over
  void shiftZ(Field3D &var, BoutReal zangle, const std::string& rgn="RGN_ALL");

protected:
  /// This method should be called in the constructor to check that if the grid
  /// has a 'parallel_transform' variable, it has the correct value
  virtual void checkInputGrid() = 0;

  Mesh &mesh; ///< The mesh this paralleltransform is part of
};

/*!
 * This class implements the simplest form of ParallelTransform
 * where the domain is a logically rectangular domain, and
 * yup() and ydown() refer to the same field.
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  ParallelTransformIdentity(Mesh& mesh_in, std::vector<BoutReal> ShiftAngle_in)
      : ParallelTransform(mesh_in), ShiftAngle(std::move(ShiftAngle_in)) {

    // check the coordinate system used for the grid data source
    ParallelTransformIdentity::checkInputGrid();
  }

  /*!
   * Merges the yup and ydown() fields of f, so that
   * f.yup() = f.ydown() = f
   */
  void calcParallelSlices(Field3D& f) override;

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D toFieldAligned(const Field3D& f, const std::string& UNUSED(region) = "RGN_ALL") override {
    Field3D result = f;
    return result.setDirectionY(YDirectionType::Aligned);
  }
  const FieldPerp toFieldAligned(const FieldPerp& f, const std::string& UNUSED(region) = "RGN_ALL") override {
    FieldPerp result = f;
    return result.setDirectionY(YDirectionType::Aligned);
  }

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D& f, const std::string& UNUSED(region) = "RGN_ALL") override {
    Field3D result = f;
    return result.setDirectionY(YDirectionType::Standard);
  }
  const FieldPerp fromFieldAligned(const FieldPerp& f, const std::string& UNUSED(region) = "RGN_ALL") override {
    FieldPerp result = f;
    return result.setDirectionY(YDirectionType::Standard);
  }

  bool canToFromFieldAligned() override { return true; }

  void applyTwistShift(Field3D& f, bool twist_shift_enabled);

protected:
  void checkInputGrid() override;

  std::vector<BoutReal> ShiftAngle;  ///< Angle for twist-shift location
};

/*!
 * Shifted metric method
 * Each Y location is shifted in Z with respect to its neighbours
 * so that the grid is orthogonal in X-Z, but requires interpolation
 * to calculate the values of points along field-lines.
 *
 * In this implementation the interpolation is done using FFTs in Z
 */
class ShiftedMetric : public ParallelTransform {
public:
  ShiftedMetric() = delete;
  ShiftedMetric(Mesh& mesh, CELL_LOC location, Field2D zShift, BoutReal zlength_in);

  /*!
   * Calculates the yup() and ydown() fields of f
   * by taking FFTs in Z and applying a phase shift.
   */
  void calcParallelSlices(Field3D& f) override;

  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually).
   *
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D& f, const std::string& region = "RGN_ALL") override;
  const FieldPerp toFieldAligned(const FieldPerp& f,
                                 const std::string& region = "RGN_ALL") override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D& f,
                                 const std::string& region = "RGN_ALL") override;
  const FieldPerp fromFieldAligned(const FieldPerp& f,
                                   const std::string& region = "RGN_ALL") override;

  bool canToFromFieldAligned() override { return true; }

  /// Save zShift to the output
  void outputVars(Datafile& file) override;

  void applyTwistShift(Field3D& f, bool twist_shift_enabled);

protected:
  void checkInputGrid() override;

private:
  CELL_LOC location{CELL_CENTRE};

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;

  /// Length of the z-domain in radians
  BoutReal zlength{0.};

  int nmodes;

  Tensor<dcomplex> toAlignedPhs;   ///< Cache of phase shifts for transforming from X-Z
                                   /// orthogonal coordinates to field-aligned coordinates
  Tensor<dcomplex> fromAlignedPhs; ///< Cache of phase shifts for transforming from
                                   /// field-aligned coordinates to X-Z orthogonal
                                   /// coordinates

  Matrix<dcomplex> twistShiftDownPhs; ///< Cache of phase shifts for twist-shift by
                                      /// ShiftAngle going down in Y
  Matrix<dcomplex> twistShiftUpPhs;   ///< Cache of phase shifts for twist-shift by
                                      /// ShiftAngle going up in Y

  /// Helper POD for parallel slice phase shifts
  struct ParallelSlicePhase {
    Tensor<dcomplex> phase_shift;
    int y_offset;
  };

  /// Cache of phase shifts for the parallel slices. Slices are stored
  /// in the following order:
  ///     {+1, ..., +n, -1, ..., -n}
  /// slice[i] stores offset i+1
  /// slice[2*i + 1] stores offset -(i+1)
  /// where i goes from 0 to (n-1), with n the number of y guard cells
  std::vector<ParallelSlicePhase> parallel_slice_phases;

  /*!
   * Shift a 2D field in Z.
   * Since 2D fields are constant in Z, this has no effect
   */
  const Field2D shiftZ(const Field2D& f, const Field2D& UNUSED(zangle),
                       const std::string UNUSED(region) = "RGN_NOX") const {
    return f;
  };
  [[gnu::deprecated("Please use shiftZ(const Field2D& f, const Field2D& zangle, "
      "const std::string& region = \"RGN_NOX\") instead")]]
  const Field2D shiftZ(const Field2D& f, const Field2D& UNUSED(zangle),
                       REGION UNUSED(region)) const {
    return f;
  };

  /*!
   * Shift a 3D field \p f in Z by the given \p zangle
   *
   * @param[in] f  The field to shift
   * @param[in] zangle   Toroidal angle (z)
   *
   */
  const Field3D shiftZ(const Field3D& f, const Field2D& zangle,
                       const std::string& region = "RGN_NOX") const;
  [[gnu::deprecated("Please use shiftZ(const Field3D& f, const Field2D& zangle, "
      "const std::string& region = \"RGN_NOX\") instead")]]
  const Field3D shiftZ(const Field3D& f, const Field2D& zangle,
                       REGION region) const {
    return shiftZ(f, zangle, toString(region));
  };

  /*!
   * Shift a 3D field or FieldPerp \p f by the given phase \p phs in Z
   *
   * Calculates FFT in Z, multiplies by the complex phase
   * and inverse FFTS.
   *
   * @param[in] f  The field to shift
   * @param[in] phs  The phase to shift by
   * @param[in] y_direction_out  The value to set yDirectionType of the result to
   */
  const Field3D shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                       const YDirectionType y_direction_out,
                       const std::string& region = "RGN_NOX") const;
  const FieldPerp shiftZ(const FieldPerp& f, const Tensor<dcomplex>& phs,
                         const YDirectionType y_direction_out,
                         const std::string& region = "RGN_NOX") const;

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length \p len
   * @param[in] len  Length of the in and out arrays
   * @param[in] zangle  The angle (z coordinate) to shift by
   * @param[out] out  A 1D array of length \p len, already allocated
   */
  void shiftZ(const BoutReal* in, int len, BoutReal zangle, BoutReal* out) const;

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length mesh.LocalNz
   * @param[in] phs Phase shift, assumed to have length (mesh.LocalNz/2 + 1) i.e. the
   * number of modes
   * @param[out] out  A 1D array of length mesh.LocalNz, already allocated
   */
  void shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out) const;

  /// Calculate and store the phases for to/from field aligned and for
  /// the parallel slices using zShift
  void cachePhases();

  /// Shift a 3D field \p f in Z to all the parallel slices in \p phases
  ///
  /// @param[in] f      The field to shift
  /// @param[in] phases The phase and offset information for each parallel slice
  /// @return The shifted parallel slices
  std::vector<Field3D> shiftZ(const Field3D& f,
                              const std::vector<ParallelSlicePhase>& phases) const;
};

#endif // __PARALLELTRANSFORM_H__
