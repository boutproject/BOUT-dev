/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef BOUT_PARALLELTRANSFORM_H
#define BOUT_PARALLELTRANSFORM_H

#include "bout/bout_types.hxx"
#include "bout/field3d.hxx"
#include "bout/options.hxx"
#include "bout/unused.hxx"

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
  ParallelTransform(Mesh& mesh_in, Options* opt = nullptr)
      : mesh(mesh_in),
        options(opt == nullptr ? Options::root()["mesh:paralleltransform"] : *opt) {}
  virtual ~ParallelTransform() = default;

  /// Given a 3D field, calculate and set the Y up down fields
  virtual void calcParallelSlices(Field3D& f) = 0;

  /// Calculate Yup and Ydown fields by integrating over mapped points
  /// This should be used for parallel divergence operators
  virtual void integrateParallelSlices(Field3D& f) { return calcParallelSlices(f); }

  /// Convert a field into field-aligned coordinates
  /// so that the y index is along the magnetic field
  virtual Field3D toFieldAligned(const Field3D& f,
                                 const std::string& region = "RGN_ALL") = 0;
  virtual FieldPerp toFieldAligned(const FieldPerp& f,
                                   const std::string& region = "RGN_ALL") = 0;

  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual Field3D fromFieldAligned(const Field3D& f,
                                   const std::string& region = "RGN_ALL") = 0;
  virtual FieldPerp fromFieldAligned(const FieldPerp& f,
                                     const std::string& region = "RGN_ALL") = 0;

  /// Field2D are axisymmetric, so transformation to or from field-aligned coordinates is
  /// a null operation.
  virtual Field2D toFieldAligned(const Field2D& f,
                                 const std::string& UNUSED(region) = "RGN_ALL") {
    return f;
  }
  virtual Field2D fromFieldAligned(const Field2D& f,
                                   const std::string& UNUSED(region) = "RGN_ALL") {
    return f;
  }

  virtual bool canToFromFieldAligned() const = 0;

  struct PositionsAndWeights {
    int i, j, k;
    BoutReal weight;
  };

  virtual std::vector<PositionsAndWeights> getWeightsForYUpApproximation(int i, int j,
                                                                         int k) {
    return getWeightsForYApproximation(i, j, k, 1);
  }
  virtual std::vector<PositionsAndWeights> getWeightsForYDownApproximation(int i, int j,
                                                                           int k) {
    return getWeightsForYApproximation(i, j, k, -1);
  }
  virtual std::vector<PositionsAndWeights>
  getWeightsForYApproximation(int UNUSED(i), int UNUSED(j), int UNUSED(k),
                              int UNUSED(yoffset)) {
    throw BoutException("ParallelTransform::getWeightsForYApproximation not implemented "
                        "in this subclass");
  }

  /// Output variables used by a ParallelTransform instance to \p output_options
  virtual void outputVars([[maybe_unused]] Options& output_options) {}

  /// If \p twist_shift_enabled is true, does a `Field3D` with Y direction \p ytype
  /// require a twist-shift at branch cuts on closed field lines?
  virtual bool requiresTwistShift(bool twist_shift_enabled, YDirectionType ytype) = 0;

protected:
  /// This method should be called in the constructor to check that if the grid
  /// has a 'parallel_transform' variable, it has the correct value
  virtual void checkInputGrid() = 0;

  Mesh& mesh;       ///< The mesh this paralleltransform is part of
  Options& options; ///< Options for this ParallelTransform
};

/*!
 * This class implements the simplest form of ParallelTransform
 * where the domain is a logically rectangular domain, and
 * yup() and ydown() refer to the same field.
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  ParallelTransformIdentity(Mesh& mesh_in, Options* opt = nullptr)
      : ParallelTransform(mesh_in, opt) {

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
  Field3D toFieldAligned(const Field3D& f,
                         const std::string& UNUSED(region) = "RGN_ALL") override {
    ASSERT2(f.getDirectionY() == YDirectionType::Standard);
    Field3D result = f;
    return result.setDirectionY(YDirectionType::Aligned);
  }
  FieldPerp toFieldAligned(const FieldPerp& f,
                           const std::string& UNUSED(region) = "RGN_ALL") override {
    ASSERT2(f.getDirectionY() == YDirectionType::Standard);
    FieldPerp result = f;
    return result.setDirectionY(YDirectionType::Aligned);
  }

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  Field3D fromFieldAligned(const Field3D& f,
                           const std::string& UNUSED(region) = "RGN_ALL") override {
    ASSERT2(f.getDirectionY() == YDirectionType::Aligned);
    Field3D result = f;
    return result.setDirectionY(YDirectionType::Standard);
  }
  FieldPerp fromFieldAligned(const FieldPerp& f,
                             const std::string& UNUSED(region) = "RGN_ALL") override {
    ASSERT2(f.getDirectionY() == YDirectionType::Aligned);
    FieldPerp result = f;
    return result.setDirectionY(YDirectionType::Standard);
  }

  virtual std::vector<PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int yoffset) override {
    return {{i, j + yoffset, k, 1.0}};
  }

  bool canToFromFieldAligned() const override { return true; }

  bool requiresTwistShift(bool twist_shift_enabled,
                          YDirectionType UNUSED(ytype)) override {
    // All Field3Ds require twist-shift, because all are effectively field-aligned, but
    // allow twist-shift to be turned off by twist_shift_enabled
    return twist_shift_enabled;
  }

protected:
  void checkInputGrid() override;
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
  ShiftedMetric(Mesh& mesh, CELL_LOC location, Field2D zShift, BoutReal zlength_in,
                Options* opt = nullptr);

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
  Field3D toFieldAligned(const Field3D& f,
                         const std::string& region = "RGN_ALL") override;
  FieldPerp toFieldAligned(const FieldPerp& f,
                           const std::string& region = "RGN_ALL") override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  Field3D fromFieldAligned(const Field3D& f,
                           const std::string& region = "RGN_ALL") override;
  FieldPerp fromFieldAligned(const FieldPerp& f,
                             const std::string& region = "RGN_ALL") override;

  std::vector<PositionsAndWeights>
  getWeightsForYApproximation(int UNUSED(i), int UNUSED(j), int UNUSED(k),
                              int UNUSED(yoffset)) override {
    throw BoutException("ParallelTransform::getWeightsForYApproximation not implemented"
                        "for `type = shifted`. Try `type = shiftedinterp`");
  }

  bool canToFromFieldAligned() const override { return true; }

  /// Save zShift to the output
  void outputVars(Options& output_options) override;

  bool requiresTwistShift(bool twist_shift_enabled, YDirectionType ytype) override {
    // Twist-shift only if field-aligned
    if (ytype == YDirectionType::Aligned and not twist_shift_enabled) {
      throw BoutException("'twistshift = true' is required to communicate field-aligned "
                          "Field3Ds when using ShiftedMetric.");
    }
    return ytype == YDirectionType::Aligned;
  }

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

  /// Helper POD for parallel slice phase shifts
  struct ParallelSlicePhase {
    Tensor<dcomplex> phase_shift;
    int y_offset;
  };

  /// Cache of phase shifts for the parallel slices. Slices are stored
  /// in the following order:
  ///     {+1, ..., +n, -1, ..., -n}
  /// slice[i] stores offset i+1
  /// slice[n + i] stores offset -(i+1)
  /// where i goes from 0 to (n-1), with n the number of y guard cells
  std::vector<ParallelSlicePhase> parallel_slice_phases;

  /*!
   * Shift a 2D field in Z.
   * Since 2D fields are constant in Z, this has no effect
   */
  Field2D shiftZ(const Field2D& f, const Field2D& UNUSED(zangle),
                 const std::string UNUSED(region) = "RGN_NOX") const {
    return f;
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
  Field3D shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                 const YDirectionType y_direction_out,
                 const std::string& region = "RGN_NOX") const;
  FieldPerp shiftZ(const FieldPerp& f, const Tensor<dcomplex>& phs,
                   const YDirectionType y_direction_out,
                   const std::string& region = "RGN_NOX") const;

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

#endif // BOUT_PARALLELTRANSFORM_H
