#ifndef __MONITOR_H__
#define __MONITOR_H__

#include "bout_types.hxx"
#include "bout/assert.hxx"
#include "utils.hxx"

#include <cmath>

class Options;
class Solver;

/// Return true if either \p a is a multiple of \p b or vice-versa
///
/// Assumes both arguments are greater than zero
inline bool isMultiple(BoutReal a, BoutReal b) {
  ASSERT2(a > 0);
  ASSERT2(b > 0);

  auto min = a>b?b:a;
  auto max = a>b?a:b;
  auto ratio = std::round(max/min);
  auto error = ratio*min - max;
  return (std::abs(error/max) < 1e-12);
}

/// Monitor baseclass for the Solver
///
/// Can be called ether with a specified frequency, or with the
/// frequency of the BOUT++ output monitor.
class Monitor{
  friend class Solver; ///< needs access to timestep and freq
public:
  /// A \p timestep_ of -1 defaults to the the frequency of the BOUT++
  /// output monitor
  Monitor(BoutReal timestep_ = -1) : timestep(timestep_) {};

  virtual ~Monitor() = default;

  /// Callback function for the solver, called after timestep_ has passed
  ///
  /// @param[in] solver The solver calling this monitor
  /// @param[in] time   The current simulation time
  /// @param[in] iter   The current simulation iteration
  /// @param[in] nout   The total number of iterations for this simulation
  ///
  /// @returns non-zero if simulation should be stopped
  virtual int call(Solver* solver, BoutReal time, int iter, int nout) = 0;

  /// Callback function for when a clean shutdown is initiated
  virtual void cleanup(){};

  virtual void outputVars(MAYBE_UNUSED(Options& options),
                          MAYBE_UNUSED(const std::string& time_dimension)) {}

protected:
  /// Get the currently set timestep for this monitor
  BoutReal getTimestep() const { return timestep; }

  /// Set the timestep for this Monitor
  ///
  /// Can only be called before the Monitor is added to a Solver
  void setTimestep(BoutReal new_timestep) {
    if (is_added) {
      throw BoutException("Monitor::set_timestep - Error: Monitor has already"
                          "been added to a Solver, so timestep cannot be changed.");
    }
    timestep = new_timestep;
  }

private:
  /// Set to true when Monitor is added to a Solver
  bool is_added{false};
  /// The desired physical timestep
  BoutReal timestep{-1};
  /// How often this monitor should be called, in internal Solver steps
  int period{1};
};

struct RunMetrics {
  public:
  /// cumulative wall clock time in seconds
  BoutReal t_elapsed = 0;
  /// time step's wall clock time in seconds
  BoutReal wtime = 0;

  /// number of RHS calls
  int ncalls = 0;
  /// number of RHS calls for fast timescale
  int ncalls_e = 0;
  /// number of RHS calls for slow timescale
  int ncalls_i = 0;

  /// wall time spent calculating RHS
  BoutReal wtime_rhs = 0;
  /// wall time spent inverting Laplacian
  BoutReal wtime_invert = 0;
  /// wall time spent communicating (part of RHS)
  BoutReal wtime_comms = 0;
  /// wall time spent on I/O
  BoutReal wtime_io = 0;

  // Derived metrics

  /// wall time per RHS evaluation
  BoutReal wtime_per_rhs = 0;
  /// wall time per fast timescale RHS evaluation
  BoutReal wtime_per_rhs_e = 0;
  /// wall time per slow timescale RHS evaluation
  BoutReal wtime_per_rhs_i = 0;

  /*!
   * Adds variables to the output file, for post-processing
   */
  void outputVars(Options& output_options);

  /*!
   * Calculates derived metrics
   */
  void calculateDerivedMetrics();

  /*!
   * Write job progress to screen
   */
  void writeProgress(BoutReal simtime, bool output_split);

};


#endif // __MONITOR_H__
