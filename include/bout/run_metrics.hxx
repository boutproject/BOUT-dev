/**************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Run Metrics class by J. T. Parker 2019
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
 *
 **************************************************************************/

/// Run Metrics class
///
/// The `Run Metrics` class holds non-physics information about the job execution,
/// such as wall clock times and number of iterations.

#ifndef __RUNMETRICS_H__
#define __RUNMETRICS_H__


class Datafile;
#include "bout_types.hxx"


class RunMetrics {
  public:
  /// cumulative wall clock time in seconds
  BoutReal t_elapsed;
  /// time step's wall clock time in seconds
  BoutReal wtime;

  /// number of RHS calls
  int ncalls;
  /// number of RHS calls for fast timescale
  int ncalls_e;
  /// number of RHS calls for slow timescale
  int ncalls_i;

  /// wall time spent calculating RHS
  BoutReal wtime_rhs;
  /// wall time spent inverting Laplacian
  BoutReal wtime_invert;
  /// wall time spent communicating (part of RHS)
  BoutReal wtime_comms;
  /// wall time spent on I/O
  BoutReal wtime_io;

  // Derived metrics

  /// wall time per RHS evaluation
  BoutReal wtime_per_rhs;
  /// wall time per fast timescale RHS evaluation
  BoutReal wtime_per_rhs_e;
  /// wall time per slow timescale RHS evaluation
  BoutReal wtime_per_rhs_i;

  /*!
   * Adds variables to the output file, for post-processing
   */
  void outputVars(Datafile &file);

  /*!
   * Calculates derived metrics
   */
  void calculateDerivedMetrics();

  /*!
   * Write job progress to screen 
   */
  void writeProgress(BoutReal simtime, bool output_split);

};

#endif // __RUNMETRICS_H__
