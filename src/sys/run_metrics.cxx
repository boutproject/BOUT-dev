/*!************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 *********************************************************/

#include <bout/run_metrics.hxx>
#include <datafile.hxx>
#include <output.hxx>

/*!
 * Adds variables to the output file, for post-processing
 */
void RunMetrics::outputVars(Datafile &file) {
  file.add(simtime, "t_array", true);
  file.add(t_elapsed, "wall_time", true);
  file.add(wtime, "wtime", true);
  file.add(ncalls, "ncalls", true);
  file.add(ncalls_e, "ncalls_e", true);
  file.add(ncalls_i, "ncalls_i", true);
  file.add(wtime_rhs, "wtime_rhs", true);
  file.add(wtime_invert, "wtime_invert", true);
  file.add(wtime_comms, "wtime_comms", true);
  file.add(wtime_io, "wtime_io", true);
}

void RunMetrics::writeProgress(bool output_split) {

  if (!output_split) {
    output_progress.write("%.3e      %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", 
               simtime, ncalls, wtime,
               100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
               100.*wtime_invert/wtime,  // Inversions
               100.0*wtime_comms/wtime,  // Communications
               100.* wtime_io / wtime,      // I/O
               100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else

  } else {
    output_progress.write("%.3e      %5d            %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n",
               simtime, ncalls_e, ncalls_i, wtime,
               100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
               100.*wtime_invert/wtime,  // Inversions
               100.0*wtime_comms/wtime,  // Communications
               100.* wtime_io / wtime,      // I/O
               100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else
  }
}
