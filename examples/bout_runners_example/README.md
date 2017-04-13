# bout_runners_example

A tutorial on how to use `bout_runners`.

Extra documentation of `bout_runners` can be found in the
docstring of `bout_runners`.

## Contents
### The program:

* `diffusion_3D.cxx` - Simulates 3D diffusion
* `make` - The corresponding make file (notice that the `bout_runner`
  calls this, so no make is necessary for the `bout_runner` to work).

### Folders:

* `data` - Contains a `BOUT.inp` file
* `MMS` - Contains the `BOUT.inp` file for the MMS runs
* `pre_and_post_processing` - Contains the grid generator and the
  post processing functions

### Examples:

* `1-basic_driver.py` - How to use `bout_runners` for a basic run
* `2-run_with_simple_post_processing.py` - How couple `bout_runners`
  to a post processing routine
* `3-override_BOUTinp.py` - Use `bout_runners` to override settings
  in `BOUT.inp`
* `4-run_with_combinations.py` - Use `bout_runners` to change
  several settings at once
* `5-run_with_grid_files.py` - Run `bout_runners` with a grid file
* `6a-run_with_MMS_post_processing_specify_numbers.py` - Use
  `bout_runners` to MMS the program
* `6b-run_with_MMS_post_processing_grid_file.py` - The same as `6a`,
  but using a gride file
* `7-basic_PBS_run.py` - Submit jobs to a cluster using
  `bout_runners`
* `8-PBS_run_extra_option.py` - Set the `PBS` option using
  `bout_runners`
* `9-PBS_with_MMS_post_processing_grid_file.py` - As 6b, but on a
  cluster
* `10-restart_with_resize.py` - Restart a run and re-size the grid
  using `bout_runners`
* `11-restart_with_scan.py` - Use `bout_runners` to restart runs
  belonging to a parameter scan
* `12-PBS_restart_with_waiting.py` - Runs where the restart waits for jobs to
  finish
*  `13-restart_w_add_noise.py` - Adds noise to a restart run
