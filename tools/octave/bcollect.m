# Collect data from BOUT++ output files
#
# Tested on Octave 3.2 
# Needs the octcdf library for NetCDF file reading
# 
# Usage:
# =====
#
#  f = bcollect()  # path is "."
#
#  f = bsetxrange(f, 1, 10) # Set ranges
#  Same for y, z, and t (NOTE: indexing from 1!)
#  
#  u = bread(f, "U")  # Finally read the variable
#


# Collect metadata from collection of data files
function desc = bcollect(path)
  error("This is currently broken for BOUT++ > v4.0.0. See issue #394")
  narg = nargin();
  if (narg < 1)
    # No path specified, so use current directory
    path="./";
  endif
  
  desc.path = path;
  
  f = netcdf(strcat(path,"/BOUT.dmp.0.nc"), 'r');
  
  desc.version = f{'BOUT_VERSION'}(:);
  desc.mxsub   = f{'MXSUB'}(:);
  desc.mysub   = f{'MYSUB'}(:);
  desc.mz      = f{'MZ'}(:);
  desc.myg     = f{'MYG'}(:);
  desc.t_array = f{'t_array'}(:);
  
  desc.nxpe = f{'NXPE'}(:);
  desc.nype = f{'NYPE'}(:);
  desc.npe  = desc.nxpe * desc.nype;
  desc.mxg  = f{'MXG'}(:);
  
  desc.nx = desc.nxpe * desc.mxsub + 2*desc.mxg;
  desc.ny = desc.mysub * desc.nype;
  
  printf("Size of the grid: %d x %d x %d\n", desc.nx, desc.ny, desc.mz);
  printf("In each file: %d x %d x %d\n", desc.mxsub, desc.mysub, desc.mz);
  
  # Set starting range
  desc.xrange = [1, desc.nx];
  desc.yrange = [1, desc.ny];
  desc.zrange = [1, desc.mz];
  desc.trange = [1, length(desc.t_array)];
  
  close(f);
endfunction
