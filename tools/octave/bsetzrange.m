# Adjust Z range to be read
function d = bsetzrange(desc, zmin, zmax)
  narg = nargin();
  if (narg < 1)
    d = 0;
    return;
  endif
  
  d = desc;
  
  if (narg < 2)
    # Reset Y range
    d.zrange = [1, desc.mz];
    return;
  elseif (narg < 3)
    # Only one number
    zmax = zmin;
  endif

  if(zmin > zmax)
    # Swap around
    tmp = zmin;
    zmin = zmax;
    zmax = tmp;
  endif

  if (zmin < 1)
    zmin = 1;
  endif
  if (zmin > desc.mz)
    zmin = desc.mz;
  endif
  
  if (zmax < zmin)
    zmax = zmin;
  endif
  if (zmax > desc.mz)
    zmax = desc.mz;
  endif
  
  d.zrange = [zmin, zmax];
endfunction
