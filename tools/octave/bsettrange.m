# Adjust T range to be read
function d = bsettrange(desc, tmin, tmax)
  narg = nargin();
  if (narg < 1)
    d = 0;
    return;
  endif
  
  d = desc;
  
  nt = length(d.t_array)
  
  if (narg < 2)
    # Reset Y range
    d.trange = [1, nt];
    return;
  elseif (narg < 3)
    # Only one number
    tmax = tmin;
  endif

  if(tmin > tmax)
    # Swap around
    tmp = tmin;
    tmin = tmax;
    tmax = tmp;
  endif

  if (tmin < 1)
    tmin = 1;
  endif
  if (tmin > nt)
    tmin = nt;
  endif
  
  if (tmax < tmin)
    tmax = tmin;
  endif
  if (tmax > nt)
    tmax = nt;
  endif
  
  d.trange = [tmin, tmax];
endfunction
