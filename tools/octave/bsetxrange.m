# Adjust X range to be read
function d = bsetxrange(desc, xmin, xmax)
  narg = nargin();
  if (narg < 1)
    d = 0;
    return;
  endif
  
  d = desc;
  
  if (narg < 2)
    # Reset X range
    d.xrange = [1, desc.nx];
    return;
  elseif (narg < 3)
    # Only one number
    xmax = xmin;
  endif

  if(xmin > xmax)
    # Swap around
    tmp = xmin;
    xmin = xmax;
    xmax = tmp;
  endif

  if (xmin < 1)
    xmin = 1;
  endif
  if (xmin > desc.nx)
    xmin = desc.nx;
  endif
  
  if (xmax < xmin)
    xmax = xmin;
  endif
  if (xmax >= desc.nx)
    xmax = desc.nx;
  endif
  
  d.xrange = [xmin, xmax];
endfunction
