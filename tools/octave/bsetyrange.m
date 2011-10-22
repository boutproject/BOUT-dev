# Adjust Y range to be read
function d = bsetyrange(desc, ymin, ymax)
  narg = nargin();
  if (narg < 1)
    d = 0;
    return;
  endif
  
  d = desc;
  
  if (narg < 2)
    # Reset Y range
    d.yrange = [1, desc.ny];
    return;
  elseif (narg < 3)
    # Only one number
    ymax = ymin;
  endif

  if(ymin > ymax)
    # Swap around
    tmp = ymin;
    ymin = ymax;
    ymax = tmp;
  endif

  if (ymin < 1)
    ymin = 1;
  endif
  if (ymin > desc.ny)
    ymin = desc.ny;
  endif
  
  if (ymax < ymin)
    ymax = ymin;
  endif
  if (ymax > desc.ny)
    ymax = desc.ny;
  endif
  
  d.yrange = [ymin, ymax];
endfunction
