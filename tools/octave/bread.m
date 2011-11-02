
# Read variable
function data = bread(desc, varname)
  # Open the first file
  f = netcdf(strcat(desc.path,"/BOUT.dmp.0.nc"), 'r');
  
  dims = size(f{varname});
  nd = length(dims);
  
  printf("Number of dimensions: %d\n", nd)

  if (nd < 2) 
    # 0 or 1-D, just read from file
    data = f{'varname'}(:);
    close(f);
    return;
  endif
  
  nx = desc.xrange(2) - desc.xrange(1) + 1
  ny = desc.yrange(2) - desc.yrange(1) + 1
  nz = desc.zrange(2) - desc.zrange(1) + 1
  nt = desc.trange(2) - desc.trange(1) + 1
  
  # Create an empty array
  data = zeros(nt, nx, ny, nz);
  
  close(f);
  
  for i = 0:(desc.npe-1)
    # Get X and Y processor indices
    pe_yind = int32( i / desc.nxpe );
    pe_xind = rem(i, desc.nxpe);
    
    # Get local ranges
    ymin = desc.yrange(1) - pe_yind*desc.mysub + desc.myg;
    ymax = desc.yrange(2) - pe_yind*desc.mysub + desc.myg;
    
    xmin = desc.xrange(1) - pe_xind*desc.mxsub;
    xmax = desc.xrange(2) - pe_xind*desc.mxsub;
    
    inrange = 1;
    
    if ((ymin >= (desc.mysub + desc.myg)) || (ymax < desc.myg))
      inrange = 0; # Out of Y range
    endif
    
    if (ymin < desc.myg)
      ymin = desc.myg;
    endif
    if (ymax >= desc.mysub + desc.myg)
      ymax = desc.myg + desc.mysub - 1
    endif
    
    # Check lower X boundary
    if (pe_xind == 0)
      # Keep inner boundary
      if (xmax < 0)
	inrange = 0;
      endif
      if (xmin < 0)
	xmin = 0;
      endif
    else
      if (xmax < desc.mxg)
	inrange = 0;
      endif
      if (xmin < desc.mxg)
	xmin = desc.mxg;
      endif
    endif
    
    # Check upper boundary
    if (pe_xind == (desc.nxpe - 1))
      # Keeping outer boundary
      if (xmin >= (desc.mxsub + 2*desc.mxg))
	inrange = 0;
      endif
      if (xmax >= (desc.mxsub + 2*desc.mxg))
	xmax = desc.mxsub + 2*desc.mxg - 1;
      endif
    else
      if (xmin >= (desc.mxsub + desc.mxg))
	inrange = 0;
      endif
      if (xmax >= (desc.mxsub + desc.mxg))
	xmax = desc.mxsub + desc.mxg - 1;
      endif
    endif
    
    if (inrange == 1)
      # Calculate global indices
      xgmin = xmin + pe_xind * desc.mxsub;
      xgmax = xmax + pe_xind * desc.mxsub;
      
      ygmin = ymin + pe_yind*desc.mysub - desc.myg;
      ygmax = ymax + pe_yind*desc.mysub - desc.myg;
      
      # Open the file
      filename = strcat(desc.path,"/BOUT.dmp.", num2str(i), ".nc")
      f = netcdf(filename, 'r');
      
      if (nd == 4)
	# Print local to global ranges
	printf("[%d:%d, %d:%d, %d:%d, %d:%d] -> [%d:%d, %d:%d, %d:%d, %d:%d]\n",
	       desc.trange(1),desc.trange(2),
	       xmin,xmax,ymin,ymax,desc.zrange(1),desc.zrange(2),
	       desc.trange(1),desc.trange(2),
	       xgmin,xgmax,ygmin,ygmax,desc.zrange(1),desc.zrange(2))
	
	tmp = f{varname}(desc.trange(1):desc.trange(2),xmin:xmax,ymin:ymax,desc.zrange(1):desc.zrange(2));
	
	data(desc.trange(1):desc.trange(2), xgmin:xgmax, ygmin:ygmax, desc.zrange(1):desc.zrange(2)) = tmp
      elseif (nd == 3)
	# Could be TXY or XYZ
	if (dims(3) == desc.mz)
	  # XYZ
	  tmp = f{varname}(xmin:xmax,ymin:ymax,desc.zrange(1):desc.zrange(2));
	  data(xgmin:xgmax, ygmin:ygmax, desc.zrange(1):desc.zrange(2)) = tmp
	else
	  # TXY
	  tmp = f{varname}(desc.trange(1):desc.trange(2),xmin:xmax,ymin:ymax);
	  data(desc.trange(1):desc.trange(2), xgmin:xgmax, ygmin:ygmax) = tmp
	endif
      elseif (nd == 2)
	# Assume XY
	tmp = f{varname}(xmin:xmax,ymin:ymax);
	data(xgmin:xgmax, ygmin:ygmax) = tmp
      endif 
      
      close(f);
    endif
  endfor
endfunction
