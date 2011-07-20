FUNCTION derivZApar, xi, yi, zi
  COMMON BDATA, bd
  COMMON griddata, g, deltaZtor, Ntor
  
  nz=bd.nz
  
  zm = ((FLOOR(zi) MOD nz) + nz) MOD nz
  zp = (zm + 1) MOD nz
  
  dz = deltaZtor / FLOAT(bd.nz)
  
  RETURN, (bd.apar[xi,yi,zp] - bd.apar[xi,yi,zm])/dz
END


