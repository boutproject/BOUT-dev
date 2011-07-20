FUNCTION interpZApar, xi, yi, zi
  COMMON BDATA, bd
  
  nz=bd.nz
  
  zm = ((FLOOR(zi) MOD nz) + nz) MOD nz
  zp = (zm + 1) MOD nz
  
  d = zi - FLOOR(zi)
  
  RETURN, (1.-d)*bd.apar[xi,yi,zm] + d*bd.apar[xi,yi,zp]
END


