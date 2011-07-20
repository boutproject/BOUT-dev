pro IntegrationStep, x0, y0, z0, x1, y1, z1
;
;
; Perform integration from y0 to y1
;============================================;


  V1=[x0,z0] ;;-dependent variable
  y=y0       ;;-independent variable
  delta_y=y1-y0


  if (0) then begin
      ;;-lsode will attempt to jump beyond the branch-cut
      V2=LSODE(V1,y,delta_y, 'differential', atol=1e-7, rtol=1e-7)
  endif else begin

      ;;-use RK4 with Nstep substeps
      Nstep=100
      dy=delta_y/Nstep

      for i=0,Nstep-1 do begin
          dvdy=DIFFERENTIAL(y,V1)
          V2=RK4(V1, dvdy, y, dy, 'differential', /DOUBLE)
          y = y + dy
          V1 = V2
      endfor

  endelse


  x1=V2[0]
  z1=V2[1]

;
;
end
