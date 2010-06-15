function differential_theta, x, y
;
; Expression of differentials
; dr/ds = Br/B
; dz/ds = Bz/B
;
; Inputs: x->s (distance along theta line)
; y[0:1] -> [r,z] (coords)
;------------------------------------------------;

  brnow=Get_Br(y[0],y[1],n=0)  
  bznow=Get_Bz(y[0],y[1],n=0)
  bnow=sqrt(brnow^2+bznow^2)
  res=[brnow/bnow, bznow/bnow]

return, res
end


function differential_rho, x, y
;
; Expression of differentials
; dr/dpsi = -Bz/B^2
; dz/dpsi = Br/B^2
;
; Inputs: x->s (distance along rho line)
; y[0:1] -> [r,z] (coords)
;------------------------------------------------;

  brnow=Get_Br(y[0],y[1])  
  bznow=Get_Bz(y[0],y[1])
  bnow=sqrt(brnow^2+bznow^2)
  res=[-bznow/bnow^2, brnow/bnow^2]

return, res
end


function theta_line, r0,z0, ds, nstep, z_end=z_end
;
; Integrate along the poloidal (theta) line
; z_end is the termination z value
;------------------------------------------------;

  res={s:dblarr(nstep), r:dblarr(nstep), z:dblarr(nstep)} ;-integration path
  
  ;-initial vals
  s=0d0
  rzold=[r0,z0]

  res.s[0]=0d0
  res.r[0]=r0
  res.z[0]=z0

  for i=1,nstep-1 do begin
      rznew=LSODE(rzold,s,ds,'differential_theta')
      rzold=rznew
      s=s+ds

      res.s[i]=s
      res.r[i]=rznew[0]
      res.z[i]=rznew[1]


      if keyword_set(Z_END) then begin
          

          ;-terminate at z_end
          sLeft=res.s[i-1]
          sRight=res.s[i]
          rLeft=res.r[i-1]
          rRight=res.r[i]
          zLeft=res.z[i-1]
          zRight=res.z[i]

          delzLeft=zLeft-z_end
          delzRight=zRight-z_end

          if (delzLeft*delzRight le 0.) then begin              
              s_end=sLeft + (sRight-sLeft)/(1d0 + abs(delzRight/delzLeft))
              r_end=rLeft + (rRight-rLeft)*(s_end-sLeft)/(sRight-sLeft)

              res.s[i]=s_end
              res.r[i]=r_end
              res.z[i]=z_end

              break ;-exit the loop
          endif

      endif

  endfor


;-truncate the rest of arrays, use ABS(s) so that s is minimum at the x-point
ilast=MIN([i,nstep-1])
return, {r:res.r[0:ilast], z:res.z[0:ilast], s:ABS(res.s[0:ilast])}
end


function rho_line, r0,z0, dpsi, nstep, psi_end=psi_end
;
; Integrate along the radial (rho) line
;------------------------------------------------;
  
  ;-initial vals
  psi0=Get_psi(r0,z0)
  rzold=[r0,z0]

  ;-make a step of size dpsi
  if (1) then begin
      rznew=LSODE(rzold,psi0,dpsi,'differential_rho')
  endif else begin
      dydx=Differential_Rho(psi0, rzold)
      rznew=RK4(rzold, dydx, psi0, dpsi, 'differential_rho')
  endelse


return, {r:rznew[0], z:rznew[1], psi:psi0+dpsi}
end
