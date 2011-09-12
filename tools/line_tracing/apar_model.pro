; Model function for Apar
; Input is x,y,z coordinates
; Output is structure ={f:Apar(x,y,z),x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_model, x, y, z 
  sFactor=1e-3

  ;xmid calculated from g.psixy[xs_s0*512], NOT xs_s0*delta_psi+psi_min !!!
  xmid=-0.164009 ;; using xs_s0=0.5
  ;xmid=-0.200792 ;; using xs_s0=0.45
  ymid=!PI

  alpx=185.637; =(1/(delta_x*xs_wd))^2, delta_x=0.733953, xs_wd=0.1
  ;alpx=20.626; =(1/(delta_x*xs_wd))^2, delta_x=0.733953, xs_wd=0.3
  alpy=2.53303; =(1/(delta_y*ys_wd))^2, delta_y=2*!PI, ys_wd=0.1
  ;alpy=0

  ;-if BOUT domain is 2PI/5 then in the perturbation
  ;only nz=5,10,15,.. etc can be present
  nz=5.0d

  Apar=SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)
  dAdx=   SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)*(-2*(x-xMid)*alpx)
  dAdy=   SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)*(-2*(y-yMid)*alpy)
  dAdz=nz*COS(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)

  res={f:sFactor*Apar,x:sFactor*dAdx, y:sFactor*dAdy, z:sFactor*dAdz}

  RETURN, res
;
END
