; collisions.pro
;
; Author:       Ilon Joseph
; Began:        2011/08/22
; Modified:     2011/09/02 
;
;
; defines useful functions associated with collisions and transport


function Collision_Time_e, N, T, loglam=loglam
 ; = 3.44e5 sec x T^1.5/N*loglam 

  common constants, c
  physical_constants
  if not keyword_set(loglam) then loglam = c.loglam
  taue = (3.0*sqrt(c.M.e)*(c.kb*T)^1.5)/(4.0*sqrt(2.0*!pi)*n*loglam*c.ee^4)
  return, taue
end

function Collision_Time_i, N, T, loglam=loglam
 ; = 2.09e7 sec x T^1.5*A^0.5/N*loglam

  common constants, c
  physical_constants
  if not keyword_set(loglam) then loglam = c.loglam
  taui = (3.0*sqrt(c.M.i)*(c.kb*T)^1.5)/(4.0*sqrt(!pi)*N*loglam*c.ee^4)
  return, taui
end

function Conductivity_Perp, T, loglam=loglam
; = 8.70e13 sec x T^1.5/loglam/Zeff

  common constants, c
  physical_constants
  if not keyword_set(loglam) then loglam = c.loglam
  sigma_perp = Collision_Time_e(1.0,T,loglam=loglam)*c.ee^2/c.M.e
  return, sigma_perp
end

function Conductivity_Par, T, loglam=loglam
; 1.71e14 sec x T^1.5/loglam/Zeff

  sigma_par = 1.96*Conductivity_Par(T,loglam)
  return, sigma_par
end

function Diffusivity_Res, T, loglam=loglam
; = 8.21e5 cm2/s loglam*Zeff/T^1.5

  common constants, c
  physical_constants
  diff_res = c.clight^2/(4.0*!pi*Conductivity_Perp(T,loglam=loglam) )
  return, diff_res
end

function Diffusivity_Bohm, T, B
; = 6.25e6 cm2/s x T/B

  common constants, c
  physical_constants
  dbohm = c.clight*c.kb*T/(16*c.ee*B)
  return, dbohm
end
