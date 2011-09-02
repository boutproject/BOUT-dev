; set_plasma.pro
;
; Author:       Ilon Joseph
; Began:        2011/08/22
; Modified:     2011/09/02 
;
; given a reference state, set_plasma calculates plasma parameters and
; creates a structure to hold important plasma property calculations

pro set_plasma, p, cyclone=cyclone
; UNITS: cgs + ev

  common constants, c
  physical_constants

  zeff = 1.0
  aeff = 2.0
  ni = 1.0e14
  te = 1.0e3
  ti = 1.0e3
  B0 = 6.4626e4 ; 2.0e4
  R0 = 5.0e2    ; 2.0e2
  a0 = 1.0e2    ; 50.0 
  eps0 = a0/R0
  q0 = 3.0
  dq = 1.0
  dr = 0.1*a0
  s0 = (dq/q0)/(dr/a0)


  etai = 1.
  etae = 1.
  Ln  = 10.0
  Lti = Ln/etai
  Lte = Ln/etae


  if keyword_set(cyclone) then begin
    zeff = 1.0
    aeff = 2.0
    ni = 1.0e14
    te = 1.0e3
    ti = 1.0e3
    etai = 3.114
    etae = etai     ;? is this correct?

    R0 = 4.0
    q0 = 1.4
    s0 = 0.776      ; = dlog(q)/dlog(r)
    eps0 = 0.18     ; = a0/R0
    a0 = eps0*R0

    Rnorm = 6.92    ; = R/Lti
    Lti = R/Rnorm
    Ln  = etai*Lti
    Lte = Ln/etae

    rho_norm = 0.01 ; = rho_ci/Lti 
    rho_ci = rho_norm*Lti 

    B0 = sqrt(2*ti*c.meff/c.ee)*1e-4
  endif

  Bp0  = B0*a0/(q0*R0)
  Lq   = s0*a0
  Lcon = q0*R0
  Lshear = Lcon/s0

  omg0 = 1.0e4
  n0 = 1.0
  m0 = q0*n0
  ky = m0/a0
  kz = n0/R0
  kprime = ky/Lshear
 
  diff_ni = 1.0e4
  diff_vi = diff_ni
  diff_ve = diff_ni
  diff_ti = diff_ni
  diff_te = diff_ni

  diff = {ni:diff_ni,vi:diff_vi,ve:diff_ve,ti:diff_ti,te:diff_te} 


  r = {major:R0, minor:a0, eps:eps0, $
       q:Lq, con:Lcon, shear:Lshear, $
       ni:Ln, ti:Lti, te:Lte, etai:etai, etae:etae}

  b = {tor:B0, pol:Bp0, mag:B0, q:q0, s:s0}


  mode = {ntor:n0, mpol:m0, ky:ky, kz:kz, kprime:kprime, omg:omg0}
 

  p = { aeff:aeff, zeff:zeff, ni:ni, ti:ti, te:te, loglam:c.loglam, $
        b:b, r:r, omg_eb:omg0, diff0:diff, mode:mode}

  plasma_params, p
  recon_params, p
end

pro plasma_params, p 
  common constants, c
  physical_constants

; need to preload 
; p.ne, p.te, p.ti, p.Rmajor, p.rminor
; p.B0, p.q0, p.s0. p.loglam

  if not keyword_set(nz) then nz=1.0
  if not keyword_set(nz) then omega=1.0


; Velocity
  valfven = sqrt(p.B.mag^2/(p.ni*c.M.i));
  vti = sqrt(p.ti*c.kb/c.M.i)
  vte = sqrt(p.te*c.kb/c.M.e)
  vsound = sqrt((p.ti+p.te)*c.kb/c.M.i)

; Microscopic frequency & distances
  omg_ci = c.Zeff*c.ee*p.B.mag/(c.M.i*c.clight)
  omg_ce =        c.ee*p.B.mag/(c.M.e*c.clight)
  rho_sound = vsound/omg_ci
  rho_ci = vti/omg_ci
  rho_ce = vte/omg_ce

  omg_pi = sqrt(4.0*!pi*c.ee^2*p.ni/c.m.i)
  omg_pe = sqrt(4.0*!pi*c.ee^2*p.ni/c.m.e)
  rho_pi = vti/omg_pi
  rho_pe = vte/omg_pe
  rho_di = c.clight/omg_pi
  rho_de = c.clight/omg_pe


; Drift frequencies in Toroidal Freq norm
  omg_ni = p.mode.ky*vti*rho_ci/p.r.ni * p.b.q
  omg_ne = p.mode.ky*vte*rho_ce/p.r.ni * p.b.q
  omg_ti = p.mode.ky*vti*rho_ci/p.r.ti * p.b.q
  omg_te = p.mode.ky*vte*rho_ce/p.r.te * p.b.q

; Diffusivity
  diff_res  = Diffusivity_Res(p.Te,loglam=p.loglam)
  diff_bohm = Diffusivity_Bohm(p.Te,p.B.mag)
  diff = create_struct(p.diff0, {res:diff_res, bohm:diff_bohm})

; Transport times
  tau_alf = p.R.major/valfven
  tau_shear = p.R.shear/valfven

  tau_ee = Collision_Time_e(p.ni,p.te,loglam=p.loglam) 
  tau_ii = Collision_Time_i(p.ni,p.ti,loglam=p.loglam) 

  tau_res = p.R.minor^2/diff_res
  tau_bohm = p.R.minor^2/diff_bohm
  tau_ni = p.R.minor^2/p.diff0.ni
  tau_vi = p.R.minor^2/p.diff0.vi
  tau_ti = p.R.minor^2/p.diff0.ti
  tau_te = p.R.minor^2/p.diff0.te

  tau = {alfven:tau_alf, shear:tau_shear, ee:tau_ee, ii:tau_ii, $
         bohm:tau_bohm, $
         ni:tau_ni, vi:tau_vi, ti:tau_ti, te:tau_te, $
         res:tau_res}

  omg = create_struct({ce:omg_ce, pe:omg_pe,  pi:omg_pi,  ci:omg_ci, $
         eb:p.omg_eb, dni:omg_ni, dne:omg_ne , dti:omg_ti, dte:omg_te  })

  rho = {ce:rho_ce, ci:rho_ci, sound:rho_sound, $
         pe:rho_pe, pi:rho_pi, $
         de:rho_de, di:rho_di  }

  vel = {light:c.clight, te:vte,  alfven:valfven, sound:vsound, ti:vti}

  cs = sqrt(c.kb*p.te/c.m.i)
  rhos = cs/omg_ci
  drift = {omg_ci:omg_ci, cs:cs, rhos:rhos, rhocs:rhos*cs}

  p = create_struct(p,{rho:rho, vel:vel, omg:omg, tau:tau, diff:diff, drift:drift})

end


pro recon_params, p

; Dimensionless #'s
  Sra = p.tau.res/p.tau.alfven
  Sva = p.tau.vi /p.tau.alfven
  Sna = p.tau.vi /p.tau.alfven
  Srs = p.tau.res/p.tau.shear
  Svs = p.tau.vi /p.tau.shear
  Sns = p.tau.ni /p.tau.shear
  Prv = p.tau.res/p.tau.vi
  Prn = p.tau.res/p.tau.ni

  Slund   = {ra:Sra, va:Sva, na:Sna, rs:Srs, vs:Svs, ns:Sns}
  Prandtl = {rv:Prv, rn:Prn}

  i = complex(0.,1.)
  mpol = p.mode.mpol
  ntor = p.mode.ntor

  qmg_eb = p.omg.eb *p.tau.shear*Srs^(1/3.)  
  qmg_ni = p.omg.dni*p.tau.shear*Srs^(1/3.)  
  qmg_ne = p.omg.dne*p.tau.shear*Srs^(1/3.)  
  qmg_ti = p.omg.dti*p.tau.shear*Srs^(1/3.)  
  qmg_te = p.omg.dte*p.tau.shear*Srs^(1/3.)  
  Qmg = {eb:qmg_eb, dni:qmg_ni, dne:qmg_ne , dti:qmg_ti, dte:qmg_te, rec:Srs^(-1/3.)}


; Visco-Resistive
  del_vr = Prv^(1/6.) * Srs^(-1/3.)
  rho_vr = abs(del_vr)/p.mode.ky
  rci_vr = p.rho.ci/rho_vr
  delta_vr = 2.104*i*qmg_eb*Prv^(1/3.)/del_vr
  trans_vr = 1./(1.+delta_vr/(2.*mpol))
  qmg_vr = 2.*mpol/abs(delta_vr/qmg_eb) 
  omg_vr = qmg_vr*Srs^(-1/3.)/p.tau.shear 

  VR = {del:del_vr, rho:rho_vr, rci:Rci_vr, delta:delta_vr, trans:trans_vr, $
        qmg:qmg_vr, omg:omg_vr}

; Resistive-Inertial
  del_ri = qmg_eb^(1/4.0) * Srs^(-1/3.)
  rho_ri = abs(del_ri)/p.mode.ky
  rci_ri = p.rho.ci/rho_ri 
  delta_ri = 2.124*exp(i*!pi*5./8.)*qmg_eb/del_ri 
  trans_ri = 1/(1.+delta_ri/(2.*mpol))
  qmg_ri = (2.*mpol/abs(delta_ri/qmg_eb^(5./4.)))^(4./5.) 
  omg_ri = qmg_ri*Srs^(-1/3.)/p.tau.shear 

  RI = {del:del_ri, rho:rho_ri, rci:Rci_ri, delta:delta_ri, trans:trans_ri, $
        qmg:qmg_ri, omg:omg_ri}
 
; Visco-Inertial
  del_vi = qmg_eb^(1/4.)*Prv^(1/4.) * Srs^(-1/3.)
  delta_vi = - 4.647*exp(-i*!pi/8.0)/del_vi
  trans_vi = 1/(1.+delta_vi/(2.0*mpol))
  rho_vi = abs(del_vi)/p.mode.ky
  rci_vi = p.rho.ci/rho_vi
  qmg_vi = ( abs(delta_vi*qmg_eb^(1/4.))/(2.*mpol) )^4 
  omg_vi = qmg_vi*Srs^(-1/3.)/p.tau.shear 

  VI = {del:del_vi, rho:rho_vi, rci:Rci_vi, delta:delta_vi, trans:trans_vi, $
        qmg:qmg_vi, omg:omg_vi}

; Inertial
  del_ii = qmg_eb * Srs^(-1/3.)
  rho_ii = abs(del_ii)/p.mode.ky
  rci_ii = p.rho.ci/rho_ii
  delta_ii = 3.142*i/del_ii 
  trans_ii = 1./(1.+delta_ii/(2.*mpol))
  qmg_ii = abs(delta_ri*qmg_eb)/(2.*mpol)
  omg_ii = qmg_ii*Srs^(-1/3.)/p.tau.shear
 
  II = {del:del_ii, rho:rho_ii, rci:Rci_ii, delta:delta_ii, trans:trans_ii, $
        qmg:qmg_ii, omg:omg_ii}


  recon = {Slund:Slund, Prandtl:Prandtl, Qmg:Qmg, RI:RI, VR:VR, VI:VI, II:II}
  p = create_struct(p,{recon:recon})

end



