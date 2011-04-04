;; Main routine to generate shifted circle equilibria

PRO generate
  set_params   ; set defaults

  REPEAT BEGIN
      ; Modify the input
      set_params, /set, /noreset
      
      ; Display the equilibrium
      GET_DELTA,/sh,/nores

  ENDREP UNTIL get_yesno("Is this ok?") EQ 1

  rmin = get_float("Minimum rho");
  rmax = get_float("Maximum rho");
  nr = get_integer("Number of radial points:")
  ntheta = get_integer("Number of poloidal points:")

  
  SET_MESH_SHCIR,/pl, /nores, nth=ntheta, nr=nr,/exp, $
    rhomin=rmin, rhomax=rmax
  
END
