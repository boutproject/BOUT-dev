; Generates an example grid for straight linear device like LAPD

set_mesh_cyl, /export, Nr=54, Nz=64, rMin=0.15, rMax=0.45, $
              ni0=2.5e18, te0=5., Bz0=0.08, ni_profile_type=1, $
              phi_profile_type=0, phi0V=0.0, /NOPLOTS


!path=!path+":../tokamak_grids/all/"

uedge2bout, /default    ; Default keyword uses default settings for everything


; Should now have a file "uedge.grd.nc" which can be input to BOUT++

