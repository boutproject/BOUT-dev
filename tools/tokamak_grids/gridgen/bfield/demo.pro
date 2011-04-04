$rm -f efitxs.so; ln -s ../efit/efitxs.so
$m -f  aeqdsk;    ln -s ../efit/aeqdsk
$m -f  neqdsk;    ln -s ../efit/neqdsk

.run dct2d.pro  gg_gint.pro  gg_read_efit.pro, gg_integrator

READ_EFIT, psi, r, z, /load

tek_color
CONTOUR, psi,r,z,/iso, nlev=10

tl=theta_line(2.,1.,3e-2, 100) ;;-integration along flux surface
oplot, tl.r, tl.z, col=2, psym=3
