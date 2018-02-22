; Create an input file for testing

nx = 12
ny = 12
nz = 5  ; Number of components

ivar = 1
rvar = !PI
f2d = RANDOMN(1, nx, ny)
f3d = RANDOMN(2, nx, ny, nz)

f = file_open('test_io.grd.nc', /create)

status = file_write(f, 'nx', nx)
status = file_write(f, 'ny', ny)
status = file_write(f, 'ivar', ivar)
status = file_write(f, 'rvar', rvar)
status = file_write(f, 'f2d', f2d)
status = file_write(f, 'f3d', f3d)

file_close, f

exit
