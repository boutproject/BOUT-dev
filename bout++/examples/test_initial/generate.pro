; Create an input file for testing. 
; Needs to have an x-point region, hence needs MPI

nx = 12 ; 4 for guard cells, so 8 in domain
ny = 16 ; Going to need 4 processors

ixseps = 6
jyseps1_1 = 3
jyseps1_2 = 6
jyseps2_1 = 6
jyseps2_2 = 11

; Angle for twist-shift location
ShiftAngle = (1.0 + FINDGEN(nx))/FLOAT(nx) * 10.*!PI

f = file_open('test_initial.grd.nc', /create)

status = file_write(f, 'nx', nx)
status = file_write(f, 'ny', ny)
status = file_write(f, 'ixseps1', ixseps)
status = file_write(f, 'ixseps2', ixseps)
status = file_write(f, 'jyseps1_1', jyseps1_1)
status = file_write(f, 'jyseps1_2', jyseps1_2)
status = file_write(f, 'jyseps2_1', jyseps2_1)
status = file_write(f, 'jyseps2_2', jyseps2_2)
status = file_write(f, 'ShiftAngle', ShiftAngle)

file_close, f

exit
