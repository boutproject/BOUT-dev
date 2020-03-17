
g = read_neqdsk("efit/neqdsk")

R = REFORM(g.r[*,0])
Z = REFORM(g.z[0,*])

boundary=DBLARR(2,4)
boundary[0,*] = [1.0D, 1.0D, 2.5D, 2.5D]
boundary[1,*] = [-1.4D, 1.4D, 1.4D, -1.4D]

;boundary = TRANSPOSE([[g.xlim], [g.ylim]])

; Find a field-aligned mesh
mesh = create_grid(g.psi, R, Z, boundary=boundary, /strict, /simple)

; Create a grid from the mesh and g-file data
grid = process_grid(g, mesh)

