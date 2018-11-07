
for ord in o2 o3 o4
do
    for bndry in dirichlet_ neumann_ free_ dirichlet_nu_ neumann_nu_ free_nu_
    do
	out=$(./BoundaryOp_timing -q -q -q f:bndry_all=${bndry}$ord) && echo $out
    done
done
