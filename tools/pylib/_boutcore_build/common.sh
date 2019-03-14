

# make a list of the format "nx, ny, nz" or something similar
# 'n$d' gets expaneded to "nx, ny, nz"
# 'f[$i]' to "f[0], f[1], f[2]"
# NB: use single quotes to prevent early expansion
# First argument        : the expression to expand, e.g. 'n$d' or
#                         'f[$i]'
# Second argument (opt) : the spacing between the
#                         elements in the list (default ', ')
makelist () {
    format="$1"
    test $# -gt 1 &&
	spacing="$2" || spacing=", "
    start="" ; i=0
    for d in ${dims[@]} ; do
	echo -n "$start" ;
	eval echo -n "\"$format\""
	start="$spacing" ; i=$(( i + 1 ))
    done ; }

# Set variables for Vectors and Fields
setvars() {
    case $1 in
	Vector2D)
	    field=Field2D
	    fdd=f2d
	    vdd=v2d
	    fheader=vector2d
	    ;;
	Vector3D)
	    field=Field3D
	    fdd=f3d
	    vdd=v3d
	    fheader=vector3d
	    ;;
	Field3D)
	    fdd="f3d"
	    ndim=3
	    dims=(x y z)
	    fheader="field3d"
	    ;;
	Field2D)
	    fdd="f2d"
	    ndim=2
	    dims=(x y)
	    fheader="field2d"
	    ;;
	*)
	    echo "$1 - not implemented"
	    exit 3
    esac
}

vecs="Vector3D Vector2D"

fields="Field3D Field2D"
