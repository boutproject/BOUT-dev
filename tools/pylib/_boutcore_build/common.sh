

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

# What is the metric type?
coords_is() {
    if grep '^#define BOUT_USE_METRIC_3D 1' ../../../include/bout/build_defines.hxx -q
    then
	echo "f3d"
    else
	echo "f2d"
    fi ; }

name_of_fdd() {
    case $1 in
	f3d)
	    echo "Field3D"
	    ;;
	f2d)
	    echo "Field2D"
	    ;;
	*)
	    echo "unexpected fdd"
	    exit 4;
    esac ; }

# Set variables for Vectors and Fields
setvars() {
    case $1 in
	Vector2D)
	    fdd=$(coords_is)
	    field=$(name_of_fdd $fdd)
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

metric_fdd=$(coords_is)
metric_field=$(name_of_fdd $metric_fdd)
