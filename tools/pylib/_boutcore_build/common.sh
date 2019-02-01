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
	*)
	    echo "$1 - not implemented"
	    exit 3
    esac
}

vecs="Vector3D Vector2D"

fields="Field3D Field2D"
