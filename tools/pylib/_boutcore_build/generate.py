import sys
import os
import jinja2

from scan_enums import enums


class Field(object):
    """Abstracts over BoutReals and Field2D/3D/Perps

    Provides some helper functions for writing function signatures and
    passing data

    """

    def __init__(self, field_type, dimensions, fdd):
        # C++ type of the field, e.g. Field3D
        self.field_type = field_type
        # array: dimensions of the field
        self.dimensions = dimensions
        # Short identifier
        self.fdd = fdd
        self.ddd = fdd

    def __eq__(self, other):
        try:
            return self.field_type == other.field_type
        except AttributeError:
            return self.field_type == other

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return self.field_type

    @property
    def header(self):
        return self.field_type.lower()

    @property
    def ndims(self):
        return len(self.dimensions)

    def makelist(self, fmt):
        """make an index list based on `fmt` and replace $d and $i"""
        return ", ".join(
            [
                fmt.replace("$d", d).replace("$i", str(i))
                for i, d in enumerate(self.dimensions)
            ]
        )


class Vector(object):
    def __init__(self, vec_type, vdd, field, fdd, header):
        self.vec_type = vec_type
        self.vdd = vdd
        self.field = Field(field, [], fdd)
        self.fdd = fdd
        self.ddd = vdd
        self.header = header

    def __repr__(self):
        return self.vec_type


ops = [("add", "+"), ("mul", "*"), ("truediv", "/"), ("div", "/"), ("sub", "-")]

field2d = Field("Field2D", ["x", "y"], "f2d")
field3d = Field("Field3D", ["x", "y", "z"], "f3d")
fieldperp = Field("FieldPerp", ["x", "z"], "fperp")
vector3d = Vector("Vector3D", "v3d", "Field3D", "f3d", "vector3d")
vector2d = Vector("Vector2D", "v2d", "Field2D", "f2d", "vector2d")
vecs = [vector3d, vector2d]
fields = [field3d, field2d]

if __name__ == "__main__":
    inf = sys.argv[1]
    outf = sys.argv[2]
    tmpf = f"{outf}.tmp"

    with open(tmpf, "w") as out:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader("."), trim_blocks=True)

        template = env.get_template(inf)

        args = dict(
            field3d=field3d,
            field2d=field2d,
            fieldperp=fieldperp,
            ops=ops,
            vecs=vecs,
            fields=fields,
            vector3d=vector3d,
            vector2d=vector2d,
            enums=enums,
        )
        out.write(template.render(**args) + "\n")

    os.rename(tmpf, outf)
