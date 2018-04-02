#!/usr/bin/python3
""" Code generator for derivs.hxx

"""

import jinja2


file_header="""\
/*!************************************************************************
 * \\file derivs.hxx
 *
 * Basic differential functions
 *
 **************************************************************************
 * Copyright 2010,2017
 *    B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu, D. Schwörer
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __DERIVS_H__
#define __DERIVS_H__

#include "field3d.hxx"
#include "field2d.hxx"
#include "vector3d.hxx"
#include "vector2d.hxx"

#include "bout_types.hxx" // See this for code
""")


class Function(object):
    def __init__(self,name,flux=None,desc=None,latex=None):
        if flux is not None:
            self.name=name
            self.flux=flux
            self.desc=desc
            self.latex=latex
        else:
            # Copy constructor
            self.name=name.name
            self.flux=name.flux
            self.desc=name.desc
            self.latex=name.latex
        if flux:
            self.in_field_desc = "/// @param[in] v       The velocity field\n"
            self.in_field_desc+= "/// @param[in] f       The field of the advected quantity"
        else:
            self.in_field_desc = "/// @param[in] f       The field to be differentiated"
        if flux:
            self.in_sig="const $f &v, const $f &f"
            self.in_field="v, f"
        else:
            self.in_sig="const $f &f"
            self.in_field="f"
    def set(self,d,field):
        copy=Function(self.name.replace('d',d),
                      self.flux,
                      self.desc.replace('$d',d),
                      self.latex.replace('$d_lower',d.lower()))
        copy.in_sig=self.in_sig.replace('$f',field)
        copy.field=field
        copy.d=d
        return copy
    def render(self):
        args=vars(self)
        args['DD']=self.name
        global function_template
        print(function_template.render(**args))

env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'),
                         trim_blocks=True)

function_template = env.get_template("derivs.hxx.in.jinja")

first=Function('DDd',False,
                desc="Calculate first partial derivative in $d",
                latex="\partial / \partial $d_lower")
second=Function('D2Dd2',False,
                desc="Calculate second partial derivative in $d",
                latex="\partial^2 / \partial $d_lower^2")
upwind=Function('VDDd',True,
                desc="For terms of form v * grad(f)",
                latex="v \cdot \partial f / \partial $d_lower")
funcs=[first,
       second,
       Function('D4Dd4',False,
                desc="Calculate forth partial derivative in $d",
                latex="\partial^4 / \partial $d_lower^4"),
       upwind,
       Function('FDDd',True,
                desc="for terms of form div(v * f)",
                latex="\partial (v f) / \partial $d_lower")]

deprecated_methods="""
// Deprecated methods
//
// Calculate first partial derivative in Z
//
//   $\partial / \partial z$
//
// @param[in] f       The field to be differentiated
// @param[in] outloc  The cell location where the result is desired.
//                    If staggered grids is not enabled then this has no effect
// @param[in] method  Differencing method to use. This overrides the default
// @param[in] inc_xbndry  DEPRECATED: use REGION flags
//                    Determines whether the derivative should be calculated in
//                    the X boundaries. This allows mixed operators (e.g.
//                    D2DXDZ) without additional communication


inline const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry) {
  return DDZ(f, outloc, method, inc_xbndry ? RGN_NOY : RGN_NOBNDRY);
}

inline const Field3D DDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc, bool inc_xbndry) {
  return DDZ(f, outloc, method, inc_xbndry ? RGN_NOY : RGN_NOBNDRY);
}

inline const Field3D DDZ(const Field3D &f, DIFF_METHOD method, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, method, inc_xbndry ? RGN_NOY : RGN_NOBNDRY);
}

inline const Field3D DDZ(const Field3D &f, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, DIFF_DEFAULT, inc_xbndry ? RGN_NOY : RGN_NOBNDRY);
}

"""

end_of_file="#endif // __DERIVS_H__"

if __name__ == "__main__":
    print(file_header)

    # Generate normal derivatives for Field3D and Field2D for the
    # various directions
    for fun in funcs:
        if fun.name == 'DDd':
            print("////////// FIRST DERIVATIVES //////////")
        elif fun.name == 'D2Dd2':
            print("////////// SECOND DERIVATIVES //////////")
        elif fun.name == 'D4Dd4':
            print("////////// FORTH DERIVATIVES //////////")
        elif fun.name == 'VDDd':
            print("///////// UPWINDING METHODS /////////////")
        elif fun.name == 'FDDd':
            print("///////// FLUX METHODS /////////////")
        else:
            print("Unhandeled case")
            exit(1)
        for d in ['X', 'Y', 'Z']:
            for field in ['Field3D', 'Field2D']:
                # get copy
                fun.set(d,field).render()


    # Generate header file for the Z derivative of Vector 3D
    first.set('Z',"Vector3D").render()

    # Generate the mixed derivative
    x='x'
    y='y'
    z='z'
    for DD in [[x,y],[x,z],[y,z]]:
        for field in ['Field2D', 'Field3D']:
            cur=second.set('error',field)
            cur.name="D2D%sD%s"%(DD[0].upper(),DD[1].upper())
            cur.desc="Calculate mixed partial derivative in %s and %s"%(DD[0],DD[1])
            cur.latex="\partial^2 / \partial %s \partial %s"%(DD[0],DD[1])
            cur.render()

    # Generate a case of mixed Field2D and Field3D for Z upwinding
    # scheeme
    cur=upwind.set('Z','Field2D')
    cur.in_sig="const Field3D &v, const Field2D &f"
    cur.render()

    print(deprecated_methods)

    print(end_of_file)
