from __future__ import print_function

# BOUT++ Library - Write fluid simulations in curviilinear geometry
# Copyright (C) 2016, 2017, 2018 David Schwörer
#
# Contact: Ben Dudson, bd512@york.ac.uk
#
# This file is part of BOUT++.
#
# BOUT++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BOUT++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.

# Code to read the stencils from BOUT++ and generate a more efficient
# version for the different directions and fields. Also generate the
# interp_to code including forward and backword interpolation assuming
# free boundaries.
# This code has to be called from gen_derivs.py with what functions
# are to be generated, and some informations on what stencils are to
# be used for what.


# Define some constants, functions
from common import *


class Stencil(object):

    def __init__(self, body):
        self.body = body
        self.valid = True
        self.mbf = "main"
        self.guards = 3
        if body[0] != '':
            self.name = self.body[0].split()[1].split("(")[0]
            self.checkFlow()
            self.checkStag()
        else:
            self.name = ''
            self.flow = False
            self.stag = False
        if self.name == 'DDX_CWENO3':
            debug("Skipping DDX_CWENO3 ...")
            self.valid = False
        if self.name[:6] == 'DDX_KT':
            debug("Skipping DDX_KT* ...")
            self.valid = False
        self.setGuards()

    def setGuards(self):
        self.guards = 1
        for l in self.body:
            for off in ['pp', 'mm', 'p2', 'm2']:
                if l.find(off) > -1:
                    self.guards = 2

    def checkStag(self):
        self.stag = (self.body[0].find("stag") > 0)

    def checkFlow(self):
        f0 = self.body[0]
        if f0.find("BoutReal") == 0 and f0.find('stencil') > 0:
            if (f0.find("BoutReal V") == -1 and f0.find("BoutReal F") == -1):
                # They have one stencil
                self.flow = False
            else:
                self.flow = True
        elif f0.find("&fRm") > -1:
            debug(f0, "is invalid")
            self.valid = False
        else:
            raise RuntimeError(
                self.body, "We did not set the stencil type - error parsing")

    def getFullName(self, direction, mode, field=None):
        if field:
            return "%s_%s_%s_%s" % (self.name, direction, field, mode)
        else:
            return "%s_%s_%s" % (self.name, direction, mode)

start_of_file = ""


def read_and_parse_stencils():
    # Read the stencils
    stencils = []
    with open("stencils_cleaned.cxx", "r") as f:
        inFunc = 0
        curr_func = []
        for line in f:
            if line[:5] == 'const':
                global start_of_file
                start_of_file += line
            line = line.strip()
            # avoid to overwrite result
            line = line.replace("result", "result_")
            if line == '':
                continue
            # start of function
            inFunc += line.count("{")
            if inFunc:
                curr_func.append(line)
            inFunc -= line.count("}")
            if inFunc == 0 and curr_func:
                sten = Stencil(curr_func)
                if sten.valid:
                    stencils.append(sten)
                curr_func = []
            ASSERT(inFunc >= 0)

    # Do some cleaning of the stencils bodies
    for sten in stencils:
        ls = ""
        for l in sten.body:
            ls += l
        inv = ls[::-1].find("}") + 1
        ls = ls[ls.find("{") + 1:-inv]
        fu = ls.split(";")
        mymax = len(fu)
        i = 0
        while i < mymax:
            fu[i] += ';'
            if fu[i].find("}") > -1:
                l = fu[i]
                t = l.split("}")
                for ti in t[:-1]:
                    fu.insert(i, ti + "}")
                    i += 1
                    mymax += 1
                fu[i] = t[-1]
            i += 1
        sten.body = fu
    return stencils


def replace_stencil(line, sten, fname, field, mode, sten_mbf, d, update=None, z0=None):
    if update is None:
        update = sten_mbf == "main"
    pos = line.find(sten)
    part_of_offset = ['p', 'm', 'c']
    for i in range(2, 9):
        part_of_offset.append(str(i))
    while pos > -1:
        end = pos + len(sten)
        while line[end] in part_of_offset:
            end += 1
        off = line[pos + 2:end]
        try:
            line_ = line
            line = line[
                :pos] + get_diff(off_diff[mode][off], fname, field, d, update, z0) + line[end:]
        except:
            debug(enable=True)
            debug("Encountered error - some infos:")
            debug(line_, mode, off, sten)
            debug(off_diff)
            debug(off_diff[mode])
            raise
        pos = line.find(sten)
    return line


def parse_body(sten, field, mode, d, z0=None):
    if sten.stag == False and mode != 'norm' and sten.name:
        debug(enable=True)
        debug("Found bug - infos:")
        debug(sten.body, 'name: ' + sten.name, sten.stag, field, mode, d)
        assert(mode == 'norm')
    body = ""
    result = ''
    result_ = [''] * 2
    for line in sten.body[:-1]:
        if sten.flow:
            try:
                line = line.replace("vc", "v.c")
                line = replace_stencil(
                    line, 'v.', "v_in", field, mode, sten.mbf, d, z0=z0)
                line = replace_stencil(
                    line, 'f.', "f_in", field, "norm", sten.mbf, d, z0=z0)
            except:
                debug("Some infos on error:", sten.name, sten.mbf, enable=True)
                raise
        else:
            line = replace_stencil(
                line, 'f.', "in", field, mode, sten.mbf, d, z0=z0)
        if line.find("return") == -1:
            if sten.mbf == 'main':
                body += "     " + line + "\n"
            else:
                toPrint = True
                for resi, res in enumerate(["result_.inner", "result_.outer"]):
                    if line.find(res) > -1:
                        tmp = line[line.index(res) + len(res):]
                        if tmp.find("=") > -1:
                            if result_[resi] != '':
                                raise RuntimeError("Did not expect another defintion of", res,
                                                   "The last one was %s = %s" % (
                                                       res, result_[resi]),
                                                   "thise one is ", line)
                            result_[resi] = tmp[tmp.index("=") + 1:]
                            toPrint = False
                if toPrint:
                    if line.find("=") > -1:
                        raise RuntimeError("While parsing function %s" % sten.name, "\n",
                                           sten.body, "\n",
                                           sten.mbf, "\n",
                                           "Failed to parse - unexpected line: ", line, "\n",
                                           result_, "\n",
                                           line)

        else:
            if sten.mbf == 'main':
                result = line[len("return") + line.index("return"):] + "\n"
            else:
                returned = line[len("return") + line.index("return"):]
    return [body, result, result_]


def get_for_loop_z(sten, field, stag):
    d = 'z'
    print('  if (LocalNz > 3) {')
    for d2 in dirs[field]:
        if d != d2:
            print("  for (int %s = 0 ; %s < LocalN%s; ++%s ){" % (d2, d2, d2, d2))
    body, result, _ = parse_body(sten, field, stag, d)
    for i in range(guards_[0]):
        print('    {')
        print('      int z=%d;' % i)
        body, result, _ = parse_body(sten, field, stag, d, z0=i + 1)
        print(body)
        print("      " + get_diff('c()', "result", field, d) + "=", end=' ')
        if result:
            print(result)
        else:
            print('result_;')
        print('    }')
    print('    for (int z=%d;z<LocalNz-%d;++z) {' % (guards_[0], guards_[1]))
    body, result, _ = parse_body(sten, field, stag, d)
    print(body)
    print("      " + get_diff('c()', "result", field, d) + "=", end=' ')
    if result:
        print(result)
    else:
        print('result_;')
    print('    }')
    for i in range(guards_[1], 0, -1):
        print('    {')
        print('      int z=LocalNz-%d;' % i)
        body, result, _ = parse_body(sten, field, stag, d, z0=-i)
        print(body)
        print("      " + get_diff('c()', "result", field, d) + "=", end=' ')
        if result:
            print(result)
        else:
            print('result_;')
        print('    }')
    for d2 in dirs[field]:
        if d != d2:
            print("  }")
    print('  } else {')
    for d2 in dirs[field]:
        if d != d2:
            print("  for (int %s = 0 ; %s < LocalN%s; ++%s ){" % (d2, d2, d2, d2))
    print('    for (int z=0;z<LocalNz;++z) {')
    body, result, _ = parse_body(sten, field, stag, d, z0='secure')
    print(body)
    print("      " + get_diff('c()', "result", field, d) + "=", end=' ')
    if result:
        print(result)
    else:
        print('result_;')
    print('    }')
    for d2 in dirs[field]:
        if d != d2:
            print("  }")
    print('  }')


def get_for_loop(d, mode, field, guards, sten_name):
    if sten_name == "main":
        print('#if CHECK > 0')
        print('  if (%sstart < %d){' % (d, max(guards)))
        print('    throw BoutException("Cannot compute derivative - need at least %d guard cells in %s direction!");' %
              (max(guards), d.upper()))
        print('  }')
        print('#endif')
        dp = guards[0]
        dm = -guards[1]
        for d2 in dirs[field]:
            if d == d2:
                print(
                    "  for (int %s = %d; %s < LocalN%s%+d; ++%s ){" % (d, dp, d, d, dm, d))
            else:
                print(
                    "  for (int %s = 0 ; %s < LocalN%s; ++%s ){" % (d2, d2, d2, d2))
    else:
        if sten_name == 'forward':
            print("    int " + d, "=%d ;" % (guards_[0] - 1))
        else:
            print("    " + d, "=LocalN%s" % d, "-%d ;" % (guards_[1]))
            if guards_[1] > 2:
                raise RuntimeError(guards_, "To many guards")
        for d2 in perp_dir[field][d]:
            print("    for (int " + d2, "=0; " +
                  d2, "< LocalN" + d2, ";++" + d2, ") {")


def get_for_end(d, field, sten_name):
    if sten_name == 'main':
        for d2 in dirs[field]:
            print("  }")
    else:
        for d2 in perp_dir[field][d]:
            print("  }")


def get_diff(diff, fname, field, d, update=False, z0=None):
    global use_field_operator
    if use_field_operator:
        ret = fname + "("
    else:
        ret = fname + "_ptr["
    for d2 in dirs[field][1:]:
        if not use_field_operator:
            ret += "("
    first = True
    for d2 in dirs[field]:
        if not first:
            if use_field_operator:
                ret += ','
            else:
                ret += ')*LocalN%s + ' % d2
        first = False
        if (d2 == d):
            # do ugly stuff
            global guards_
            diffd = diff2[diff]
            if update:
                if abs(diffd) > 2:
                    raise RuntimeError(
                        "We do not expect to have more then 2 guard cells")
                if diffd > guards_[1]:
                    guards_[1] = diffd
                elif -diffd > guards_[0]:
                    guards_[0] = -diffd
                # else:
                    # print diffd,guards_;
            if d == 'z' and z0 is not None:
                # We want to generate code for the interp_to case with wrapping
                if z0 == 'secure':
                    ret += "+((" + d + "%+d" % (diffd)
                    if diffd < 0:
                        ret += '+%d*LocalNz' % (-diffd)
                    ret += ')%LocalNz)'
                elif z0 > 0 and -diffd >= z0:
                    ret += d2 + "%+d+LocalNz" % (diffd)
                elif z0 < 0 and -diffd <= z0:
                    ret += d2 + "%+d-LocalNz" % (diffd)
                else:
                    ret += d2 + "%+d" % (diffd)
            else:
                ret += d2 + "%+d" % (diffd)
        else:
            ret += d2
    if use_field_operator:
        ret += ")\n    \t\t\t\t"
    else:
        ret += "]\n    \t\t\t\t"
    return ret


def get_pointer(field, field_type, const):
    if const:
        print('  checkData(%s);' % field)
        print("  const BoutReal * __restrict__", end=' ')
    else:
        print("  BoutReal * __restrict__", end=' ')
    print("%s_ptr = &" % field, field, "(", end=' ')
    first = True
    for d in dirs[field_type]:
        if not first:
            print(',', end=' ')
        first = False
        print("0", end=' ')
    print(");")

stencils = read_and_parse_stencils()


def print_stencil_implementations(header_only=False):
    print(start_of_file)
    for stencil in stencils:
        ASSERT(stencil.valid)
        flow = stencil.flow
        if stencil.stag:
            modes = ['CtoL', 'LtoC']
        else:
            modes = ['norm']
        for mode in modes:
            for field in fields:
                for d in dirs[field]:
                    if not header_only:
                        text="void AiolosMesh::"
                    else:
                        text="void "
                    text += stencil.getFullName(direction=d, field=field, mode=mode)
                    text += "(BoutReal * __restrict__ result_ptr"
                    if flow:
                        text += ", const BoutReal * __restrict__ v_in_ptr"
                        text += ", const BoutReal * __restrict__ f_in_ptr"
                    else:
                        text += ", const BoutReal * __restrict__ in_ptr"
                    text+=") const"
                    if header_only:
                        print(text+";")
                    else:
                        print(text+"{")
                        global guards_
                        guards_ = [0, 0]
                        if d == 'z':
                            get_for_loop_z(stencil, field, mode)
                        else:
                            body, result, result_ = parse_body(
                                stencil, field, mode, d)
                            get_for_loop(d, mode, field, guards_, "main")
                            if body + result == '':
                                raise RuntimeError(stencil.body)
                            print(body)
                            print("      " + get_diff('c()',
                                                      "result", field, d) + "=", end=' ')
                            if result:
                                print(result)
                            else:
                                print('result_;')
                            get_for_end(d, field, "main")

                        print("}")
                    warn()
                    text=field+" "
                    if not header_only:
                        text+="AiolosMesh::"
                    text += stencil.getFullName(direction=d, mode=mode)
                    text += "(const " + field
                    if flow:
                        text += " &v_in, const " + field + " &f_in) const"
                    else:
                        text += " &in) const"
                    if header_only:
                        print(text+";")
                    else:
                        print(text+"{")
                        if flow:
                            print("  ASSERT1(this == v_in.getMesh());")
                            print("  ASSERT1(this == f_in.getMesh());")
                        else:
                            print("  ASSERT1(this == in.getMesh());")
                        print('  output_debug.write("Using method %s!\\n");' % stencil.getFullName(
                            direction=d, field=field, mode=mode))
                        print('#if CHECK > 0')
                        print(
                            'if (LocalN%s < %d) {' % (d, sum(guards_) + 1))
                        print('  throw BoutException("%s AiolosMesh::%s'
                              ' - Not enough guards cells to take derivative!");' % (
                                  field, stencil.getFullName(direction=d, mode=mode)))
                        print('}')
                        print('#endif')
                        print(" ", field, "result((AiolosMesh*)this);")
                        print("  result.allocate();")
                        get_pointer("result", field, False)
                        if flow:
                            get_pointer("v_in", field, True)
                            get_pointer("f_in", field, True)
                        else:
                            get_pointer("in", field, True)
                        print("  " + stencil.getFullName(direction=d, field=field, mode=mode) +
                              "(result_ptr,")
                        if flow:
                            print("v_in_ptr,f_in_ptr")
                        else:
                            print("in_ptr")
                        print(");")
                        if mode == "CtoL":
                            print("  result.setLocation(CELL_%sLOW);" %
                                  d.upper())
                        elif mode == "LtoC":
                            print("  result.setLocation(CELL_CENTRE);")
                        else:
                            if flow:
                                print("  result.setLocation(f_in.getLocation());")
                            else:
                                print("  result.setLocation(in.getLocation());")
                        print("""  checkData(result);
                                   return result;
                                 }
                        """)

use_field_operator = False


def get_interp_vals(order, pos):
    import numpy as np
    rhs = np.zeros(order)
    rhs[0] = 1
    mat = np.zeros((order, order))
    x0 = -pos + (order - 1.) / 2.
    for i in range(order):
        x = x0 - i
        for j in range(order):
            mat[j, i] = x**j * np.math.factorial(j)
    # debug(mat)
    facs = np.dot(np.linalg.inv(mat), rhs)
    return facs

useFloat = False
staticCastFloat = False


def float_to_string_with_sign(val):
    global useFloat
    if useFloat:
        return "%+.5e" % val
    else:
        from fractions import Fraction
        val = Fraction(val).limit_denominator()
        ret = ""
        if val > 0:
            ret += "+"
        else:
            val *= -1
            ret += "-"
        if staticCastFloat:
            ret += "static_cast<const BoutReal>"
        ret += "("
        ret += "%d." % val.numerator
        ret += "/%d." % val.denominator
        ret += ")"
        # ret+="std::ratio<"
        # ret+="%d,"%val.numerator
        # ret+="%d"%val.denominator
        # ret+=">"
        return ret


def get_interp_sten(order, pos):
    if order % 2:
        raise RuntimeError("Cannot handle uneven order!")
    oh = order / 2
    if pos:
        oh -= pos / abs(pos)
    vals = get_interp_vals(order, pos)
    ret = ""
    first = True
    for i in range(order):
        ret += float_to_string_with_sign(vals[i])
        ret += "*"
        # debug(ret)
        if i < oh:
            ret += 'f.m' + (str(int(oh - i)) if oh - i > 1 else "")
        else:
            ret += 'f.p' + (str(int(i - oh + 1)) if i - oh + 1 > 1 else "")
    # debug(vals)
    # debug(ret)
    return ret + " ;"


def print_interp_to_code(header_only=False):
    interp = ['', "return " + get_interp_sten(4, 0), '']
    # Hard coded version:
    #["return ( 9.*(s.m + s.p) - s.mm - s.pp ) / 16.;"]
    field = fields[0]  # 3d
    for mode in ['CtoL', 'LtoC']:
        for order in [4]:
            for d in dirs[field]:
                global guards_
                guards_ = [0, 0]
                line = interp[1]
                sten_name = "main"
                line = replace_stencil(
                    line, 'f.', "in", field, mode, sten_name, d)
                text="void "
                if not header_only:
                    text+="AiolosMesh::"
                text+="interp_to_%s_%s_%s(" % (mode, field, d)
                if use_field_operator:
                    text+=field + "& result, const " + field + " & in"
                else:
                    text+="BoutReal * __restrict__ result_ptr, const BoutReal * __restrict__ in_ptr"
                text+=") const"
                print(text)
                if header_only:
                    print(";")
                    continue
                print("{")
                if d == 'z':
                    sten = Stencil(['', interp[1], ''])
                    get_for_loop_z(sten, field, mode)
                else:
                    body = "    " + get_diff('c()', "result", field, d, update=True) + "= " + line[
                        len("return") + line.index("return"):] + "\n"
                    get_for_loop(d, mode, field, guards_, sten_name)
                    print(body)
                    guards__ = guards_
                    get_for_end(d, field, sten_name)
                    sten_names = ["forward", "backward"]
                    for sten_name in sten_names:
                        sten_name_index = sten_names.index(sten_name)
                        sign = -1
                        if sten_name == "forward":
                            sign = 1
                        _sign = sign
                        if d == 'z':
                            sign = 0
                        get_for_loop(d, mode, field, guards_, sten_name)
                        print("      " + get_diff('c()',
                                                  "result", field, d) + "=", end=' ')
                        print(replace_stencil(get_interp_sten(4, sign), 'f.', "in", field,
                                              mode, sten_name, d, False, z0=((guards_[sten_name_index]) * _sign)))
                        guards_ = guards__
                        if order / 2 > 1:
                            if (sten_name == 'backward' and mode == 'CtoL') or \
                               (sten_name == 'forward' and mode == 'LtoC'):
                                pass  # dont do anything ...
                            else:
                                if sten_name == 'forward':
                                    print("        " + get_diff('m()',
                                                                "result", field, d) + "=", end=' ')
                                else:
                                    print("        " + get_diff('p()',
                                                                "result", field, d) + "=", end=' ')
                                print(replace_stencil(get_interp_sten(4, sign * 2), 'f.', "in", field,
                                                      mode, sten_name, d, False, z0=((guards_[sten_name_index]) * _sign)))
                                guards_ = guards__
                        get_for_end(d, field, sten_name)
                print("}")
    if header_only:
        return
    print(
        "const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc, REGION region) const {")
    print("  Field3D result((AiolosMesh*)this);")
    print("  result.allocate();")
    if not use_field_operator:
        print("  Indices i0{0,0,0};")
    print("  if (f.getLocation() != CELL_CENTRE){")
    print("    // we first need to go back to centre before we can go anywhere else")
    print("    switch (f.getLocation()){")
    for d in dirs[field]:
        print("    case CELL_%sLOW:" % d.upper())
        if use_field_operator:
            print("      interp_to_LtoC_%s_%s(result,f,this);" % (field, d))
        else:
            print(
                "      interp_to_LtoC_%s_%s(&result[i0],&f[i0]);" % (field, d))
        print("      result.setLocation(CELL_CENTRE);")
        print("      // return or interpolate again")
        print("      return interp_to(result, loc, region);")
        print("      break;")
    print("    default:")
    print('      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",strLocation(loc));')
    print("    }")
    print("  }")
    print("  // we are in CELL_CENTRE and need to go somewhere else ...")
    print("  switch (loc){")
    for d in dirs[field]:
        print("    case CELL_%sLOW:" % d.upper())
        if use_field_operator:
            print("      interp_to_CtoL_%s_%s(result,f);" % (field, d))
        else:
            print(
                "      interp_to_CtoL_%s_%s(&result[i0],&f[i0]);" % (field, d))
        print("      result.setLocation(CELL_%sLOW);" % d.upper())
        print("      // return or interpolate again")
        print("      return interp_to(result, loc, region);")
        print("      break;")
    print("    default:")
    print('      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",strLocation(loc));')
    print("    }")
    print("}")

use_field_operator = False
