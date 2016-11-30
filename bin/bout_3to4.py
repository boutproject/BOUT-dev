#!/usr/bin/env python3

import argparse
import re
import fileinput

nonmembers = {
    'DC': ['DC', 1],
    'slice': ['sliceXZ', 2],
}

coordinates = [
    "outputVars",
    "dx", "dy", "dz",
    "zlength",
    "non_uniform",
    "d1_dx", "d1_dy",
    "J", "Bxy",
    "g11", "g22", "g33", "g12", "g13", "g23",
    "g_11", "g_22", "g_33", "g_12", "g_13", "g_23",
    "G1_11", "G1_22", "G1_33", "G1_12", "G1_13",
    "G2_11", "G2_22", "G2_33", "G2_12", "G2_23",
    "G3_11", "G3_22", "G3_33", "G3_13", "G3_23",
    "G1", "G2", "G3",
    "ShiftTorsion", "IntShiftTorsion",
    "geometry", "calcCovariant", "calcContravariant", "jacobian"
]


def fix_nonmembers(line_text, filename, line_num, replace=False):
    """Replace member functions with nonmembers
    """

    old_line_text = line_text

    for old, (new, num_args) in nonmembers.items():
        pattern = re.compile("(\w*)\.{}\(".format(old))
        matches = re.findall(pattern, line_text)
        for match in matches:
            replacement = "{func}({object}".format(func=new, object=match)
            if num_args > 1:
                replacement += ", "
            line_text = re.sub(pattern, replacement, line_text)
            if not replace:
                name_num = "{name}:{num}:".format(name=filename, num=line_num)
                print("{name_num}{line}".format(name_num=name_num, line=old_line_text), end='')
                print(" "*len(name_num) + line_text)
    if replace:
        return line_text


def fix_subscripts(line_text, filename, line_num, replace=False):
    """Replace triple square brackets with round brackets

    Should also check that the variable is a Field3D/Field2D - but doesn't
    """

    old_line_text = line_text
    # Catch both 2D and 3D arrays
    pattern = re.compile(r"\[([^[]*)\]\[([^[]*)\](?:\[([^[]*)\])?")
    matches = re.findall(pattern, line_text)
    for match in matches:
        # If the last group is non-empty, then it was a 3D array
        if len(match[2]):
            replacement = r"(\1, \2, \3)"
        else:
            replacement = r"(\1, \2)"
        line_text = re.sub(pattern, replacement, line_text)
        if not replace:
            name_num = "{name}:{num}:".format(name=filename, num=line_num)
            print("{name_num}{line}".format(name_num=name_num, line=old_line_text), end='')
            print(" "*len(name_num) + line_text)
    if replace:
        return line_text


def fix_coordinates(line_text, filename, line_num, replace=False):
    """Fix variables that have moved from mesh to coordinates
    """

    old_line_text = line_text

    for var in coordinates:
        pattern = re.compile("mesh->{}".format(var))
        matches = re.findall(pattern, line_text)
        for match in matches:
            line_text = re.sub(pattern, "mesh->coordinates()->{}".format(var), line_text)
            if not replace:
                name_num = "{name}:{num}:".format(name=filename, num=line_num)
                print("{name_num}{line}".format(name_num=name_num, line=old_line_text), end='')
                print(" "*len(name_num) + line_text)
    if replace:
        return line_text


if __name__ == '__main__':

    epilog = """
    Currently bout_3to4 can detect the following transformations are needed:
        - Triple square brackets instead of round brackets for subscripts
        - Field member functions that are now non-members
        - Variables/functions that have moved from Mesh to Coordinates

    Note that in the latter case of transformations, you will still need to manually add
        Coordinates *coords = mesh->coordinates();
    to the correct scopes
    """

    parser = argparse.ArgumentParser(description="A little helper for upgrading from BOUT++ version 3 to version 4",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=epilog)
    parser.add_argument("-r", "--replace", action="store_true",
                        help="Actually make the fix")
    parser.add_argument("files", nargs='+',
                        help="Files to process")

    args = parser.parse_args()

    # Loops over all lines across all files
    for line in fileinput.input(files=args.files, inplace=args.replace):
        filename = fileinput.filename()
        line_num = fileinput.filelineno()

        # Apply the transformations and then update the line if we're doing a replacement
        new_line = fix_nonmembers(line, filename, line_num, args.replace)
        line = new_line if args.replace else line

        new_line = fix_subscripts(line, filename, line_num, args.replace)
        line = new_line if args.replace else line

        new_line = fix_coordinates(line, filename, line_num, args.replace)
        line = new_line if args.replace else line

        # If we're doing a replacement, then we need to print all lines, without a newline
        if args.replace:
            print(line, end='')
