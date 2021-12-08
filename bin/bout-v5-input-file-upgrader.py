#!/usr/bin/env python3

import argparse
import copy
import difflib
import itertools
import textwrap
import warnings

from boutdata.data import BoutOptionsFile, BoutOptions
from boututils.boutwarnings import AlwaysWarning


def case_sensitive_init(self, name="root", parent=None):
    self._sections = dict()
    self._keys = dict()
    self._name = name
    self._parent = parent
    self.comments = dict()
    self.inline_comments = dict()
    self._comment_whitespace = dict()


# Monky-patch BoutOptions to make sure it's case sensitive
BoutOptions.__init__ = case_sensitive_init


# This should be a list of dicts, each containing "old", "new" and optionally "new_values".
# The values of "old"/"new" keys should be the old/new names of input file values or
# sections. The value of "new_values" is a dict containing replacements for values of the
# option. "old_type" optionally specifies the type of the old value of the option; for
# example this is needed for special handling of boolean values.
REPLACEMENTS = [
    {"old": "mesh:paralleltransform", "new": "mesh:paralleltransform:type"},
    {"old": "fci", "new": "mesh:paralleltransform"},
    {"old": "interpolation", "new": "mesh:paralleltransform:xzinterpolation"},
    {
        "old": "fft:fft_measure",
        "new": "fft:fft_measurement_flag",
        "old_type": bool,
        "new_values": {False: "estimate", True: "measure"},
    },
    {"old": "TIMESTEP", "new": "timestep"},
    {"old": "NOUT", "new": "nout"},
    {"old": "ddx", "new": "mesh:ddx"},
    {"old": "ddy", "new": "mesh:ddy"},
    {"old": "ddz", "new": "mesh:ddz"},
    {"old": "laplace:laplace_nonuniform", "new": "laplace:nonuniform"},
    {"old": "mesh:dump_format", "new": "dump_format"},
    {"old": "solver:ATOL", "new": "solver:atol"},
    {"old": "solver:RTOL", "new": "solver:rtol"},
    # This was inconsistent in the library
    {"old": "All", "new": "all"},
    # The following haven't been changed, but are frequently spelt with the wrong case
    {"old": "mxg", "new": "MXG"},
    {"old": "myg", "new": "MYG"},
    {"old": "nxpe", "new": "NXPE"},
    {"old": "nype", "new": "NYPE"},
    {"old": "mesh:NX", "new": "mesh:nx"},
    {"old": "mesh:NY", "new": "mesh:ny"},
    {"old": "mesh:shiftangle", "new": "mesh:ShiftAngle"},
    {"old": "mesh:shiftAngle", "new": "mesh:ShiftAngle"},
    {"old": "mesh:zshift", "new": "mesh:zShift"},
    {"old": "mesh:StaggerGrids", "new": "mesh:staggergrids"},
    {"old": "output:shiftOutput", "new": "output:shiftoutput"},
    {"old": "output:ShiftOutput", "new": "output:shiftoutput"},
    {"old": "output:shiftInput", "new": "output:shiftinput"},
    {"old": "output:ShiftInput", "new": "output:shiftinput"},
    {"old": "output:flushFrequency", "new": "output:flushfrequency"},
    {"old": "output:FlushFrequency", "new": "output:flushfrequency"},
    {"old": "TwistShift", "new": "twistshift"},
    {"old": "zmin", "new": "ZMIN"},
    {"old": "zmax", "new": "ZMAX"},
    {"old": "ZPERIOD", "new": "zperiod"},
    # 'restart' can be either a section or a value, so move all the
    # section:values instead
    {"old": "restart:parallel", "new": "restart_files:parallel"},
    {"old": "restart:flush", "new": "restart_files:flush"},
    {"old": "restart:guards", "new": "restart_files:guards"},
    {"old": "restart:floats", "new": "restart_files:floats"},
    {"old": "restart:openclose", "new": "restart_files:openclose"},
    {"old": "restart:enabled", "new": "restart_files:enabled"},
    {"old": "restart:init_missing", "new": "restart_files:init_missing"},
    {"old": "restart:shiftOutput", "new": "restart_files:shiftOutput"},
    {"old": "restart:shiftInput", "new": "restart_files:shiftInput"},
    {"old": "restart:flushFrequency", "new": "restart_files:flushFrequency"},
]

for section, derivative in itertools.product(
    ["ddx", "ddy", "ddz", "diff"], ["First", "Second", "Fourth", "Flux", "Upwind"]
):
    REPLACEMENTS.append(
        {
            "old": f"mesh:{section}:{derivative}",
            "new": f"mesh:{section}:{derivative.lower()}",
        }
    )


def parse_bool(bool_expression):
    try:
        bool_expression_lower = bool_expression.lower()
    except AttributeError:
        # bool_expression was not a str: no need to lower
        bool_expression_lower = bool_expression

    if bool_expression_lower in ["true", "y", "t", 1, True]:
        return True
    elif bool_expression_lower in ["false", "n", "f", 0, False]:
        return False
    else:
        raise RuntimeError(
            f"Expected boolean option. Could not parse {bool_expression}"
        )


def already_fixed(replacement, options_file):
    """Check if the options_file already has already had this particular fix applied"""
    # The old key is there and the new one isn't, then it's definitely not fixed
    if replacement["old"] in options_file and replacement["new"] not in options_file:
        return False
    # If the new isn't there, there's nothing to fix
    if replacement["new"] not in options_file:
        return True
    # If we don't need to fix values, we're done
    if "new_values" not in replacement:
        return True
    # Check if the current value is acceptable
    return options_file[replacement["new"]] in replacement["new_values"].values()


def fix_replacements(replacements, options_file):
    """Change the names of options in options_file according to the list
    of dicts replacements

    """
    for replacement in replacements:
        try:
            if already_fixed(replacement, options_file):
                continue
            options_file.rename(replacement["old"], replacement["new"])
        except KeyError:
            pass
        except TypeError as e:
            raise RuntimeError(
                "Could not apply transformation: '{old}' -> '{new}' to file '{0}', due to error:"
                "\n\t{1}".format(options_file.filename, e.args[0], **replacement)
            ) from e
        else:
            if "old_type" in replacement:
                # Special handling for certain types, replicating what BOUT++ does
                if replacement["old_type"] is bool:
                    # The original value must be something that BOUT++ recognises as a
                    # bool.
                    # replacement["new_values"] must contain both True and False keys.
                    old_value = parse_bool(options_file[replacement["new"]])
                    options_file[replacement["new"]] = replacement["new_values"][
                        old_value
                    ]
                else:
                    raise ValueError(
                        f"Error in REPLACEMENTS: type {replacement['type']} is not handled"
                    )
            else:
                # Option values are just a string
                if "new_values" in replacement:
                    old_value = options_file[replacement["new"]]
                    try:
                        old_value = old_value.lower()
                    except AttributeError:
                        # old_value was not a str, so no need to convert to lower case
                        pass

                    try:
                        options_file[replacement["new"]] = replacement["new_values"][
                            old_value
                        ]
                    except KeyError:
                        # No replacement given for this value: keep the old one
                        pass


def apply_fixes(replacements, options_file):
    """Apply all fixes in this module"""

    modified = copy.deepcopy(options_file)

    fix_replacements(replacements, modified)

    return modified


def yes_or_no(question):
    """Convert user input from yes/no variations to True/False"""
    while True:
        reply = input(question + " [y/N] ").lower().strip()
        if not reply or reply[0] == "n":
            return False
        if reply[0] == "y":
            return True


def create_patch(filename, original, modified):
    """Create a unified diff between original and modified"""

    patch = "\n".join(
        difflib.unified_diff(
            original.splitlines(),
            modified.splitlines(),
            fromfile=filename,
            tofile=filename,
            lineterm="",
        )
    )

    return patch


def possibly_apply_patch(patch, options_file, quiet=False, force=False):
    """Possibly apply patch to options_file. If force is True, applies the
    patch without asking, overwriting any existing file. Otherwise,
    ask for confirmation from stdin

    """
    if not quiet:
        print("\n******************************************")
        print("Changes to {}\n".format(options_file.filename))
        print(patch)
        print("\n******************************************")

    if force:
        make_change = True
    else:
        make_change = yes_or_no("Make changes to {}?".format(options_file.filename))
    if make_change:
        options_file.write(overwrite=True)
    return make_change


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Fix input files for BOUT++ v5+

            Please note that this will only fix input options in sections with
            standard or default names. You may also need to fix options in custom
            sections.

            Warning! Even with no fixes, there may still be changes as this script
            will "canonicalise" the input files:

              * nested sections are moved to be under their parent section, while
                preserving relative order

              * empty sections are removed

              * floating point numbers may have their format changed, although the
                value will not change

              * consecutive blank lines will be reduced to a single blank line

              * whitespace around equals signs will be changed to exactly one space

              * trailing whitespace will be removed

              * comments will always use '#'

            Files that change in this way will have the "canonicalisation" patch
            presented first. If you choose not to apply this patch, the "upgrade
            fixer" patch will still include it."""
        ),
    )

    parser.add_argument("files", action="store", nargs="+", help="Input files")

    force_patch_group = parser.add_mutually_exclusive_group()
    force_patch_group.add_argument(
        "--force", "-f", action="store_true", help="Make changes without asking"
    )
    force_patch_group.add_argument(
        "--patch-only", "-p", action="store_true", help="Print the patches and exit"
    )

    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Don't print patches"
    )
    parser.add_argument(
        "--accept-canonical",
        "-c",
        action="store_true",
        help="Automatically accept the canonical patch",
    )
    parser.add_argument(
        "--canonical-only",
        "-k",
        action="store_true",
        help="Only check/fix canonicalisation",
    )

    args = parser.parse_args()

    warnings.simplefilter("ignore", AlwaysWarning)

    for filename in args.files:
        with open(filename, "r") as f:
            original_source = f.read()

        try:
            original = BoutOptionsFile(filename)
        except ValueError:
            pass

        canonicalised_patch = create_patch(filename, original_source, str(original))
        if canonicalised_patch and not args.patch_only:
            print(
                "WARNING: original input file '{}' not in canonical form!".format(
                    filename
                )
            )
            applied_patch = possibly_apply_patch(
                canonicalised_patch,
                original,
                args.quiet,
                args.force or args.accept_canonical,
            )
            # Re-read input file
            if applied_patch:
                original_source = str(original)

        if args.canonical_only:
            continue

        try:
            modified = apply_fixes(REPLACEMENTS, original)
        except RuntimeError as e:
            print(e)
            continue
        patch = create_patch(filename, original_source, str(modified))

        if args.patch_only:
            print(patch)
            continue

        if not patch:
            if not args.quiet:
                print("No changes to make to {}".format(filename))
            continue

        possibly_apply_patch(patch, modified, args.quiet, args.force)
