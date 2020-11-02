#!/usr/bin/env python3

import argparse
import copy
import difflib
import textwrap
import warnings

from boutdata.data import BoutOptionsFile
from boututils.boutwarnings import AlwaysWarning

# This should be a list of dicts, each containing "old", "new" and optionally "values".
# The values of "old"/"new" keys should be the old/new names of input file values or
# sections. The value of "values" is a dict containing replacements for values of the
# option.
REPLACEMENTS = [
    {"old": "mesh:paralleltransform", "new": "mesh:paralleltransform:type"},
    {"old": "fci", "new": "mesh:paralleltransform"},
    {"old": "interpolation", "new": "mesh:paralleltransform:xzinterpolation"},
    {"old": "fft:fft_measure", "new": "fft:fft_measurement_flag",
        "values": {"false": "estimate", "true": "measure"}}
]


def fix_replacements(replacements, options_file):
    """Change the names of options in options_file according to the list
    of dicts replacements

    """
    for replacement in replacements:
        try:
            options_file.rename(replacement["old"], replacement["new"])
        except KeyError:
            pass
        except TypeError as e:
            raise RuntimeError(
                "Could not apply transformation: '{old}' -> '{new}' to file '{0}', due to error:"
                "\n\t{1}".format(options_file.filename, e.args[0], **replacement)
            ) from e
        else:
            if "values" in replacement:
                old_value = options_file[replacement["new"]]
                options_file[replacement["new"]] = replacement["values"].get(
                    old_value, old_value
                )


def apply_fixes(replacements, options_file):
    """Apply all fixes in this module
    """

    modified = copy.deepcopy(options_file)

    fix_replacements(replacements, modified)

    return modified


def yes_or_no(question):
    """Convert user input from yes/no variations to True/False

    """
    while True:
        reply = input(question + " [y/N] ").lower().strip()
        if not reply or reply[0] == "n":
            return False
        if reply[0] == "y":
            return True


def create_patch(filename, original, modified):
    """Create a unified diff between original and modified
    """

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
