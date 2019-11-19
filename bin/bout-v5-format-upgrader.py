#!/usr/bin/env python3

import argparse
import copy
import difflib
import re


format_replacements = {
    "%c": "{:c}",
    "%d": "{:d}",
    "%e": "{:e}",
    "%f": "{:f}",
    "%g": "{:g}",
    "%i": "{:d}",
    "%ld": "{:d}",
    "%le": "{:e}",
    "%lu": "{:d}",
    "%p": "{:p}",
    "%s": "{:s}",
    "%zu": "{:d}",
}


def fix_format_replacement(format_replacement, source):
    """Replace printf format with fmt format

    """
    return re.sub(format_replacement[0], format_replacement[1], source)


def fix_trivial_format(source):
    """Reduce trivial formatting of strings to just the string

    """

    def trivial_replace(match):
        if match.group(2):
            return "{}{}{}".format(match.group(1), match.group(2), match.group(4))
        if match.group(3):
            return "{}{}{}".format(match.group(1), match.group(3), match.group(4))
        raise ValueError("Found an unexpected match: {}".format(match))

    return re.sub(
        r"""
        (.*)?
        "{:s}",\s*                  # Entire format is just a string
        (?:([\w_]+)\.c_str\(\)      # And replacement is std::string::c_str
        |(".*?"))
        (.*)?
        """,
        trivial_replace,
        source,
        flags=re.VERBOSE,
    )


def fix_string_c_str(source):
    """Fix formats that use {:s} where the replacement is using std::string::c_str

    """
    return re.sub(
        r"""
        (".*{:s}[^;]*?",)       # A format string containing {:s}
        \s*([^);]+?)\.c_str\(\) # Replacement of std::string::c_str
        """,
        r"\1 \2",
        source,
        flags=re.DOTALL | re.VERBOSE,
    )


def fix_trace(source):
    """Fix TRACE macros where fix_string_c_str has failed for some reason

    """
    return re.sub(
        r"""
        (TRACE\(".*{:s}.*",)
        \s*([\w_]+)\.c_str\(\)\); # Replacement of std::string::c_str
        """,
        r"\1 \2);",
        source,
        flags=re.VERBOSE,
    )


def apply_fixes(format_replacements, source):
    """Apply the various fixes for each factory to source. Returns
    modified source

    Parameters
    ----------
    factories
        Dictionary of factory properties
    source
        Text to update
    """

    modified = source

    for format_replacement in format_replacements.items():
        modified = fix_format_replacement(format_replacement, modified)

    modified = fix_trivial_format(modified)
    modified = fix_string_c_str(modified)
    modified = fix_trace(modified)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix format specifiers")

    parser.add_argument("files", action="store", nargs="+", help="Input files")
    parser.add_argument(
        "--force", "-f", action="store_true", help="Make changes without asking"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Don't print patches"
    )
    parser.add_argument(
        "--patch-only", "-p", action="store_true", help="Print the patches and exit"
    )

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()
        original = copy.deepcopy(contents)

        modified = apply_fixes(format_replacements, contents)
        patch = create_patch(filename, original, modified)

        if args.patch_only:
            print(patch)
            continue

        if not patch:
            if not args.quiet:
                print("No changes to make to {}".format(filename))
            continue

        if not args.quiet:
            print("\n******************************************")
            print("Changes to {}\n".format(filename))
            print(patch)
            print("\n******************************************")

        if args.force:
            make_change = True
        else:
            make_change = yes_or_no("Make changes to {}?".format(filename))

        if make_change:
            with open(filename, "w") as f:
                f.write(modified)
