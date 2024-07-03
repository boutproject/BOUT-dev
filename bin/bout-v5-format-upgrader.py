#!/usr/bin/env python3

import boutupgrader

import argparse
import copy
import re


format_replacements = {
    "c": "c",
    "d": "d",
    "e": "e",
    "f": "f",
    "g": "g",
    "i": "d",
    "ld": "d",
    "le": "e",
    "lu": "d",
    "p": "p",
    "s": "s",
    "zu": "d",
}


def fix_format_replacement(format_replacement, source):
    """Replace printf format with fmt format"""
    return re.sub(
        r"%([0-9]*\.?[0-9]*){}".format(format_replacement[0]),
        r"{{:\1{}}}".format(format_replacement[1]),
        source,
    )


def fix_trivial_format(source):
    """Reduce trivial formatting of strings to just the string"""

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
    """Fix formats that use {:s} where the replacement is using std::string::c_str"""
    return re.sub(
        r"""
        (".*{:s}[^;]*?",)       # A format string containing {:s}
        \s*([^;]+?)\.c_str\(\) # Replacement of std::string::c_str
        """,
        r"\1 \2",
        source,
        flags=re.VERBOSE,
    )


def fix_trace(source):
    """Fix TRACE macros where fix_string_c_str has failed for some reason"""
    return re.sub(
        r"""
        (TRACE\(".*{:s}.*",)
        \s*([\w_]+)\.c_str\(\)\); # Replacement of std::string::c_str
        """,
        r"\1 \2);",
        source,
        flags=re.VERBOSE,
    )


def fix_toString_c_str(source):
    """Fix formats that call toString where the replacement is using std::string::c_str"""
    return re.sub(
        r"""
        (".*{:s}[^;]*?",.*?)         # A format string containing {:s}
        (toString\(.*?\))\.c_str\(\) # Replacement of std::string::c_str
        """,
        r"\1\2",
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


if __name__ == "__main__":
    parser = boutupgrader.default_args(
        argparse.ArgumentParser(description="Fix format specifiers")
    )
    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()
        original = copy.deepcopy(contents)

        modified = apply_fixes(format_replacements, contents)

        boutupgrader.apply_or_display_patch(
            filename, original, modified, args.patch_only, args.quiet, args.force
        )
