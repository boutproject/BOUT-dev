#!/usr/bin/env python3

import boutupgrader

import argparse
import copy
import re
import textwrap
from pathlib import Path
from typing import List
from subprocess import run


header_shim_sentinel = "// BOUT++ header shim"

header_warning = f"""\
#pragma once
{header_shim_sentinel}
#warning Header "{{0}}" has moved to "bout/{{0}}". Run `bin/bout-v5-header-upgrader.py` to fix
#include "bout/{{0}}"
"""


def header_needs_moving(header: Path) -> bool:
    """Check if `header` has not yet been moved"""
    with open(header, "r") as f:
        return header_shim_sentinel not in f.read()


def deprecated_header_list(include_path: Path = Path("./include")):
    """List of deprecated header paths (that is, those in bare
    ``include`` directory)

    """
    return include_path.glob("*.hxx")


def write_header_shim(header: Path):
    """Write 'shim' for header, that ``include``s new location"""
    with open(header, "w") as f:
        f.write(header_warning.format(header.name))


def fix_library_header_locations(
    include_path: Path = Path("./include"), quiet: bool = False
):
    unmoved_headers = list(
        filter(header_needs_moving, deprecated_header_list(include_path))
    )
    include_bout_path = include_path / "bout"

    if unmoved_headers == []:
        print("No headers to move!")
        return

    out = run("git diff-index --cached HEAD --quiet", shell=True)
    if out.returncode:
        raise RuntimeError(
            "git index not clean! Please commit or discard any changes before continuing"
        )

    # First we move any existing headers and commit this change, so
    # that history is preserved
    for header in unmoved_headers:
        new_path = include_bout_path / header.name
        if not quiet:
            print(f"Moving '{header}' to '{new_path}'")
        run(f"git mv {header} {new_path}", shell=True, check=True)

    run(r"git commit -m 'Move headers to `include/bout`'", shell=True, check=True)

    # Now we can write the compatibility shim
    for header in unmoved_headers:
        write_header_shim(header)

    run(f"git add {' '.join(map(str, unmoved_headers))}", shell=True, check=True)
    run(
        r"git commit -m 'Add compatibility shims for old header locations'",
        shell=True,
        check=True,
    )


def make_header_regex(deprecated_headers: List[str]) -> re.Pattern:
    """Create a regular expression to match deprecated header locations"""
    deprecated_header_joined = "|".join((header.name for header in deprecated_headers))
    return re.compile(rf'(#include\s+<|")(?:\.\./)?({deprecated_header_joined})(>|")')


def apply_fixes(header_regex, source):
    """Apply all fixes in this module"""
    modified = copy.deepcopy(source)

    return header_regex.sub(r"\1bout/\2\3", modified)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Fix deprecated header locations for BOUT++ v4 -> v5

            All BOUT++ headers are now under ``include/bout`` and
            should be included as ``#include <bout/header.hxx>``. This
            tool will fix such includes.

            For developers: the option ``--move-deprecated-headers``
            will move the headers from ``include`` to
            ``include/bout``, and add a compatibility shim in the old
            location. This option is mutually exclusive with
            ``--files``, and should be used after running this tool
            over the library files.

            WARNING: If any files do need moving, this will create a
            git commit in order to preserve history across file moves.
            If you have staged changes, this tool will not work, so to
            avoid committing undesired or unrelated changes.

            """
        ),
    )
    parser = boutupgrader.default_args(parser, pos_files=False)

    parser.add_argument(
        "--include-path",
        "-i",
        help="Path to `include` directory",
        default="./include",
        type=Path,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--files", nargs="*", action="store", help="Input files")
    group.add_argument(
        "--move-deprecated-headers",
        action="store_true",
        help="Move the deprecated headers",
    )

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    deprecated_headers = deprecated_header_list(args.include_path)

    if args.move_deprecated_headers:
        fix_library_header_locations(args.include_path, args.quiet)
        exit(0)

    header_regex = make_header_regex(deprecated_headers)

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()
        original = copy.deepcopy(contents)

        modified = apply_fixes(header_regex, contents)

        boutupgrader.apply_or_display_patch(
            filename, original, modified, args.patch_only, args.quiet, args.force
        )
