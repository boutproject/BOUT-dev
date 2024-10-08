#!/usr/bin/env python3
from dataclasses import dataclass
from pathlib import Path
import argparse
import difflib
import copy
import textwrap
import os
import re


@dataclass
class VersionNumber:
    major: int
    minor: int
    patch: int

    def __str__(self):
        return f"{self.major}.{self.minor}.{self.patch}"

    @property
    def short(self) -> str:
        return f"{self.major}.{self.minor}"


def get_full_filepath(filepath):
    main_directory = Path(os.path.abspath(__file__)).parent.parent
    return Path(main_directory) / filepath


def update_version_number_in_file(relative_filepath, pattern, new_version_number):
    full_filepath = get_full_filepath(relative_filepath)

    with open(full_filepath, "r", encoding="UTF-8") as file:
        file_contents = file.read()
        original = copy.deepcopy(file_contents)

        modified = apply_fixes(pattern, new_version_number, file_contents)
        patch = create_patch(str(full_filepath), original, modified)

        if args.patch_only:
            print(patch)
            return

        if not patch:
            if not args.quiet:
                print("No changes to make to {}".format(full_filepath))
            return

        if not args.quiet:
            print("\n******************************************")
            print("Changes to {}\n".format(full_filepath))
            print(patch)
            print("\n******************************************")

        if args.force:
            make_change = True
        else:
            make_change = yes_or_no("Make changes to {}?".format(full_filepath))

        if make_change:
            with open(full_filepath, "w", encoding="UTF-8") as file:
                file.write(modified)


def bump_version_numbers(
    new_version_number: VersionNumber, next_version_number: VersionNumber
):
    update_version_number_in_file(
        "configure.ac",
        r"^AC_INIT\(\[BOUT\+\+\],\[(\d+\.\d+\.\d+)\]",
        new_version_number,
    )
    update_version_number_in_file(
        "CITATION.cff", r"^version: (\d+\.\d+\.\d+)", new_version_number
    )
    update_version_number_in_file(
        "manual/sphinx/conf.py", r"^version = \"(\d+\.\d+)\"", new_version_number.short
    )
    update_version_number_in_file(
        "manual/sphinx/conf.py", r"^release = \"(\d+\.\d+\.\d+)\"", new_version_number
    )
    update_version_number_in_file(
        "manual/doxygen/Doxyfile_readthedocs",
        r"^PROJECT_NUMBER\s*=\s*(\d+\.\d+\.\d+)",
        new_version_number,
    )
    update_version_number_in_file(
        "manual/doxygen/Doxyfile",
        r"^PROJECT_NUMBER\s*=\s*(\d+\.\d+\.\d+)",
        new_version_number,
    )
    update_version_number_in_file(
        "CMakeLists.txt",
        r"^set\(_bout_previous_version \"(\d+\.\d+\.\d+)\"\)",
        new_version_number,
    )
    update_version_number_in_file(
        "CMakeLists.txt",
        r"^set\(_bout_next_version \"(\d+\.\d+\.\d+)\"\)",
        next_version_number,
    )
    update_version_number_in_file(
        "tools/pylib/_boutpp_build/backend.py",
        r"_bout_previous_version = \"v(\d+\.\d+\.\d+)\"",
        new_version_number,
    )
    update_version_number_in_file(
        "tools/pylib/_boutpp_build/backend.py",
        r"_bout_next_version = \"v(\d+\.\d+\.\d+)\"",
        next_version_number,
    )


def apply_fixes(pattern, new_version_number, source):
    """Apply the various fixes for each factory to source. Returns
    modified source

    Parameters
    ----------
    pattern
        Regex pattern to apply for replacement
    new_version_number
        New version number to use in replacement
    source
        Text to update
    """

    def get_replacement(match):
        return match[0].replace(match[1], str(new_version_number))

    modified = re.sub(pattern, get_replacement, source, flags=re.MULTILINE)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Update the software version number to the specified version, 
            to be given in the form major.minor.patch, 
            e.g. 5.10.3
            
            Where the 3rd ('patch') part of the version is omitted, 
            only the 'major' and 'minor' parts will be used, 
            e.g. 5.10.3 -> 5.10 
            
            For the 'bout-next' version number, 
            the 'minor' version number of the provided version will be incremented by 1, 
            e.g. 5.10.3 -> 5.11.3
            
            """
        ),
    )

    parser.add_argument(
        "--force", "-f", action="store_true", help="Make changes without asking"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Don't print patches"
    )
    parser.add_argument(
        "--patch-only", "-p", action="store_true", help="Print the patches and exit"
    )
    parser.add_argument("new_version", help="New (current) version number")
    parser.add_argument(
        "next_version", help="Next version number", nargs="?", default=None
    )

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    new_version = VersionNumber(*map(int, args.new_version.split(".")))
    next_version = (
        VersionNumber(new_version.major, new_version.minor + 1, 0)
        if args.next_version is None
        else VersionNumber(*map(int, args.next_version.split(".")))
    )

    bump_version_numbers(
        new_version_number=new_version, next_version_number=next_version
    )
