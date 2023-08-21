from pathlib import Path
import argparse
import difflib
import copy
import textwrap
import os
import re


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


def bump_version_numbers(new_version_number):
    short_version_number = ShortVersionNumber(
        new_version_number.major_version, new_version_number.minor_version
    )
    bout_next_version_number = VersionNumber(
        new_version_number.major_version,
        new_version_number.minor_version + 1,
        new_version_number.patch_version,
    )

    update_version_number_in_file(
        "configure.ac",
        r"^AC_INIT\(\[BOUT\+\+\],\[(\d+\.\d+\.\d+)\]",
        new_version_number,
    )
    update_version_number_in_file(
        "CITATION.cff", r"^version: (\d+\.\d+\.\d+)", new_version_number
    )
    update_version_number_in_file(
        "manual/sphinx/conf.py", r"^version = \"(\d+\.\d+)\"", short_version_number
    )
    update_version_number_in_file(
        "manual/sphinx/conf.py", r"^release = \"(\d+\.\d+\.\d+)\"", new_version_number
    )
    update_version_number_in_file(
        "manual/doxygen/Doxyfile_readthedocs",
        r"^PROJECT_NUMBER         = (\d+\.\d+\.\d+)",
        new_version_number,
    )
    update_version_number_in_file(
        "manual/doxygen/Doxyfile",
        r"^PROJECT_NUMBER         = (\d+\.\d+\.\d+)",
        new_version_number,
    )
    update_version_number_in_file(
        "CMakeLists.txt",
        r"^set\(_bout_previous_version \"v(\d+\.\d+\.\d+)\"\)",
        new_version_number,
    )
    update_version_number_in_file(
        "CMakeLists.txt",
        r"^set\(_bout_next_version \"(\d+\.\d+\.\d+)\"\)",
        bout_next_version_number,
    )
    update_version_number_in_file(
        "tools/pylib/_boutpp_build/backend.py",
        r"_bout_previous_version = \"v(\d+\.\d+\.\d+)\"",
        new_version_number,
    )
    update_version_number_in_file(
        "tools/pylib/_boutpp_build/backend.py",
        r"_bout_next_version = \"v(\d+\.\d+\.\d+)\"",
        bout_next_version_number,
    )


class VersionNumber:
    major_version: int
    minor_version: int
    patch_version: int

    def __init__(self, major_version, minor_version, patch_version):
        self.major_version = major_version
        self.minor_version = minor_version
        self.patch_version = patch_version

    def as_string(self):
        return "%d.%d.%d" % (self.major_version, self.minor_version, self.patch_version)


class ShortVersionNumber:
    def __init__(self, major_version, minor_version):
        self.major_version = major_version
        self.minor_version = minor_version

    def as_string(self):
        return "%d.%d" % (self.major_version, self.minor_version)


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
        return match[0].replace(match[1], new_version_number.as_string())

    modified = re.sub(pattern, get_replacement, source, flags=re.MULTILINE)

    return modified

    # """Apply all fixes in this module"""
    # modified = copy.deepcopy(source)
    #
    # for replacement in replacements:
    #     if replacement["new"] is None:
    #         print(
    #             "'%s' has been removed, please delete from your code"
    #             % replacement["old"]
    #         )
    #         continue
    #
    #     modified = fix_include_version_header(
    #         replacement["old"], replacement["headers"], modified
    #     )
    #     if replacement["macro"] and replacement["always_defined"]:
    #         modified = fix_always_defined_macros(
    #             replacement["old"], replacement["new"], modified
    #         )
    #     elif replacement["always_defined"]:
    #         modified = fix_ifdefs(replacement["old"], modified)
    #     modified = fix_replacement(replacement["old"], replacement["new"], modified)
    #
    # return modified


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
            TODO: Description here...
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

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    bump_version_numbers(new_version_number=VersionNumber(63, 15, 12))