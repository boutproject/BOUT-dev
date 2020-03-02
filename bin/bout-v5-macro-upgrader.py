#!/usr/bin/env python3

import argparse
import copy
import difflib
import re
import textwrap

# List of macros, their replacements and what header to find them
# in. Each element should be a dict with "old", "new" and "headers"
# keys, with "old" and "new" values being strings, and "headers" being a
# list of strings
MACRO_REPLACEMENTS = [
    {
        "old": "REVISION",
        "new": "bout::version::revision",
        "headers": ["bout/revision.hxx"],
    },
    {
        "old": "BOUT_VERSION_DOUBLE",
        "new": "bout::version::as_double",
        "headers": ["bout/version.hxx", "bout.hxx"],
    },
    {
        "old": "BOUT_VERSION_STRING",
        "new": "bout::version::full",
        "headers": ["bout/version.hxx", "bout.hxx"],
    },
    # Next one is not technically a macro, but near enough
    {
        "old": "BOUT_VERSION",
        "new": "bout::version::full",
        "headers": ["bout/version.hxx", "bout.hxx"],
    },
]


def fix_include_version_header(old, headers, source):
    """Make sure version.hxx header is included
    """

    # If header is already included, we can skip this fix
    for header in headers:
        if (
            re.search(r'^#\s*include.*(<|"){}(>|")'.format(header), source, flags=re.M)
            is not None
        ):
            return source

    # If the old macro isn't in the file, we can skip this fix
    if re.search(r"\b{}\b".format(old), source) is None:
        return source

    # Now we want to find a suitable place to stick the new include
    # Good candidates are includes of BOUT++ headers
    includes = []
    source_lines = source.splitlines()
    for linenumber, line in enumerate(source_lines):
        if re.match(r"^#\s*include.*bout/", line):
            includes.append(linenumber)
        if re.match(r"^#\s*include.*physicsmodel", line):
            includes.append(linenumber)

    if includes:
        last_include = includes[-1] + 1
    else:
        # No suitable includes, so just stick at the top of the file
        last_include = 0
    source_lines.insert(last_include, '#include "{}"'.format(headers[0]))

    return "\n".join(source_lines)


def fix_ifdefs(old, source):
    """Remove any code inside #ifdef/#ifndef blocks that would now not be compiled

    """
    source_lines = source.splitlines()

    # Something to keep track of nested sections
    in_ifdef = None
    # List of (#ifdef or #ifndef, dict of start/else/end lines)
    macro_blocks = []
    for linenumber, line in enumerate(source_lines):
        if_def = re.match(r"#\s*(ifn?def)\s*(.*)", line)
        else_block = re.match(r"#\s*else", line)
        endif = re.match(r"#\s*endif", line)
        if not (if_def or else_block or endif):
            continue
        # Now we need to keep track of whether we're inside an
        # interesting #ifdef/ifndef, as they might be nested, and we
        # want to find the matching #endif and #else
        if endif:
            if in_ifdef is not None:
                in_ifdef -= 1
            if in_ifdef == 0:
                in_ifdef = None
                macro_blocks[-1]["end"] = linenumber
            continue
        if else_block:
            if in_ifdef == 1:
                macro_blocks[-1]["else"] = linenumber
            continue
        if if_def.group(2) == old:
            in_ifdef = 1
            macro_blocks.append({"start": linenumber, "if_def_type": if_def.group(1)})
        elif in_ifdef is not None:
            in_ifdef += 1

    if macro_blocks == []:
        return source

    # Get all of the lines to be removed
    lines_to_remove = set()
    for block in macro_blocks:
        lines_to_remove |= set(block.values())
        if block["if_def_type"] == "ifdef":
            if "else" in block:
                # Delete the #else block for #ifdef
                lines_to_remove |= set(range(block["else"], block["end"]))
        else:
            # Keep the #else block for #ifndef if there is one, otherwise remove the
            # whole block
            lines_to_remove |= set(
                range(block["start"], block.get("else", block["end"]))
            )

    # Apparently this is actually the best way of removing a bunch of (possibly)
    # non-contiguous indices
    modified_lines = [
        line for num, line in enumerate(source_lines) if num not in lines_to_remove
    ]

    return "\n".join(modified_lines)


def fix_replacement(old, new, source):
    """Straight replacements
    """
    return re.sub(r'([^"])\b{}\b([^"])'.format(old), r"\1{}\2".format(new), source)


def apply_fixes(replacements, source):
    """Apply all fixes in this module
    """
    modified = copy.deepcopy(source)

    for replacement in replacements:
        modified = fix_include_version_header(
            replacement["old"], replacement["headers"], modified
        )
        modified = fix_ifdefs(replacement["old"], modified)
        modified = fix_replacement(replacement["old"], replacement["new"], modified)

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
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Fix macro defines for BOUT++ v4 -> v5

            Please note that this is only slightly better than dumb text replacement. It
            will fix the following:

            * replacement of macros with variables
            * inclusion of correct headers for new variables
            * removal of #if(n)def/#endif blocks that do simple checks for the old
              macro, keeping the appriopriate part

            It will try not to replace quoted macro names.

            Please check the diff output carefully!
            """
        ),
    )

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

        modified = apply_fixes(MACRO_REPLACEMENTS, contents)
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
