#!/usr/bin/env python3


import argparse
import copy
import difflib
import re

try:
    import clang.cindex

    has_clang = True
except ImportError:
    has_clang = False


headers = {"interpolation": {"old": "interpolation.hxx", "new": "interpolation_xz.hxx"}}

interpolations = {
    "Hermite": {"old": "HermiteSpline", "new": "XZHermiteSpline"},
    "Interpolation": {"old": "Interpolation", "new": "XZInterpolation"},
    "MonotonicHermite": {
        "old": "MonotonicHermiteSpline",
        "new": "XZMonotonicHermiteSpline",
    },
    "Bilinear": {"old": "Bilinear", "new": "XZBilinear"},
    "Lagrange4pt": {"old": "Lagrange4pt", "new": "XZLagrange4pt"},
}

factories = {
    "InterpolationFactory": {
        "old": "InterpolationFactory",
        "new": "XZInterpolationFactory",
    }
}


def fix_header_includes(old_header, new_header, source):
    """Replace old_header with new_header in source

    Parameters
    ----------
    old_header: str
        Name of header to be replaced
    new_header: str
        Name of replacement header
    source: str
        Text to search

    """
    return re.sub(
        r"""
        (\s*\#\s*include\s*)     # Preprocessor include
        (<|")
        ({header})              # Header name
        (>|")
        """.format(
            header=old_header
        ),
        r"\1\2{header}\4".format(header=new_header),
        source,
        flags=re.VERBOSE,
    )


def fix_interpolations(old_interpolation, new_interpolation, source):

    return re.sub(
        r"""
        \b{}\b
        """.format(
            old_interpolation
        ),
        r"{}".format(new_interpolation),
        source,
        flags=re.VERBOSE,
    )


def clang_parse(filename, source):
    index = clang.cindex.Index.create()
    return index.parse(filename, unsaved_files=[(filename, source)])


def clang_find_interpolations(node, typename, nodes=None):
    if nodes is None:
        nodes = []
    if node.kind == clang.cindex.CursorKind.TYPE_REF:
        if node.type.spelling == typename:
            nodes.append(node)
    for child in node.get_children():
        clang_find_interpolations(child, typename, nodes)
    return nodes


def clang_fix_single_interpolation(
    old_interpolation, new_interpolation, source, location
):
    modified = source
    line = modified[location.line - 1]
    new_line = (
        line[: location.column - 1]
        + new_interpolation
        + line[location.column + len(old_interpolation) - 1 :]
    )
    modified[location.line - 1] = new_line
    return modified


def clang_fix_interpolation(old_interpolation, new_interpolation, node, source):
    nodes = clang_find_interpolations(node, old_interpolation)
    modified = source
    for node in nodes:
        modified = clang_fix_single_interpolation(
            old_interpolation, new_interpolation, modified, node.location
        )
    return modified


def fix_factories(old_factory, new_factory, source):

    return re.sub(
        r"""
        \b{}\b
        """.format(
            old_factory
        ),
        r"{}".format(new_factory),
        source,
        flags=re.VERBOSE,
    )


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


def yes_or_no(question):
    """Convert user input from yes/no variations to True/False

    """
    while True:
        reply = input(question + " [y/N] ").lower().strip()
        if not reply or reply[0] == "n":
            return False
        if reply[0] == "y":
            return True


def apply_fixes(headers, interpolations, factories, source):
    """Apply all Interpolation fixes to source

    Parameters
    ----------
    headers
        Dictionary of old/new headers
    interpolations
        Dictionary of old/new Interpolation types
    source
        Text to update

    """

    modified = copy.deepcopy(source)

    for header in headers.values():
        modified = fix_header_includes(header["old"], header["new"], modified)
    for interpolation in interpolations.values():
        modified = fix_interpolations(
            interpolation["old"], interpolation["new"], modified
        )
    for factory in factories.values():
        modified = fix_factories(factory["old"], factory["new"], modified)

    return modified


def clang_apply_fixes(headers, interpolations, factories, filename, source):

    # translation unit
    tu = clang_parse(filename, source)

    modified = source

    for header in headers.values():
        modified = fix_header_includes(header["old"], header["new"], modified)

    modified = modified.split("\n")
    for interpolation in interpolations.values():
        modified = clang_fix_interpolation(
            interpolation["old"], interpolation["new"], tu.cursor, modified
        )
    modified = "\n".join(modified)
    for factory in factories.values():
        modified = fix_factories(factory["old"], factory["new"], modified)

    return modified


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix types of Interpolation objects")

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
    parser.add_argument(
        "--clang", action="store_true", help="Use libclang if available"
    )

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    if args.clang and not has_clang:
        raise RuntimeError(
            "libclang is not available. Please install libclang Python bindings"
        )

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()
        original = copy.deepcopy(contents)

        if args.clang and has_clang:
            modified = clang_apply_fixes(
                headers, interpolations, factories, filename, contents
            )
        else:
            modified = apply_fixes(headers, interpolations, factories, contents)
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
