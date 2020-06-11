#!/usr/bin/env python3

import argparse
import copy
import difflib
import pathlib
import re
import textwrap
import warnings


PHYSICS_MODEL_INCLUDE = '#include "bout/physicsmodel.hxx"'

PHYSICS_MODEL_SKELETON = """
class {name} : public PhysicsModel {{
protected:
  {methods}
}};
"""

PHYSICS_MODEL_RHS_SKELETON = "int {function}({arg_type}{time}) override;"

BOUTMAIN = "\n\nBOUTMAIN({})\n"

# Regular expression for a PhysicsModel
PHYSICS_MODEL_RE = re.compile(r":\s*public\s*PhysicsModel")

# Regular expressions for a legacy physics model
LEGACY_MODEL_RE_TEMPLATE = r"""int\s+{}\s*\({}
    (\s+                           # Require spaces only if the argument is named
    (?P<unused>UNUSED\()?          # Possible UNUSED macro
    [a-zA-Z_0-9]*                  # Argument name
    (?(unused)\))                  # If UNUSED macro was present, we need an extra closing bracket
    )?
    \)"""

LEGACY_MODEL_INCLUDE_RE = re.compile(
    r'^#\s*include.*(<|")boutmain.hxx(>|")', re.MULTILINE
)

SPLIT_OPERATOR_RE = re.compile(
    r"solver\s*->\s*setSplitOperator\(([a-zA-Z0-9_]+),\s*([a-zA-Z0-9_]+)\s*\)"
)


def has_split_operator(source):
    """Return the names of the split operator functions if set, otherwise False

    """

    match = SPLIT_OPERATOR_RE.search(source)
    if not match:
        return False

    return match.group(1), match.group(2)


def is_legacy_model(source):
    """Return true if the source is a legacy physics model

    """
    return LEGACY_MODEL_INCLUDE_RE.search(source) is not None


def find_last_include(source_lines):
    """Return the line number after the last #include (or 0 if no includes)
    """
    for number, line in enumerate(reversed(source_lines)):
        if line.startswith("#include"):
            return len(source_lines) - number
    return 0


def fix_model_operator(source, model_name, operator_name, operator_type, new_name):
    """Fix any definitions of the operator, and return the new declaration

    May modify source
    """

    operator_re = re.compile(
        LEGACY_MODEL_RE_TEMPLATE.format(operator_name, operator_type),
        re.VERBOSE | re.MULTILINE,
    )

    matches = list(operator_re.finditer(source))
    if matches == []:
        warnings.warn(
            f"Could not find {operator_name}; is it defined in another file? If so, you will need to fix it manually"
        )
        return source, False

    if len(matches) > 1:
        source = re.sub(
            LEGACY_MODEL_RE_TEMPLATE.format(operator_name, operator_type) + r"\s*;",
            "",
            source,
            flags=re.VERBOSE | re.MULTILINE,
        )

        arg_name = operator_re.search(source).group(1)
    else:
        arg_name = matches[0].group(1)

    # Fix definition and any existing declarations
    modified = operator_re.sub(
        fr"int {model_name}::{new_name}({operator_type}\1)", source
    )

    # Create the declaration
    return (
        modified,
        PHYSICS_MODEL_RHS_SKELETON.format(
            function=new_name, arg_type=operator_type, time=arg_name
        ),
    )


def convert_legacy_model(source, name):
    """Convert a legacy physics model to a PhysicsModel
    """

    if not is_legacy_model(source):
        return source

    # Replace legacy header
    source = LEGACY_MODEL_INCLUDE_RE.sub(r"#include \1bout/physicsmodel.hxx\2", source)

    method_decls = []

    source, decl = fix_model_operator(source, name, "physics_init", "bool", "init")
    if decl:
        method_decls.append(decl)

    split_operators = has_split_operator(source)
    if split_operators:
        source = SPLIT_OPERATOR_RE.sub(r"setSplitOperator(true)", source)

        convective, diffusive = split_operators
        # Fix the free functions
        source, decl = fix_model_operator(
            source, name, convective, "BoutReal", "convective"
        )
        if decl:
            method_decls.append(decl)
        source, decl = fix_model_operator(
            source, name, diffusive, "BoutReal", "diffusive"
        )
        if decl:
            method_decls.append(decl)
    else:
        # Fix the rhs free function
        source, decl = fix_model_operator(
            source, name, "physics_run", "BoutReal", "rhs"
        )
        if decl:
            method_decls.append(decl)

    source_lines = source.splitlines()
    last_include = find_last_include(source_lines)

    methods = "\n  ".join(method_decls)
    physics_model = PHYSICS_MODEL_SKELETON.format(methods=methods, name=name)

    source_lines.insert(last_include, physics_model)
    source_lines.append(BOUTMAIN.format(name))

    return "\n".join(source_lines)


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
            Upgrade legacy physics models to use the PhysicsModel class

            This will do the bare minimum required to compile, and
            won't make global objects (like Field3Ds) members of the
            new class, or free functions other than
            `physics_init`/`physics_run` methods of the new class.

            By default, this will use the file name stripped of file
            extensions as the name of the new class. Use '--name=<new
            name>' to give a different name.
            """
        ),
    )

    parser.add_argument("files", action="store", nargs="+", help="Files to fix")
    parser.add_argument(
        "--force", "-f", action="store_true", help="Make changes without asking"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Don't print patches"
    )
    parser.add_argument(
        "--patch-only", "-p", action="store_true", help="Print the patches and exist"
    )
    parser.add_argument(
        "--name",
        "-n",
        action="store",
        nargs="?",
        type=str,
        help="Name for new PhysicsModel class, default is from filename",
    )

    args = parser.parse_args()

    if args.force and args.patch_only:
        raise ValueError("Incompatible options: --force and --patch")

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()

        if not is_legacy_model(contents):
            if not args.quiet:
                print("No changes to make to {}".format(filename))
            continue

        original = copy.deepcopy(contents)

        new_name = args.name or pathlib.Path(filename).stem.capitalize().replace(
            " ", "_"
        )

        if re.match(r"^[0-9]+.*", new_name):
            raise ValueError(
                f"Invalid name: '{new_name}'. Use --name to specify a valid C++ identifier"
            )

        modified = convert_legacy_model(original, new_name)
        patch = create_patch(filename, original, modified)

        if args.patch_only:
            print(patch)
            continue

        if not args.quiet:
            print("\n******************************************")
            print("Changes to {}\n".format(filename))
            print(patch)
            print("\n******************************************")

        make_change = args.force or yes_or_no("Make changes to {}?".format(filename))

        if make_change:
            with open(filename, "w") as f:
                f.write(modified)
