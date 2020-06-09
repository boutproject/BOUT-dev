#!/usr/bin/env python3

import argparse
import copy
import difflib
import pathlib
import re
import textwrap


PHYSICS_MODEL_INCLUDE = '#include "bout/physicsmodel.hxx"'

PHYSICS_MODEL_SKELETON = """
class {name} : public PhysicsModel {{
protected:
  int init(bool{restarting}) override;
  int rhs(BoutReal{time}) override;
}};
"""

BOUTMAIN = "\n\nBOUTMAIN({})\n"

# Regular expression for a PhysicsModel
PHYSICS_MODEL_RE = re.compile(r":\s*public\s*PhysicsModel")

# Regular expressions for a legacy physics model
LEGACY_MODEL_RE_TEMPLATE = r"""int\s+physics_{}\s*\({}
    (\s+                           # Require spaces only if the argument is named
    (?P<unused>UNUSED\()?          # Possible UNUSED macro
    [a-zA-Z_0-9]*                  # Argument name
    (?(unused)\))                  # If UNUSED macro was present, we need an extra closing bracket
    )?
    \)"""

LEGACY_MODEL_INIT_RE = re.compile(
    LEGACY_MODEL_RE_TEMPLATE.format("init", "bool"), re.VERBOSE | re.MULTILINE
)
LEGACY_MODEL_RUN_RE = re.compile(
    LEGACY_MODEL_RE_TEMPLATE.format("run", "BoutReal"), re.VERBOSE | re.MULTILINE
)

LEGACY_MODEL_INCLUDE_RE = re.compile(
    r'^#\s*include.*(<|")boutmain.hxx(>|")', re.MULTILINE
)

BOUT_CONSTRAIN_RE = re.compile(r"bout_constrain\(([^,)]+,\s*[^,)]+,\s*[^,)]+)\)")
BOUT_SOLVE_RE = re.compile(r"bout_solve\(([^,)]+,\s*[^,)]+)\)")


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


def convert_legacy_model(source, name):
    """Convert a legacy physics model to a PhysicsModel
    """

    if not is_legacy_model(source):
        return source

    # Replace legacy header
    replaced_header = LEGACY_MODEL_INCLUDE_RE.sub(
        r"#include \1bout/physicsmodel.hxx\2", source
    )

    source_lines = replaced_header.splitlines()
    last_include = find_last_include(source_lines)

    init_function = LEGACY_MODEL_INIT_RE.search(source)
    if init_function is not None:
        restarting = init_function.group(1)
    else:
        restarting = ""

    run_function = LEGACY_MODEL_RUN_RE.search(source)
    if run_function is not None:
        time = run_function.group(1)
    else:
        time = ""

    source_lines.insert(
        last_include,
        PHYSICS_MODEL_SKELETON.format(name=name, restarting=restarting, time=time),
    )

    added_class = "\n".join(source_lines)

    fixed_init = LEGACY_MODEL_INIT_RE.sub(
        r"int {}::init(bool\1)".format(name), added_class
    )
    fixed_run = LEGACY_MODEL_RUN_RE.sub(
        r"int {}::rhs(BoutReal\1)".format(name), fixed_init
    )

    fixed_constraint = BOUT_CONSTRAIN_RE.sub(r"solver->constraint(\1)", fixed_run)
    fixed_solve = BOUT_CONSTRAIN_RE.sub(r"solver->add(\1)", fixed_constraint)

    added_main = fixed_solve + BOUTMAIN.format(name)
    return added_main


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
