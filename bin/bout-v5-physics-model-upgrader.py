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

PHYSICS_MODEL_RHS_SKELETON = "int {function}({arguments}){override};"

BOUTMAIN = "\n\nBOUTMAIN({})\n"

# Regular expression for a PhysicsModel
PHYSICS_MODEL_RE = re.compile(
    r"""class\s+(?P<name>[a-zA-Z0-9_]+)\s*:   # Class name
    \s*(?:public)?\s*PhysicsModel[\n\s]*{     # Inherits from PhysicsModel
    """,
    re.VERBOSE | re.MULTILINE,
)

FUNCTION_SIGNATURE_ARGUMENT_RE = r"""({arg_type}
    \s+                            # Require spaces only if the argument is named
    (?P<unused{arg_num}>UNUSED\()? # Possible UNUSED macro
    [a-zA-Z_0-9]*                  # Argument name
    (?(unused{arg_num})\))         # If UNUSED macro was present, we need an extra closing bracket
    )?
"""


def create_function_signature_re(function_name, argument_types):
    """Create a regular expression for a legacy physics model function"""

    if not isinstance(argument_types, list):
        argument_types = [argument_types]

    arguments = r",\s*".join(
        [
            FUNCTION_SIGNATURE_ARGUMENT_RE.format(arg_type=argument, arg_num=num)
            for num, argument in enumerate(argument_types)
        ]
    )

    return fr"int\s+{function_name}\s*\({arguments}\)"


LEGACY_MODEL_INCLUDE_RE = re.compile(
    r'^#\s*include.*(<|")boutmain.hxx(>|")', re.MULTILINE
)

BOUT_SOLVE_RE = re.compile(
    r"bout_solve\(([^,)]+,\s*[^,)]+(,\s*[^,)]+)?)\)", re.MULTILINE
)

RHS_RE = re.compile(r"solver\s*->\s*setRHS\(\s*([a-zA-Z0-9_]+)\s*\)")

PRECON_RE = re.compile(r"solver\s*->\s*setPrecon\(\s*([a-zA-Z0-9_]+)\s*\)")

JACOBIAN_RE = re.compile(r"solver\s*->\s*setJacobian\(\s*([a-zA-Z0-9_]+)\s*\)")

SPLIT_OPERATOR_RE = re.compile(
    r"solver\s*->\s*setSplitOperator\(\s*([a-zA-Z0-9_]+),\s*([a-zA-Z0-9_]+)\s*\)"
)


def has_split_operator(source):
    """Return the names of the split operator functions if set, otherwise False"""

    match = SPLIT_OPERATOR_RE.search(source)
    if not match:
        return False

    return match.group(1), match.group(2)


def is_legacy_model(source):
    """Return true if the source is a legacy physics model"""
    return LEGACY_MODEL_INCLUDE_RE.search(source) is not None


def find_last_include(source_lines):
    """Return the line number after the last #include (or 0 if no includes)"""
    for number, line in enumerate(reversed(source_lines)):
        if line.startswith("#include"):
            return len(source_lines) - number
    return 0


def fix_model_operator(
    source, model_name, operator_name, operator_type, new_name, override
):
    """Fix any definitions of the operator, and return the new declaration

    May modify source

    Parameters
    ----------
    source: str
        Source code to fix
    model_name: str
        Name of the PhysicsModel class to create or add methods to
    operator_name: str
        Name of the free function to fix
    operator_type: str, [str]
        Function argument types
    new_name: str
        Name of the PhysicsModel method
    override: bool
        Is `new_name` overriding a virtual?
    """

    # Make sure we have a list of types
    if not isinstance(operator_type, list):
        operator_type = [operator_type]

    # Get a regex for the function signature
    operator_re = re.compile(
        create_function_signature_re(operator_name, operator_type),
        re.VERBOSE | re.MULTILINE,
    )

    # Find any declarations of the free function
    matches = list(operator_re.finditer(source))
    if matches == []:
        warnings.warn(
            f"Could not find {operator_name}; is it defined in another file? If so, you will need to fix it manually"
        )
        return source, False

    # If we found more than one, remove the first one as it's probably
    # a declaration and not a definition
    if len(matches) > 1:
        source = re.sub(
            create_function_signature_re(operator_name, operator_type) + r"\s*;",
            "",
            source,
            flags=re.VERBOSE | re.MULTILINE,
        )

        # Get the names of the function arguments. Every other group
        # from the regex, as the other groups match the `UNUSED` macro
        arg_names = operator_re.search(source).groups()[::2]
    else:
        arg_names = matches[0].groups()[::2]

    # Fix definition and any existing declarations
    arguments = ", ".join(arg_names)

    # Modify the definition: it's out-of-line so we need the qualified name
    modified = operator_re.sub(fr"int {model_name}::{new_name}({arguments})", source)

    # Create the declaration
    return (
        modified,
        PHYSICS_MODEL_RHS_SKELETON.format(
            function=new_name,
            arguments=arguments,
            override=" override" if override else "",
        ),
    )


def fix_bout_constrain(source, error_on_warning):
    """Fix uses of bout_constrain. This is complicated because it might be
    in a conditional, and Solver::constraint returns void

    """

    if "bout_constrain" not in source:
        return source

    # The bout_constrain free function returns False if the Solver
    # doesn't have constraints, but the Solver::constraint method does
    # the checking itself, so we don't need to repeat it
    modified = re.sub(
        r"""if\s*\(\s*(?:!|not)\s*                   # in a conditional, checking for false
            bout_constrain\(([^;]+,[^;]+,[^;]+)\)        # actual function call
            \s*\)                                        # end of conditional
            (?P<brace>\s*{\s*)?                          # possible open brace
            (?:\s*\n)?                                   # possible newline
            \s*throw\s+BoutException\(.*\);(?:\s*\n)?    # throwing an exception
            (?(brace)\s*})?                              # consume matching closing brace
            """,
        r"solver->constraint(\1);\n",
        source,
        flags=re.VERBOSE | re.MULTILINE,
    )

    # The above might not fix everything, so best check if there are any uses left
    remaining_matches = list(re.finditer("bout_constrain", modified))
    if remaining_matches == []:
        # We fixed everything!
        return modified

    # Construct a useful error message
    source_lines = source.splitlines()
    lines_context = []
    for match in remaining_matches:
        bad_line = source[: match.end()].count("\n")
        line_range = range(max(0, bad_line - 1), min(len(source_lines), bad_line + 2))
        lines_context.append(
            "\n  ".join(["{}:{}".format(i, source_lines[i]) for i in line_range])
        )

    message = textwrap.dedent(
        """\
        Some uses of `bout_constrain` remain, but we could not automatically
        convert them to use `Solver::constraint`. Please fix them before
        continuing:
        """
    )
    message += "  " + "\n  ".join(lines_context)

    if error_on_warning:
        raise RuntimeError(message)
    print(message)
    return modified


def convert_old_solver_api(source, name):
    """Fix or remove old Solver API calls

    Parameters
    ----------
    source: str
        The source code to modify
    name: str
        The PhysicsModel class name
    """

    # Fixing `bout_solve` is a straight replacement, easy
    source = BOUT_SOLVE_RE.sub(r"solver->add(\1)", source)

    # Completely remove calls to Solver::setRHS
    source = RHS_RE.sub("", source)

    # List of old free functions that now need declarations inside the
    # class definition
    method_decls = []

    # Fix uses of solver->setPrecon
    # Get the name of any free functions passed as arguments to setPrecon
    precons = PRECON_RE.findall(source)
    for precon in precons:
        source, decl = fix_model_operator(
            source,
            name,
            precon,
            ["BoutReal", "BoutReal", "BoutReal"],
            precon,
            override=False,
        )
        if decl:
            method_decls.append(decl)
    # Almost a straight replacement, but it's now a member-function pointer
    source = PRECON_RE.sub(fr"setPrecon(&{name}::\1)", source)

    # Fix uses of solver->setJacobian, basically the same as for setPrecon
    jacobians = JACOBIAN_RE.findall(source)
    for jacobian in jacobians:
        source, decl = fix_model_operator(
            source, name, jacobian, "BoutReal", jacobian, override=False
        )
        if decl:
            method_decls.append(decl)
    source = JACOBIAN_RE.sub(fr"setJacobian(&{name}::\1)", source)

    # If we didn't find any free functions that need to be made into
    # methods, we're done
    if not method_decls:
        return source

    # We need to find the class defintion
    class_def = PHYSICS_MODEL_RE.search(source)
    if class_def is None:
        warnings.warn(
            f"Could not find the '{name}' class to add"
            "preconditioner and/or Jacobian declarations; is it defined"
            "in another file? If so, you will need to fix it manually"
        )
        return source, False

    # The easiest place to stick the method declaration is on the line
    # immediately following the open brace of the class def, and the
    # easiest way to insert it is to split the source into a list,
    # insert in the list, then join the list back into a string.
    # The regex from above finds the offset in the source which we
    # need to turn into a line number
    first_line_of_class = source[: class_def.end() + 1].count("\n")
    methods = "\n  ".join(method_decls)
    source_lines = source.splitlines()
    source_lines.insert(first_line_of_class, f"  {methods}")

    return "\n".join(source_lines)


def convert_legacy_model(source, name, error_on_warning):
    """Convert a legacy physics model to a PhysicsModel"""

    if not is_legacy_model(source):
        return source

    source = fix_bout_constrain(source, error_on_warning)

    # Replace legacy header
    source = LEGACY_MODEL_INCLUDE_RE.sub(r"#include \1bout/physicsmodel.hxx\2", source)

    method_decls = []

    source, decl = fix_model_operator(
        source, name, "physics_init", "bool", "init", override=True
    )
    if decl:
        method_decls.append(decl)

    split_operators = has_split_operator(source)
    if split_operators:
        source = SPLIT_OPERATOR_RE.sub(r"setSplitOperator(true)", source)

        convective, diffusive = split_operators
        # Fix the free functions
        source, decl = fix_model_operator(
            source, name, convective, "BoutReal", "convective", override=True
        )
        if decl:
            method_decls.append(decl)
        source, decl = fix_model_operator(
            source, name, diffusive, "BoutReal", "diffusive", override=True
        )
        if decl:
            method_decls.append(decl)
    else:
        # Fix the rhs free function
        source, decl = fix_model_operator(
            source, name, "physics_run", "BoutReal", "rhs", override=True
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
            Upgrade legacy physics models to use the PhysicsModel class

            This will do the bare minimum required to compile, and
            won't make global objects (like Field3Ds) members of the
            new class, or free functions (other than
            `physics_init`/`physics_run`, preconditioners, and
            Jacobians) methods of the new class. Comments may also be
            left behind.

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
        "--patch-only", "-p", action="store_true", help="Print the patches and exit"
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

        original = copy.deepcopy(contents)

        match = PHYSICS_MODEL_RE.search(original)
        if match is not None:
            new_name = match.group("name")
        else:
            new_name = args.name or pathlib.Path(filename).stem.capitalize().replace(
                " ", "_"
            )

        try:
            if re.match(r"^[0-9]+.*", new_name) and not args.force:
                raise ValueError(
                    f"Invalid name: '{new_name}'. Use --name to specify a valid C++ identifier"
                )

            modified = convert_legacy_model(
                original, new_name, not (args.force or args.patch_only)
            )

            modified = convert_old_solver_api(modified, new_name)
        except (RuntimeError, ValueError) as e:
            error_message = textwrap.indent(f"{e}", " ")
            print(
                f"There was a problem applying automatic fixes to {filename}:\n\n{error_message}"
            )
            continue

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

        make_change = args.force or yes_or_no("Make changes to {}?".format(filename))

        if make_change:
            with open(filename, "w") as f:
                f.write(modified)
