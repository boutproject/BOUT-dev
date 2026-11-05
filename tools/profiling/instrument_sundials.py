#!/usr/bin/env python3
"""Add balanced NVTX scopes to a SUNDIALS 7.9.0 source tree before rebuilding.

Scopes use the NVTX bridge exported by BOUT++ built with USE_NVTX=ON.
No NVTX dependency or solver algorithm changes are introduced in SUNDIALS.
Requires GCC/Clang on ELF platforms. Run again safely to check an existing patch.
"""

import argparse
import re
import shutil
from pathlib import Path

HEADER = "sundials_bout_nvtx.h"
MARKER = "SUNDIALS_BOUT_NVTX_SCOPE"
# Ignore braces inside comments and literals while locating function bodies.
C_TOKEN = re.compile(
    r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}]', re.DOTALL
)


def function_body(source, name):
    match = re.search(
        r"\b(?:static\s+)?(?:int|void|sunrealtype|sunbooleantype|SUNErrCode)\s+"
        + re.escape(name)
        + r"\s*\([^;{}]*\)\s*\{",
        source,
    )
    if match is None:
        raise ValueError(f"function definition not found: {name}")
    start = match.end()
    depth = 1
    for token in C_TOKEN.finditer(source, start):
        if token.group() == "{":
            depth += 1
        elif token.group() == "}":
            depth -= 1
            if not depth:
                return start, token.start()
    raise ValueError(f"unterminated function: {name}")


def instrument(source, function, label, loop=None):
    start, end = function_body(source, function)
    body = source[start:end]
    if loop is None:
        body = f'\n  {MARKER}(bout_nvtx_function, "SUNDIALS::{label}");' + body
    else:
        pattern, occurrence = loop
        matches = list(re.finditer(pattern, body))
        if len(matches) <= occurrence:
            raise ValueError(f"expected loop not found in {function}: {pattern}")
        pos = matches[occurrence].end()
        body = (
            body[:pos] + f"\n    {MARKER}(bout_nvtx_iteration, "
            f'"SUNDIALS::{label}");' + body[pos:]
        )
    return source[:start] + body + source[end:]


def patch(root):
    version = (root / "CMakeLists.txt").read_text()
    for part, value in (("MAJOR", "7"), ("MINOR", "9"), ("PATCH", "0")):
        if not re.search(
            r"set\(PACKAGE_VERSION_" + part + r'\s+"?' + value + r'"?\s*\)', version
        ):
            raise ValueError("this patch supports SUNDIALS 7.9.0 only")
    specs = {
        "src/cvode/cvode.c": [
            ("CVode", "CVode", None),
            ("cvStep", "CVStep", None),
            ("cvStep", "StepAttempt", (r"for\s*\(;;\)\s*\{", 0)),
            ("cvNls", "CVNls", None),
        ],
        "src/cvode/cvode_ls.c": [
            ("cvLsSolve", "CVLsSolve", None),
            ("cvLsATimes", "CVLsATimes", None),
            ("cvLsDQJtimes", "CVLsDQJtimes", None),
            ("cvLsSetup", "CVLsSetup", None),
            ("cvLsPSolve", "CVLsPSolve", None),
        ],
        "src/sunnonlinsol/newton/sunnonlinsol_newton.c": [
            ("SUNNonlinSolSolve_Newton", "NewtonSolve", None),
            ("SUNNonlinSolSolve_Newton", "NewtonIteration", (r"for\s*\(;;\)\s*\{", 1)),
        ],
        "src/sunlinsol/spgmr/sunlinsol_spgmr.c": [
            ("SUNLinSolSolve_SPGMR", "SPGMRSolve", None),
            (
                "SUNLinSolSolve_SPGMR",
                "GMRESIteration",
                (r"for\s*\(l = 0; l < l_max; l\+\+\)\s*\{", 0),
            ),
        ],
    }
    # Instrument the public dispatch layer, including fused fallbacks and local
    # reductions. This covers CPU and CUDA/MPIPlusX without duplicate backend
    # scopes. Keep lifecycle/accessor and buffer-transfer functions separate.
    vector_file = "src/sundials/sundials_nvector.c"
    vector_source = (root / vector_file).read_text()
    arithmetic = vector_source[
        vector_source.index("void N_VLinearSum(") : vector_source.index(
            "SUNErrCode N_VBufSize("
        )
    ]
    vector_functions = re.findall(
        r"\b(?:void|sunrealtype|sunbooleantype|SUNErrCode)\s+"
        r"(N_V\w+)\s*\([^;{}]*\)\s*\{",
        arithmetic,
    )
    if not vector_functions or len(vector_functions) != len(set(vector_functions)):
        raise ValueError("unexpected N_Vector arithmetic definitions")
    specs[vector_file] = [(name, name, None) for name in vector_functions]
    # Validate every file before modifying the tree.
    updates = {}
    for filename, scopes in specs.items():
        path = root / filename
        source = path.read_text()
        if MARKER in source:
            if (
                filename == "src/cvode/cvode.c"
                and "SUNBoutNVTXInstrumented" not in source
            ):
                raise ValueError(f"missing runtime detection export in {path}")
            for _, label, _ in scopes:
                if f'"SUNDIALS::{label}"' not in source:
                    raise ValueError(f"partial instrumentation in {path}")
            continue
        for function, label, loop in scopes:
            source = instrument(source, function, label, loop)
        relative_header = "../" * (len(Path(filename).parts) - 2) + "sundials/" + HEADER
        first_include = source.index("#include")
        source = (
            source[:first_include]
            + f'#include "{relative_header}"\n'
            + source[first_include:]
        )
        if filename == "src/cvode/cvode.c":
            source += (
                "\n/* Runtime detection for BOUT++ profiling diagnostics. */\n"
                "SUNDIALS_EXPORT int SUNBoutNVTXInstrumented(void)\n{\n"
                "  return bout_nvtx_range_push && bout_nvtx_range_pop;\n}\n"
            )
        updates[path] = source
    shutil.copyfile(Path(__file__).with_name(HEADER), root / "src/sundials" / HEADER)
    for path, source in updates.items():
        path.write_text(source)
    print(f"Instrumented {len(updates)} source files in {root}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="SUNDIALS 7.9.0 source directory")
    args = parser.parse_args()
    try:
        patch(args.source.resolve())
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
