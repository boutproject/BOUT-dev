#!/bin/bash

PROJECT_ROOT=${PROJECT_BINARY_DIR}

export PYTHONPATH="$PROJECT_ROOT/tools/pylib"

# Pre-build the project to prevent concurrent CMake race conditions
cmake --build "$PROJECT_ROOT"

pytest -m "not serial" --cache-clear -n auto --dist=loadgroup -q "$PROJECT_ROOT/tests/integrated" "$@"
pytest -m serial "$PROJECT_ROOT/tests/integrated" "$@"
