#!/bin/bash

PROJECT_ROOT="@PROJECT_BINARY_DIR@"

export PYTHONPATH="@PROJECT_BINARY_DIR@/tools/pylib:@PROJECT_SOURCE_DIR@/tools/pylib:$PYTHONPATH"

# Pre-build the project to prevent concurrent CMake race conditions
cmake --build "$PROJECT_ROOT"

python3 -m pytest -m "not serial" --cache-clear -n auto --dist=loadgroup -q "$PROJECT_ROOT/tests/integrated" "$@"
python3 -m pytest -m serial "$PROJECT_ROOT/tests/integrated" "$@"
