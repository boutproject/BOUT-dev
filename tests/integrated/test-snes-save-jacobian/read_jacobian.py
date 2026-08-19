#!/usr/bin/env python3

"""Read Jacobian diagnostics written by the SNES save-jacobian test.

This helper is intentionally lightweight:

- PETSc binary matrices are read using the bundled ``PetscBinaryIO.py``
- ``jacobian_index_base`` is read with ``boutdata.collect`` so MPI-written dump
  files can be reconstructed in serial
- NumPy is the only required array dependency
- Pandas and SciPy are optional convenience layers
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent


def _discover_repo_root() -> Path | None:
    """Find the source-tree root from either the source or copied build script."""

    for candidate in (SCRIPT_DIR, *SCRIPT_DIR.parents):
        if (candidate / "tools" / "pylib").exists():
            return candidate
        if (
            candidate.name == "build"
            and (candidate.parent / "tools" / "pylib").exists()
        ):
            return candidate.parent
    return None


REPO_ROOT = _discover_repo_root()

# Allow this script to run directly from the test directory or build directory.
for extra_path in (
    SCRIPT_DIR,
    None if REPO_ROOT is None else REPO_ROOT / "tools" / "pylib",
    None if REPO_ROOT is None else REPO_ROOT / "build" / "tools" / "pylib",
):
    if extra_path is not None and extra_path.exists():
        sys.path.insert(0, str(extra_path))

import PetscBinaryIO  # noqa: E402


RawPetscSparse = tuple[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]


def _to_native_endian(array: np.ndarray) -> np.ndarray:
    """Return an array with native byte order for downstream Python libraries."""

    result = np.asarray(array)
    if result.dtype.byteorder in {"=", "|"}:
        return result
    return result.byteswap().view(result.dtype.newbyteorder("="))


def load_metadata(path: str | Path) -> dict[str, Any]:
    """Load Jacobian variable metadata from JSON."""

    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_index_base(data_dir: str | Path) -> np.ndarray:
    """Load ``jacobian_index_base`` from BOUT++ dump files using ``boutdata``."""

    try:
        from boutdata.collect import collect
    except ImportError as exc:
        raise ImportError(
            "Could not import boutdata.collect. Make sure tools/pylib and "
            "build/tools/pylib are on PYTHONPATH, and that the Python "
            "NetCDF dependencies required by boutdata (for example "
            "'netCDF4') are installed."
        ) from exc

    index_base = np.asarray(
        collect("jacobian_index_base", path=str(data_dir), info=False, xguards=True)
    )

    if index_base.ndim == 2:
        return index_base[:, :, np.newaxis]
    if index_base.ndim != 3:
        raise ValueError(
            f"Expected jacobian_index_base to have 2 or 3 dimensions, got {index_base.ndim}"
        )
    return index_base


def _find_matrix_file(data_dir: Path) -> Path:
    matrix_files = sorted(data_dir.glob("jacobian_*.dat"))
    if not matrix_files:
        raise FileNotFoundError(f"No PETSc binary matrix files found in {data_dir}")
    if len(matrix_files) > 1:
        raise ValueError(
            f"Found multiple Jacobian matrix files in {data_dir}; specify one explicitly"
        )
    return matrix_files[0]


def load_petsc_matrix(
    path: str | Path, matrix_format: str = "dense"
) -> np.ndarray | RawPetscSparse:
    """Load a PETSc binary matrix using the bundled ``PetscBinaryIO`` helper.

    Parameters
    ----------
    path:
        Path to a PETSc binary matrix file.
    matrix_format:
        Either ``'dense'`` or ``'sparse'``.
    """

    if matrix_format not in {"dense", "sparse"}:
        raise ValueError(
            f"matrix_format must be 'dense' or 'sparse', got {matrix_format!r}"
        )

    io = PetscBinaryIO.PetscBinaryIO(
        precision="double", indices="32bit", complexscalars=False
    )
    objects = io.readBinaryFile(str(path), mattype=matrix_format)
    if not objects:
        raise ValueError(f"No PETSc objects found in {path}")

    matrix = objects[0]
    if matrix_format == "dense":
        return _to_native_endian(np.asarray(matrix))
    return matrix


def sparse_to_dense(matrix: RawPetscSparse) -> np.ndarray:
    """Convert a raw PETSc CSR tuple into a dense NumPy array."""

    (nrows, ncols), (indptr, indices, values) = matrix
    dense = np.zeros((nrows, ncols), dtype=np.asarray(values).dtype)
    for row in range(nrows):
        start = int(indptr[row])
        end = int(indptr[row + 1])
        dense[row, np.asarray(indices[start:end], dtype=int)] = values[start:end]
    return _to_native_endian(dense)


def to_numpy_dense(matrix: np.ndarray | RawPetscSparse) -> np.ndarray:
    """Return a dense NumPy matrix regardless of the original representation."""

    if isinstance(matrix, np.ndarray):
        return np.asarray(matrix)
    return sparse_to_dense(matrix)


def to_scipy_csr(matrix: np.ndarray | RawPetscSparse):
    """Convert a PETSc sparse matrix representation to ``scipy.sparse.csr_matrix``."""

    try:
        from scipy.sparse import csr_matrix
    except ImportError as exc:
        raise ImportError("SciPy is required for CSR conversion") from exc

    if isinstance(matrix, np.ndarray):
        return csr_matrix(matrix)

    (nrows, ncols), (indptr, indices, values) = matrix
    return csr_matrix((values, indices, indptr), shape=(nrows, ncols))


def _normalise_dofs(dofs: Any) -> list[dict[str, Any]]:
    """Return DOF records as a plain list of dicts.

    This accepts either the native ``list[dict]`` returned by ``load_jacobian()``
    or a Pandas DataFrame created by ``to_pandas()``.
    """

    if hasattr(dofs, "to_dict"):
        return list(dofs.to_dict(orient="records"))
    return list(dofs)


def variable_indices(dofs: Any, variable_name: str) -> list[int]:
    """Return matrix row/column positions corresponding to one variable."""

    records = _normalise_dofs(dofs)
    positions = [
        position
        for position, record in enumerate(records)
        if record["name"] == variable_name
    ]
    if not positions:
        raise KeyError(f"Variable {variable_name!r} was not found in the DOF metadata")
    return positions


def _slice_raw_petsc_sparse(
    matrix: RawPetscSparse, row_indices: list[int], col_indices: list[int]
) -> RawPetscSparse:
    """Slice a raw PETSc CSR tuple and return the same sparse representation."""

    (_, _), (indptr, indices, values) = matrix
    col_lookup = {old_col: new_col for new_col, old_col in enumerate(col_indices)}
    sliced_indptr = np.zeros(len(row_indices) + 1, dtype=np.asarray(indptr).dtype)
    sliced_indices: list[int] = []
    sliced_values: list[Any] = []

    nnz = 0
    for new_row, old_row in enumerate(row_indices):
        start = int(indptr[old_row])
        end = int(indptr[old_row + 1])
        for entry_index in range(start, end):
            old_col = int(indices[entry_index])
            new_col = col_lookup.get(old_col)
            if new_col is None:
                continue
            sliced_indices.append(new_col)
            sliced_values.append(values[entry_index])
            nnz += 1
        sliced_indptr[new_row + 1] = nnz

    return (
        (len(row_indices), len(col_indices)),
        (
            sliced_indptr,
            np.asarray(sliced_indices, dtype=np.asarray(indices).dtype),
            np.asarray(sliced_values, dtype=np.asarray(values).dtype),
        ),
    )


def extract_block(matrix: Any, dofs: Any, row_variable: str, col_variable: str) -> Any:
    """Extract a variable-to-variable Jacobian block.

    Parameters
    ----------
    matrix:
        A dense NumPy array, raw PETSc CSR tuple, or Pandas DataFrame.
    dofs:
        DOF metadata as returned by ``load_jacobian()`` or the DataFrame returned
        by ``to_pandas()``.
    row_variable, col_variable:
        Variable names to select from the row and column spaces respectively.

    Returns
    -------
    A sliced object in the same representation family as the input ``matrix``.
    """

    row_indices = variable_indices(dofs, row_variable)
    col_indices = variable_indices(dofs, col_variable)

    if isinstance(matrix, np.ndarray):
        return matrix[np.ix_(row_indices, col_indices)]

    if isinstance(matrix, tuple) and len(matrix) == 2:
        return _slice_raw_petsc_sparse(matrix, row_indices, col_indices)

    if hasattr(matrix, "iloc"):
        return matrix.iloc[row_indices, col_indices]

    raise TypeError(
        "extract_block() supports NumPy dense arrays, raw PETSc CSR tuples, "
        "and Pandas DataFrames"
    )


def _make_label(name: str, x: int, y: int, z: int) -> str:
    return f"{name}[x={x},y={y},z={z}]"


def build_dof_table(
    index_base: np.ndarray, metadata: dict[str, Any]
) -> list[dict[str, Any]]:
    """Expand ``jacobian_index_base`` plus variable metadata into one record per DOF."""

    records: list[dict[str, Any]] = []
    variables_2d = metadata.get("variables_2d", [])
    variables_3d = metadata.get("variables_3d", [])
    n2d = int(metadata.get("n2d", len(variables_2d)))

    nx, ny, nz = index_base.shape
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                base = int(round(float(index_base[x, y, z])))
                if base < 0:
                    continue

                if z == 0:
                    for variable in variables_2d:
                        global_index = base + int(variable["offset"])
                        records.append(
                            {
                                "global_index": global_index,
                                "name": variable["name"],
                                "field_rank": "2d",
                                "offset": int(variable["offset"]),
                                "x": x,
                                "y": y,
                                "z": 0,
                                "label": _make_label(variable["name"], x, y, 0),
                            }
                        )
                    base_3d = base + n2d
                else:
                    base_3d = base

                for variable in variables_3d:
                    global_index = base_3d + int(variable["offset"])
                    records.append(
                        {
                            "global_index": global_index,
                            "name": variable["name"],
                            "field_rank": "3d",
                            "offset": int(variable["offset"]),
                            "x": x,
                            "y": y,
                            "z": z,
                            "label": _make_label(variable["name"], x, y, z),
                        }
                    )

    records.sort(key=lambda record: record["global_index"])
    return records


def load_jacobian(
    data_dir: str | Path = "data",
    matrix_filename: str | Path | None = None,
    matrix_format: str = "dense",
) -> tuple[np.ndarray | RawPetscSparse, list[dict[str, Any]], list[str]]:
    """Load the Jacobian matrix and its expanded row/column metadata.

    Returns
    -------
    matrix
        Dense ``np.ndarray`` if ``matrix_format='dense'``, otherwise the raw PETSc CSR
        tuple ``((M, N), (indptr, indices, values))``.
    dofs
        One dict per matrix row/column, sorted by ``global_index``.
    labels
        Convenience list of labels, ordered to match rows and columns of ``matrix``.
    """

    data_path = Path(data_dir)
    metadata = load_metadata(data_path / "jacobian_metadata.json")
    index_base = load_index_base(data_path)

    matrix_path = (
        Path(matrix_filename)
        if matrix_filename is not None
        else _find_matrix_file(data_path)
    )
    matrix = load_petsc_matrix(matrix_path, matrix_format=matrix_format)
    dofs = build_dof_table(index_base, metadata)
    labels = [record["label"] for record in dofs]
    return matrix, dofs, labels


def to_pandas(matrix: np.ndarray | RawPetscSparse, dofs: list[dict[str, Any]]):
    """Return ``(matrix_df, dof_df)`` using Pandas.

    Pandas is imported lazily so the core reader still works when Pandas is not installed.
    """

    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError("Pandas is required for the DataFrame adapter") from exc

    dense_matrix = to_numpy_dense(matrix)
    labels = [record["label"] for record in dofs]
    matrix_df = pd.DataFrame(dense_matrix, index=labels, columns=labels)
    dof_df = pd.DataFrame(dofs)
    return matrix_df, dof_df


def _main() -> None:
    parser = argparse.ArgumentParser(description="Read SNES Jacobian test output")
    parser.add_argument(
        "--data-dir", default="data", help="Directory containing Jacobian outputs"
    )
    parser.add_argument(
        "--matrix-format",
        choices=("dense", "sparse"),
        default="dense",
        help="Read the PETSc matrix as a dense NumPy array or raw sparse CSR tuple",
    )
    parser.add_argument(
        "--matrix-file",
        default=None,
        help="Explicit path to the PETSc binary matrix file. Defaults to the only jacobian_*.dat file",
    )
    parser.add_argument(
        "--print-matrix",
        action="store_true",
        help="Print the dense Jacobian matrix after loading",
    )
    args = parser.parse_args()

    matrix, dofs, labels = load_jacobian(
        data_dir=args.data_dir,
        matrix_filename=args.matrix_file,
        matrix_format=args.matrix_format,
    )

    if isinstance(matrix, np.ndarray):
        shape = matrix.shape
    else:
        shape = matrix[0]

    print(f"Loaded Jacobian with shape {shape}")
    print(f"Loaded {len(dofs)} row/column labels")
    print("First few DOFs:")
    for record in dofs[: min(10, len(dofs))]:
        print(f"  {record['global_index']:>4}: {record['label']}")

    if args.print_matrix:
        print(to_numpy_dense(matrix))


if __name__ == "__main__":
    _main()
