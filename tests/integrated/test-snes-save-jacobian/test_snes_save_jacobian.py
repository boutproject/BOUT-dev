from pathlib import Path

import numpy as np
from read_jacobian import extract_block, load_jacobian, to_numpy_dense, to_pandas

TEST_DIR = Path(".")


def test_runtest(assert_success_in_shell):
    assert_success_in_shell("./test_snes_save_jacobian")

    data_dir = TEST_DIR / "data"
    metadata_file = data_dir / "jacobian_metadata.json"
    matrix_files = sorted(data_dir.glob("jacobian_*.dat"))

    assert metadata_file.exists()
    assert len(matrix_files) == 1

    dense_matrix, dofs, labels = load_jacobian(data_dir=data_dir, matrix_format="dense")
    sparse_matrix, sparse_dofs, sparse_labels = load_jacobian(
        data_dir=data_dir, matrix_format="sparse"
    )

    np.testing.assert_allclose(to_numpy_dense(sparse_matrix), dense_matrix)
    assert sparse_dofs == dofs
    assert sparse_labels == labels

    assert dense_matrix.shape == (32, 32)
    assert len(dofs) == 32
    assert labels[:4] == [
        "f[x=0,y=0,z=0]",
        "g[x=0,y=0,z=0]",
        "f[x=0,y=1,z=0]",
        "g[x=0,y=1,z=0]",
    ]

    expected_block = np.array([[-1.1, 0.0], [0.5, -2.0]])
    expected_matrix = np.kron(np.eye(16), expected_block)
    np.testing.assert_allclose(dense_matrix, expected_matrix)

    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "f", "f"), -1.1 * np.eye(16)
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "f", "g"), np.zeros((16, 16))
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "f"), 0.5 * np.eye(16)
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "g"), -2.0 * np.eye(16)
    )

    sparse_ff = extract_block(sparse_matrix, dofs, "f", "f")
    sparse_gf = extract_block(sparse_matrix, dofs, "g", "f")
    np.testing.assert_allclose(to_numpy_dense(sparse_ff), -1.1 * np.eye(16))
    np.testing.assert_allclose(to_numpy_dense(sparse_gf), 0.5 * np.eye(16))

    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None:
        matrix_df, dof_df = to_pandas(dense_matrix, dofs)
        ff_df = extract_block(matrix_df, dof_df, "f", "f")
        assert isinstance(ff_df, pd.DataFrame)
        assert list(ff_df.index[:2]) == ["f[x=0,y=0,z=0]", "f[x=0,y=1,z=0]"]
        assert list(ff_df.columns[:2]) == ["f[x=0,y=0,z=0]", "f[x=0,y=1,z=0]"]
        np.testing.assert_allclose(ff_df.to_numpy(), -1.1 * np.eye(16))


def test_rhs_jacobian(assert_success_in_shell):
    assert_success_in_shell(
        "./test_snes_save_jacobian "
        "solver:jacobian_export_kind=rhs "
        "solver:jacobian_export_prefix=jacobian_rhs_test"
    )

    data_dir = TEST_DIR / "data"
    metadata_file = data_dir / "jacobian_metadata.json"
    matrix_file = data_dir / "jacobian_rhs_test_rhs_000000.dat"

    assert metadata_file.exists()
    assert matrix_file.exists()

    dense_matrix, dofs, labels = load_jacobian(
        data_dir=data_dir, matrix_filename=matrix_file, matrix_format="dense"
    )
    sparse_matrix, sparse_dofs, sparse_labels = load_jacobian(
        data_dir=data_dir, matrix_filename=matrix_file, matrix_format="sparse"
    )

    np.testing.assert_allclose(to_numpy_dense(sparse_matrix), dense_matrix)
    assert sparse_dofs == dofs
    assert sparse_labels == labels

    assert dense_matrix.shape == (32, 32)
    assert len(dofs) == 32
    assert labels[:4] == [
        "f[x=0,y=0,z=0]",
        "g[x=0,y=0,z=0]",
        "f[x=0,y=1,z=0]",
        "g[x=0,y=1,z=0]",
    ]

    expected_block = np.array([[-0.1, 0.0], [0.5, -1.0]])
    expected_matrix = np.kron(np.eye(16), expected_block)
    np.testing.assert_allclose(dense_matrix, expected_matrix)

    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "f", "f"), -0.1 * np.eye(16)
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "f"), 0.5 * np.eye(16)
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "g"), -1.0 * np.eye(16)
    )

    sparse_gf = extract_block(sparse_matrix, dofs, "g", "f")
    np.testing.assert_allclose(to_numpy_dense(sparse_gf), 0.5 * np.eye(16))
