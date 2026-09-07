from pathlib import Path

import numpy as np
from boutdata.collect import collect
from read_jacobian import extract_block, load_jacobian, to_numpy_dense

TEST_DIR = Path(".")


def test_rhs_jacobian_output_trigger(assert_success_in_shell):
    assert_success_in_shell("./test_cvode_save_jacobian")

    data_dir = TEST_DIR / "data"
    metadata_file = data_dir / "jacobian_metadata.json"
    matrix_file = data_dir / "jacobian_rhs_000000.dat"

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


def test_system_jacobian_linear_setup_trigger(assert_success_in_shell):
    assert_success_in_shell(
        "./test_cvode_save_jacobian "
        "solver:cvode_precon_method=petsc "
        "solver:jacobian_export_kind=system "
        "solver:jacobian_export_trigger=linear_setup "
        "solver:jacobian_export_prefix=jacobian_system_test"
    )

    data_dir = TEST_DIR / "data"
    metadata_file = data_dir / "jacobian_metadata.json"
    matrix_file = data_dir / "jacobian_system_test_system_000000.dat"

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

    last_step = float(
        np.ravel(
            np.asarray(collect("cvode_last_step", path=str(data_dir), info=False))
        )[-1]
    )
    expected_block = np.array(
        [[1.0 + 0.1 * last_step, 0.0], [-0.5 * last_step, 1.0 + last_step]]
    )
    expected_matrix = np.kron(np.eye(16), expected_block)
    np.testing.assert_allclose(dense_matrix, expected_matrix)

    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "f", "f"),
        (1.0 + 0.1 * last_step) * np.eye(16),
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "f"),
        (-0.5 * last_step) * np.eye(16),
    )
    np.testing.assert_allclose(
        extract_block(dense_matrix, dofs, "g", "g"),
        (1.0 + last_step) * np.eye(16),
    )
