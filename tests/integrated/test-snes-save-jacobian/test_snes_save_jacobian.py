from pathlib import Path
import numpy as np
from read_jacobian import load_jacobian, to_numpy_dense

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
