import subprocess
import sys

import numpy as np
from boutdata.collect import collect, create_cache
from boututils.run_wrapper import launch_safe


def analyse_path(path, plot=False, tol=1e-12):

    cache = create_cache("data", "BOUT.dmp")

    def _collect(var):
        return collect(var, path=path, datafile_cache=cache, xguards=False, info=False)

    forward_op = _collect("forward_op")
    forward = _collect("forward")

    grad_par_op = _collect("grad_par_op")
    grad_par_yud = _collect("grad_par_yud")

    div_par_op = _collect("div_par_op")
    div_par_yud = _collect("div_par_yud")

    dV = _collect("dV")

    # Integrate divergence
    op_div_int = np.sum(dV * div_par_op)
    yud_div_int = np.sum(dV * div_par_yud)

    if np.abs(op_div_int) > tol:
        raise ValueError(f"Div_par operator has integral {op_div_int} > {tol}")

    print(f"Div_par volume integral Op: {op_div_int} Yup/down: {yud_div_int}")

    div_par_grad_par_op = _collect("div_par_grad_par_op")
    div_par_grad_par_yud = _collect("div_par_grad_par_yud")

    # Integrate divergence
    op_div2_int = np.sum(dV * div_par_grad_par_op)
    yud_div2_int = np.sum(dV * div_par_grad_par_yud)

    if np.abs(op_div2_int) > tol:
        raise ValueError(
            f"Div_par_Grad_par operator has integral {op_div2_int} > {tol}"
        )

    if plot:
        import matplotlib.pyplot as plt

        def plot_comparison(a, alabel, b, blabel):
            fig, axs = plt.subplots(1, 3, figsize=(10, 3))
            vmax = max([np.amax(a), np.amax(b)])
            vmin = min([np.amin(a), np.amin(b)])

            cs = axs[0].contourf(a, vmin=vmin, vmax=vmax)
            fig.colorbar(cs, ax=axs[0])
            axs[0].set_title(alabel)

            axs[1].contourf(b, vmin=vmin, vmax=vmax)
            axs[1].set_title(blabel)

            cs = axs[2].contourf(a - b)
            axs[2].set_title("Difference")
            fig.colorbar(cs, ax=axs[2])

            return fig

        fig = plot_comparison(
            forward_op[:, 2, :], "Operator", forward[:, 2, :], "Yup/down"
        )
        fig.suptitle("Forward")
        fig.tight_layout()

        fig = plot_comparison(
            grad_par_op[:, 2, :], "Operator", grad_par_yud[:, 2, :], "Yup/down"
        )
        fig.suptitle("Parallel Gradient")
        fig.tight_layout()
        fig.savefig("gradient.pdf")
        fig.savefig("gradient.png")

        fig = plot_comparison(
            div_par_op[:, 2, :],
            f"Operator [integral {op_div_int:.2e}]",
            div_par_yud[:, 2, :],
            f"Yup/down [integral {yud_div_int:.2e}]",
        )
        fig.suptitle("Parallel Divergence")
        fig.tight_layout()
        fig.savefig("divergence.pdf")
        fig.savefig("divergence.png")

        fig = plot_comparison(
            div_par_grad_par_op[:, 2, :],
            f"Operator [integral {op_div2_int:.2e}]",
            div_par_grad_par_yud[:, 2, :],
            f"Yup/down [integral {yud_div2_int:.2e}]",
        )
        fig.suptitle("Parallel diffusion")
        fig.tight_layout()
        fig.savefig("diffusion.pdf")
        fig.savefig("diffusion.png")


def test_petsc_operators():
    grid_path = "grid.fci.nc"

    # Generate grid file and store in temporary path
    subprocess.check_call([sys.executable, "rotating-ellipse.py", "-o", grid_path])

    nproc = 1

    _, out = launch_safe(
        "./test_petsc_operators",
        nproc=nproc,
        pipe=True,
    )

    analyse_path("data")
