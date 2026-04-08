from boutdata import collect
import matplotlib.pyplot as plt
import numpy as np

path = "data"

forward_op = collect("forward_op", path=path)
forward = collect("forward", path=path)

grad_par_op = collect("grad_par_op", path=path)
grad_par_yud = collect("grad_par_yud", path=path)

div_par_op = collect("div_par_op", path=path)
div_par_yud = collect("div_par_yud", path=path)

dV = collect("dV", path=path)  # Cell volume

# Integrate divergence
op_div_int = np.sum((dV * div_par_op)[2:-2, :, :])
yud_div_int = np.sum((dV * div_par_yud)[2:-2, :, :])

print(f"Div_par volume integral Op: {op_div_int} Yup/down: {yud_div_int}")

div_par_grad_par_op = collect("div_par_grad_par_op", path=path)
div_par_grad_par_yud = collect("div_par_grad_par_yud", path=path)

# Integrate divergence
op_div2_int = np.sum((dV * div_par_grad_par_op)[2:-2, :, :])
yud_div2_int = np.sum((dV * div_par_grad_par_yud)[2:-2, :, :])

print(f"Div_par_Grad_par volume integral Op: {op_div2_int} Yup/down: {yud_div2_int}")


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
    forward_op[2:-2, 2, :], "Operator", forward[2:-2, 2, :], "Yup/down"
)
fig.suptitle("Forward")
fig.tight_layout()

fig = plot_comparison(
    grad_par_op[2:-2, 2, :], "Operator", grad_par_yud[2:-2, 2, :], "Yup/down"
)
fig.suptitle("Parallel Gradient")
fig.tight_layout()
fig.savefig("gradient.pdf")
fig.savefig("gradient.png")

fig = plot_comparison(
    div_par_op[2:-2, 2, :],
    f"Operator [integral {op_div_int:.2e}]",
    div_par_yud[2:-2, 2, :],
    f"Yup/down [integral {yud_div_int:.2e}]",
)
fig.suptitle("Parallel Divergence")
fig.tight_layout()
fig.savefig("divergence.pdf")
fig.savefig("divergence.png")

fig = plot_comparison(
    div_par_grad_par_op[2:-2, 2, :],
    f"Operator [integral {op_div2_int:.2e}]",
    div_par_grad_par_yud[2:-2, 2, :],
    f"Yup/down [integral {yud_div2_int:.2e}]",
)
fig.suptitle("Parallel diffusion")
fig.tight_layout()
fig.savefig("diffusion.pdf")
fig.savefig("diffusion.png")


plt.show()
