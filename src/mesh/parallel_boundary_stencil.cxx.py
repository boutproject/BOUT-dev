import os
from tempfile import NamedTemporaryFile as tmpf
from stencils_sympy import dirichlet, neumann, simp, Symbol, Matrix, ccode


def gen_code(order, matrix_type):
    x = [Symbol("spacing%d" % i) for i in range(order)]
    matrix = matrix_type(x)
    A = Matrix(order, order, matrix)

    try:
        iA = A.inv()
    except:
        import sys

        print(A, matrix, file=sys.stderr)
        raise
    return ccode(simp(sum([iA[0, i] * Symbol("value%d" % i) for i in range(order)])))


def run(cmd):
    print(cmd)
    out = os.system(cmd)
    assert out == 0


if __name__ == "__main__":
    with tmpf("w", dir=".", delete=False) as f:
        f.write("namespace {\n")
        f.write(
            """
inline BoutReal pow(BoutReal val, int exp) {
  //constexpr int expval = exp;
  //static_assert(expval == 2 or expval == 3, "This pow is only for exponent 2 or 3");
  if (exp == 2) {
    return val * val;
  }
  ASSERT3(exp == 3);
  return val * val * val;
}
"""
        )

        for order in range(1, 4):
            for matrix in dirichlet, neumann:
                if order == 1 and matrix == neumann:
                    continue
                print(f"generating {matrix.name}_o{order}")
                args = ", ".join(
                    [
                        "BoutReal spacing%d, BoutReal value%d" % (i, i)
                        for i in range(order)
                    ]
                )
                f.write(
                    f"inline BoutReal stencil_{matrix.name}_o{order}({args}) {{\n  return "
                )
                f.write(gen_code(order, matrix))
                f.write(";\n}\n")
        f.write("}\n")
    run("clang-format -i " + f.name)
    run(f"mv {f.name} {__file__[:-3]}")
