#!/usr/bin/env python3

from sympy import Symbol, Eq
from sympy.matrices import Matrix
from sympy.printing import ccode
from sympy.simplify import combsimp as simp
from sympy.utilities.codegen import codegen


def pow(a,b):
    if b == 0:
        return "1"
    if b == 1:
        return a
    else:
        return "%s**%d"%(a,b)

def factorial(a):
    if a==0 or a==1:
        return 1
    else:
        assert(a>0)
        return a*factorial(a-1)

def gen_code(order,matrix):
    x=[]
    for i in range(order):
        x.append(Symbol("spacing.f%d"%i))

    matrix = matrix(x)
    A=Matrix(order,order,matrix)

    try:
        iA=A.inv()
    except:
        import sys
        print(A,matrix,file=sys.stderr)
        raise
    ret=""
    for i in range(order):
        ret+=ccode(simp(iA[0,i]),assign_to="facs.f%d"%i)
        ret+="\n"
    return ret

def taylor(x, i, j):
    if j >=0:
        return x[i]**j/factorial(j)
    else:
        return 0

class dirichlet:
    def __init__(self, x):
        self.x = x
    def __call__(self, i, j):
        return taylor(self.x, i ,j)

class neumann:
    def __init__(self, x):
        self.x = x
    def __call__(self, i, j):
        if i ==0:
            return taylor(self.x, i, j-1)
        else:
            return taylor(self.x, i, j)

if __name__ == "__main__":
    print(gen_code(3,dirichlet))
