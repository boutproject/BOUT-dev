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

def fac(a):
    if a==0 or a==1:
        return 1
    else:
        assert(a>0)
        return a*fac(a-1)

def gen_code(order,matrix):
    global x
    x=[]
    for i in range(order):
        x.append(Symbol("x%d"%i))
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

def taylor(i,j):
    if j >=0:
        return x[i]**j/fac(j)
    else:
        return 0

def dirichlet(i,j):
    return taylor(i,j)

def neumann(i,j):
    if i ==0:
        return taylor(i,j-1)
    else:
        return taylor(i,j)

if __name__ == "__main__":
    print(gen_code(3,dirichlet))
