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
    #print(A)
    try:
        iA=A.inv()
    except:
        import sys
        print(A,matrix,file=sys.stderr)
        raise
    #for i in range(order):
    #    for j in range(order):
    #        print(i,j, simp(iA[i,j]))
    #print(iA)
    #r=Matrix(order,1, lambda i,j: 1 if i==0 else 0)
    #print(r)
    #s=iA*r
    ret=""
    for i in range(order):
        ret+=ccode(simp(iA[0,i]),assign_to="fac%d"%i)
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
