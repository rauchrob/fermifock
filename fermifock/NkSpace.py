from sympy import sympify, Range, Add, expand
from fermifock.FockSpace import FockSpace, scalarProduct
from fermifock.GramSchmidt import gramMatrix, gramSchmidt

def ONBofNk(n,k):
    F = FockSpace(n)
    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I)) 

    return gramSchmidt([F.N(I) for I in Is], scalarProduct)

def ONB_guess(n,k):
    F = FockSpace(n)

    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I)) 
    return [expand(2**(-(sympify(n)/2))*F.N_tilde2(I)) for I in Is]

def verify_ONB_guess(n, k):
    return gramMatrix(ONB_guess(n,k), scalarProduct).is_Identity

from sympy import pprint
from sympy.functions.combinatorial.factorials import binomial

def verify_ONB_guesses():
    for n in Range(1,7):
        for k in Range(0,n+1):
            dim = Add(*[binomial(n,l) for l in Range(0,k+1)])
            if dim > 100:
                continue

            if verify_ONB_guess(n,k):
                print("n=%d, k=%d (dim N_k=%d): OK" % (n, k, dim))
            else:
                print("n=%d, k=%d (dim N_k=%d): FAILED" % (n, k, dim))
                pprint(G)

