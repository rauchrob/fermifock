from sympy import Add, Expr, FiniteSet, Integer, Matrix, Mul, sqrt, expand

class FermionicFockSpace:
    def __init__(self, n):
        self.one_particle_dimension = n

    def basisIndicies(self):
        return FiniteSet(*range(1, self.one_particle_dimension+1)).powerset()

    def N(self, state):
        return N_Operator(state, self)

    def N_tilde(self, state):
        N_tilde = sympify(0)
        for K in state.powerset():
            N_tilde += Rational(-1,2)**(len(K.complement(state)))*self.N(K)
        return N_tilde

    def N_tilde2(self, state):
        N_tilde = sympify(0)
        for K in state.powerset():
            N_tilde += Integer(-2)**(len(K))*self.N(K)
        return N_tilde

    @property
    def one_particle_dimension(self):
        return self.one_particle_dimension

class N_Operator(Expr):
    is_commutative = False
    is_number = False

    @property
    def state(self):
        return self.args[0]

    @property
    def fockspace(self):
        return self.args[1]

    def __str__(self):
        return 'N(%s)' % self.state

def scalarProduct(A,B):
    if A.func == Add:
        return Add(*[scalarProduct(arg, B) for arg in A.args])

    if B.func == Add:
        return Add(*[scalarProduct(A, arg) for arg in B.args])

    if A.func == Mul:
        number_args = [arg for arg in A.args if arg.is_number]
        non_number_args = [arg for arg in A.args if not arg.is_number]
        return Mul(Mul(*number_args), scalarProduct(Mul(*non_number_args), B))
    
    if B.func == Mul:
        number_args = [arg for arg in B.args if arg.is_number]
        non_number_args = [arg for arg in B.args if not arg.is_number]
        return Mul(Mul(*number_args), scalarProduct(A, Mul(*non_number_args)))

    if A.func == N_Operator and B.func == N_Operator:
        n = A.fockspace.one_particle_dimension
        return Integer(2)**(n-len(A.state.union(B.state)))

def gramMatrix(basis, scalarproduct):
    n = len(basis)
    return Matrix(n,n, lambda i,j: scalarproduct(basis[i],basis[j]))

def gramSchmidt(basis, scalarproduct):
    for b_i in basis:
        i = basis.index(b_i)
        if i == 0:
            basis[0] = b_i / sqrt(scalarproduct(b_i, b_i))
            continue
         
        for n_i in basis[:i]:
            b_i -= scalarproduct(b_i,n_i)*n_i

        basis[i]=expand(b_i/sqrt(scalarproduct(b_i,b_i)))
    return basis


def ONBofHk(n,k):
    F = FermionicFockSpace(n)
    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I)) 

    return gramSchmidt([F.N(I) for I in Is], scalarProduct)

from sympy import sympify, Rational, Range

def ONB_guess(n,k):
    F = FermionicFockSpace(n)

    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I)) 
    return [expand(2**(-(sympify(n)/2))*F.N_tilde2(I)) for I in Is]

from sympy import pprint
from sympy.functions.combinatorial.factorials import binomial

def verify_ONB_guess():
    for n in Range(1,7):
        for k in Range(0,n+1):
            dim = Add(*[binomial(n,l) for l in Range(0,k+1)])
            if dim > 100:
                continue

            G = gramMatrix(ONB_guess(n,k), scalarProduct)
            if G.is_Identity:
                print("n=%d, k=%d (dim N_k=%d): OK" % (n, k, dim))
            else:
                print("n=%d, k=%d (dim N_k=%d): FAILED" % (n, k, dim))
                pprint(G)

