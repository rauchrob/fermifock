from sympy import Add, Expr, FiniteSet, Integer, Matrix, Mul, sqrt, expand

class FockSpace:
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
