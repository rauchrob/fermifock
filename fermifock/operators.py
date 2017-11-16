import sympy


class N_Operator(sympy.Expr):
    is_commutative = False
    is_number = False

    @property
    def state(self):
        return self.args[0]

    @property
    def fockspace(self):
        return self.args[1]

    def _sympystr(self, *args, **kwargs):
        return 'N(%s)' % ','.join(map(str, sorted(self.state)))


N = N_Operator


class NCdC(sympy.Expr):
    is_commutative = False
    is_number = False

    @property
    def state(self):
        return (self.args[1], self.args[2], self.args[3])

    @property
    def fockspace(self):
        return self.args[0]

    def __str__(self):
        return 'NCdC(%s)' % self.state

    def sp(self, other):
        if not isinstance(other, NCdC):
            raise NotImplementedError

        K, A, B = self.state
        L, C, D = other.state

        if (A == C) and (B == D):
            n = self.fockspace.one_particle_dimension
            return sympy.sympify(2 ** (n - len(K.union(L).union(A).union(B))))
        else:
            return sympy.sympify(0)

    def _sympystr(self, *args, **kwargs):
        return 'N(%s)' % ','.join(map(str, sorted(self.state)))
