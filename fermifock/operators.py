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


class NCdC_Operator(sympy.Expr):
    @property
    def hilbert_space(self):
        if len(self.args) >= 4:
            return self.args[3]
        return None

    def _sympystr(self, *args, **kwargs):
        return "[%s;%s|%s]" % (self.args[1], self.args[2], self.args[3])

    def states(self):
        return self.args[0:3]

    def sp(self, other):
        if isinstance(other, NCdC_Operator):
            dimension = self.hilbert_space.dimension

            K, A, B = self.states()
            L, C, D = other.states()

            if (A == C) and (B == D):
                return dimension * (sympy.sympify(2) ** (-len(K.union(L).union(A).union(B))))
            else:
                return sympy.sympify(0)


NCdC = NCdC_Operator
