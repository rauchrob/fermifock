import sympy.physics.quantum

from sympy import Add, Mul, FiniteSet, sympify, Matrix
from sympy.physics.quantum import Operator, qapply


class FermionicFockBra(sympy.physics.quantum.Bra):
    def dual_class(self):
        return FermionicFockKet

    @classmethod
    def default_args(cls):
        return ()

Bra = FermionicFockBra

class FermionicFockKet(sympy.physics.quantum.Ket):
    def __new__(cls, *args, **kwargs):
        if _has_duplicates(args):
            return 0
        return super(FermionicFockKet, cls).__new__(cls, *args, **kwargs)

    def dual_class(self):
        return FermionicFockBra

    def _eval_innerproduct_FermionicFockBra(self, bra, **hints):
        if self.label == bra.label:
            return 1
        return 0

    def __mul__(self, other):
        return FermionicFockKet(*list(self.label + other.label))

    @classmethod
    def default_args(cls):
        return ()

Ket = FermionicFockKet

class N(Operator):
    def __new__(cls, *args, **kwargs):
        return Operator.__new__(cls, *_deduplicate_tuple(args), **kwargs)

    def _apply_operator_FermionicFockKet(self, ket, **options):
        if FiniteSet(*self.label).is_subset(FiniteSet(*ket.label)):
            return ket
        return 0 * ket

    def _eval_trace(self, **hints):
        return sympify(2) ** len(self.label)

    def _represent_default_basis(self, **options):
        n = options['one_particle_hilbertspace_dimension']
        basis = [FermionicFockKet(*tuple(s)) for s in FiniteSet(*range(1, n + 1)).powerset()]

        result = Matrix.zeros(2 ** n)
        return Matrix(2 ** n, 2 ** n, lambda i, j: (qapply(basis[i].dual * (self * basis[j])).doit()))

def trace(expr, **hints):
    if expr.func == Add:
        return Add(*[trace(arg) for arg in expr.args])
    if expr.func == Mul:
        number_args = [arg for arg in expr.args if arg.is_Number]
        non_number_args = [arg for arg in expr.args if not arg.is_Number]
        if len(non_number_args) == 1:
            return Mul(Mul(*number_args), trace(non_number_args[0]))

    if expr.func.is_Number:
        return expr

    if hasattr(expr.func, '_eval_trace'):
        return expr._eval_trace(**hints)

    raise NotImplementedError

def _deduplicate_tuple(t):
    seen_entries = []
    result = ()
    for entry in t:
        if entry in seen_entries:
            continue
        seen_entries.append(entry)
        result += (entry,)
    return result


def _has_duplicates(t):
    seen_entries = []
    for entry in t:
        if entry in seen_entries:
            return True
        seen_entries.append(entry)
    return False
