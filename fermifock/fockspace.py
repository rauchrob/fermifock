import sympy.physics.quantum
from sympy import Add, Mul, FiniteSet, sympify, Matrix
from sympy.physics.quantum import Operator, qapply, StateBase

from fermifock.helpers import deduplicate_tuple, has_duplicates, signed_sort


class FermionicFockBasis(StateBase):
    def __init__(self, n):
        self.n = n
        l = [I for I in FiniteSet(*range(1, n + 1)).powerset()]
        l.sort(key=(lambda s: (len(s), s)))
        self.indicies = l

    def __len__(self):
        return 2 ** self.n

    def __getitem__(self, item):
        return FermionicFockKet(*tuple(self.indicies[item]))


class FermionicFockKet(sympy.physics.quantum.Ket):
    def __new__(cls, *args, **kwargs):
        if has_duplicates(args):
            return 0
        (sign, sorted_args) = signed_sort(args)
        return sign * super(FermionicFockKet, cls).__new__(cls, *sorted_args, **kwargs)

    @classmethod
    def default_args(cls):
        return ()

    def dual_class(self):
        return FermionicFockBra

    def _eval_innerproduct_FermionicFockBra(self, bra, **hints):
        if self.label == bra.label:
            return 1
        return 0

    def _represent_FermionicFockBasis(self, basis, **kwargs):
        return Matrix(len(basis), 1, lambda i, j: (self.dual * basis[i]).doit())

    def __mul__(self, other):
        return FermionicFockKet(*list(self.label + other.label))


Ket = FermionicFockKet


class FermionicFockBra(sympy.physics.quantum.Bra):
    def dual_class(self):
        return FermionicFockKet

    @classmethod
    def default_args(cls):
        return ()


Bra = FermionicFockBra


# TODO: Add C and Cd Operators

class N(Operator):
    def __new__(cls, *args, **kwargs):
        return Operator.__new__(cls, *deduplicate_tuple(args), **kwargs)

    def _apply_operator_FermionicFockKet(self, ket, **options):
        if FiniteSet(*self.label).is_subset(FiniteSet(*ket.label)):
            return ket
        return 0 * ket

    def _eval_trace(self, basis, **kwargs):
        if basis == None:
            raise ValueError('computation of trace requires specifying a basis')
        if isinstance(basis, FermionicFockBasis):
            return (sympify(2) ** (-len(self.label))) * len(basis)

    def _represent_FermionicFockBasis(self, basis, **kwargs):
        return Matrix(len(basis), len(basis), lambda i, j: (qapply(basis[i].dual * (self * basis[j])).doit()))


def trace(expr, basis=None, **kwargs):
    if expr.func == Add:
        return Add(*[trace(arg, basis=basis) for arg in expr.args])
    if expr.func == Mul:
        number_args = [arg for arg in expr.args if arg.is_Number]
        non_number_args = [arg for arg in expr.args if not arg.is_Number]
        if len(non_number_args) == 1:
            return Mul(Mul(*number_args), trace(non_number_args[0], basis=basis))

    if expr.func.is_Number:
        return expr

    if hasattr(expr.func, '_eval_trace'):
        return expr._eval_trace(basis, **kwargs)

    raise NotImplementedError
