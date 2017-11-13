import sympy.physics.quantum
from sympy import Add, Mul, FiniteSet, sympify, Matrix, S
from sympy.physics.quantum import Operator, qapply, StateBase, represent

from fermifock.helpers import deduplicate_tuple, has_duplicates, signed_sort

__all__ = [
    'FermionicFockBasis',
    'FermionicFockKet',
    'FermionicFockBra',
    'Ket',
    'Bra',
    'N',
    'Cd',
    'C',
    'Vacuum',
    'represent',
    'trace',
    'qapply',
]

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

    @property
    def label_set(self):
        return FiniteSet(*self.label)

    def _eval_innerproduct_FermionicFockBra(self, bra, **hints):
        if self.label == bra.label:
            return S.One
        return S.Zero

    def _represent_FermionicFockBasis(self, basis, **kwargs):
        return Matrix(len(basis), 1, lambda i, j: (self.dual * basis[i]).doit())

    def __mul__(self, other):
        return FermionicFockKet(*list(self.label + other.label))


Ket = FermionicFockKet
Vacuum = Ket()


class FermionicFockBra(sympy.physics.quantum.Bra):
    def dual_class(self):
        return FermionicFockKet

    @classmethod
    def default_args(cls):
        return ()


Bra = FermionicFockBra


class TraceError(StandardError):
    pass


class Cd(Operator):
    def __new__(cls, *args, **kwargs):
        if has_duplicates(args):
            return S.Zero
        (sign, sorted_args) = signed_sort(args)
        return sign * Operator.__new__(cls, *sorted_args, **kwargs)

    @property
    def label_set(self):
        return FiniteSet(*self.label)

    @classmethod
    def default_args(cls):
        return ()

    def _apply_operator_FermionicFockKet(self, ket, **kwargs):
        if self.label_set.is_disjoint(ket.label_set):
            return FermionicFockKet(*self.label_set.union(ket.label_set))
        else:
            return S.Zero

    def _eval_trace(self, basis, **kwargs):
        if self.label == ():
            if basis == None:
                raise TraceError('computation of trace requires specifying a basis')
            else:
                return len(basis)
        else:
            return S.Zero

    def _represent_FermionicFockBasis(self, basis, **kwargs):
        return Matrix(len(basis), len(basis), lambda i, j: (qapply(basis[i].dual * (self * basis[j])).doit()))

    def _eval_adjoint(self):
        return C(*self.label)


class C(Operator):
    def __new__(cls, *args, **kwargs):
        if has_duplicates(args):
            return S.Zero
        (sign, sorted_args) = signed_sort(args)
        return sign * Operator.__new__(cls, *sorted_args, **kwargs)

    @classmethod
    def default_args(cls):
        return ()

    @property
    def label_set(self):
        return FiniteSet(*self.label)

    def _eval_trace(self, basis, **kwargs):
        return trace(self.adjoint(), basis, **kwargs)

    def _apply_operator_FermionicFockKet(self, ket, **kwargs):
        if self.label_set.is_subset(ket.label_set):
            K = ket.label_set - self.label_set
            return (qapply(Cd(*self.label) * Ket(*tuple(K))).dual * Ket(*ket.label)).doit() * Ket(*tuple(K))
        else:
            return S.Zero

    def _represent_FermionicFockBasis(self, basis, **kwargs):
        return Matrix(len(basis), len(basis), lambda i, j: qapply(basis[i].dual * self * basis[j]).doit())

    def _eval_adjoint(self):
        return Cd(*self.label)


class N(Operator):
    def __new__(cls, *args, **kwargs):
        return Operator.__new__(cls, *deduplicate_tuple(args), **kwargs)

    @classmethod
    def default_args(cls):
        return ()

    def _apply_operator_FermionicFockKet(self, ket, **kwargs):
        if FiniteSet(*self.label).is_subset(FiniteSet(*ket.label)):
            return ket
        return S.Zero

    def _eval_trace(self, basis, **kwargs):
        if basis == None:
            raise TraceError('computation of trace requires specifying a basis')
        if isinstance(basis, FermionicFockBasis):
            return (sympify(2) ** (-len(self.label))) * len(basis)
        else:
            return represent(self, basis=basis).trace()

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
