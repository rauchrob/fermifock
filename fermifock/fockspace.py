from sympy import FiniteSet
from sympy.physics.quantum import Bra, Ket, Operator


class FermionicFockBra(Bra):
    def __new__(cls, *args, **kwargs):
        # if (len(args) == 1):
        #     if isinstance(args[0], FiniteSet):
        #         return Bra.__new__(cls, args, **kwargs)
        #     return Bra.__new__(cls, args, **kwargs)
        if len(args) == 1 and isinstance(args[0], FiniteSet):
            return Bra.__new__(cls, *args, **kwargs)
        if len(args) >= 1:
            return Bra.__new__(cls, FiniteSet(*list(args)), **kwargs)

    def dual_class(self):
        return FermionicFockKet

    @property
    def label(self):
        return super(FermionicFockBra, self).label[0]

    def __mul__(self, other):
        if self.label.is_disjoint(other.label):
            return FermionicFockBra(self.label.union(other.label))
        return 0 * self


class FermionicFockKet(Ket):
    # TODO: ensure label is a single FiniteSet
    def dual_class(self):
        return FermionicFockBra

    def _eval_innerproduct_FermionicFockBra(self, bra, **hints):
        if self.label == bra.label:
            return 1
        return 0


class TraceClassOperator(Operator):
    def trace(self):
        raise NotImplementedError


class N(TraceClassOperator):
    def _apply_operator_FermionicFockKet(self, ket, **options):
        if self.label[0].is_subset(ket.label[0]):
            return ket
        return 0 * ket

        # print(FermionicFockBra(FiniteSet(1, 2)) * FermionicFockKet(FiniteSet(1, 3)))

        # A = N(FiniteSet(1,2))
        # ket = FermionicFockKet(FiniteSet(1,2))
        # print qapply((A+A)*ket)
        # print qapply((A)*ket)
