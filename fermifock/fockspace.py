import sympy.physics.quantum

from sympy import FiniteSet
from sympy.physics.quantum import Operator


class FermionicFockBra(sympy.physics.quantum.Bra):
    def dual_class(self):
        return FermionicFockKet


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


class N(Operator):
    def __new__(cls, *args, **kwargs):
        return Operator.__new__(cls, *_deduplicate_tuple(args), **kwargs)

    def _apply_operator_FermionicFockKet(self, ket, **options):
        if FiniteSet(*self.label).is_subset(FiniteSet(*ket.label)):
            return ket
        return 0 * ket


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
