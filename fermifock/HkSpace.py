import sympy
from fermifock.FockSpace import *

from itertools import product

class HkSpace:
    """The space of even, self-adjoint k-body operators on fermionic Fock space
    of specific dimension"""

    def __init__(self, n, k):
        self.fock = FockSpace(n)
        self.n = n
        self.k = k

    def basis(self):
        """A generator for a basis of H_k"""
        for I in self.fock.basisIndicies():
            for J in self.fock.basisIndicies():
                total_length = len(I) + len(J)
                 
                if (total_length % 2 == 0) and (total_length <= 2*self.k):
                    if list(I) <= list(J):
                        yield (I,J)
                    if list(I) < list(J):
                        yield (I,J)

    def dimension(self):
        """Compute the dimension of H_k"""
        return len(list(self.basis()))

class HCkSpace:
    """The space of even, k-body operators on fermionic Fock space
    of specific dimension"""

    def __init__(self, n, k):
        self.fock = FockSpace(n)
        self.n = n
        self.k = k

    def basis(self):
        """A generator for a basis of H_k"""
        Is = self.fock.basisIndicies()

        for K, I, J in product(Is, repeat=3):
            if not self._isPairwiseDisjoint(K, I, J):
                continue

            total_length = len(I) + len(J) + 2*len(K)
            if (total_length % 2 == 0) and (total_length <= 2*self.k):
                yield NCdC(self.fock, K, I, J)

    def dimension(self):
        """Compute the dimension of H_k"""
        return len(list(self.basis()))

    def _isPairwiseDisjoint(self, *args):
        for i, j in product(range(len(args)), repeat=2):
            if (i != j) and not args[i].is_disjoint(args[j]):
                return False
        return True

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
            return sympy.sympify(2**(n - len(K.union(L).union(A).union(B))))
        else:
            return sympy.sympify(0)
