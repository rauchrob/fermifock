from itertools import product

from sympy import FiniteSet, sympify, Integer

from fermifock.operators import *


class FockSpace:
    def __init__(self, n):
        self.one_particle_dimension = sympify(n)

    def basisIndicies(self):
        return FiniteSet(*range(1, self.one_particle_dimension + 1)).powerset()

    def N(self, state):
        return N_Operator(state, self)

    def N_tilde2(self, state):
        N_tilde = sympify(0)
        for K in state.powerset():
            N_tilde += Integer(-2) ** (len(K)) * self.N(K)
        return N_tilde

    def NCdC(self, K, I, J):
        return NCdC_Operator(K, I, J, self)

    @property
    def one_particle_dimension(self):
        return self.one_particle_dimension

    @property
    def dimension(self):
        return 2 ** self.one_particle_dimension


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

                if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                    if list(I) <= list(J):
                        yield (I, J)
                    if list(I) < list(J):
                        yield (I, J)

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

            total_length = len(I) + len(J) + 2 * len(K)
            if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                yield NCdC(self.fock, K, I, J)

    def dimension(self):
        """Compute the dimension of H_k"""
        return len(list(self.basis()))

    def _isPairwiseDisjoint(self, *args):
        for i, j in product(range(len(args)), repeat=2):
            if (i != j) and not args[i].is_disjoint(args[j]):
                return False
        return True