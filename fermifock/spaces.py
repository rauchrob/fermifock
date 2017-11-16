from itertools import product

from sympy import FiniteSet, sympify, Integer, S, binomial

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


class NkSpace:
    def __init__(self, n, k):
        self.fock = FockSpace(n)
        self.n = n
        self.k = k

    def basis(self, name):
        if name == 'default':
            return [self.fock.N(I) for I in FiniteSet(*range(n + 1)).powerset() if len(I) <= self.k]
        if name == 'orthogonal':
            return [self.fock.N_tilde2(I) for I in FiniteSet(*range(n + 1)).powerset() if len(I <= self.k)]

    def dimension(self):
        return binomial(self.n, self.k)


class HCkSpace:
    """The space of even, k-body operators on fermionic Fock space
    of specific dimension"""

    def __init__(self, n, k):
        self.fock = FockSpace(n)
        self.n = n
        self.k = k

    def basis(self, name):
        if name == 'default':
            orbits = FiniteSet(*range(1, self.n + 1))
            result = []

            for K in orbits.powerset():
                for I in (orbits - K).powerset():
                    for J in ((orbits - K) - I).powerset():
                        basis_element = 0
                        total_length = len(I) + len(J) + 2 * len(K)
                        if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                            result.append(self.fock.NCdC(K, I, J))
            return result

        if name == 'orthogonal':
            orbits = FiniteSet(*range(1, self.n + 1))
            result = []

            for K in orbits.powerset():
                for I in (orbits - K).powerset():
                    for J in ((orbits - K) - I).powerset():
                        basis_element = 0
                        total_length = len(I) + len(J) + 2 * len(K)
                        if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                            for L in K.powerset():
                                basis_element += (sympify(-2)) ** (len(L)) * self.fock.NCdC(L, I, J)

                            result.append(basis_element)

            return result

    @property
    def dimension(self):
        result = 0
        for l in range(self.k + 1):
            for i in range(2 * l + 1):
                result += binomial(self.n, i) * binomial(self.n, 2 * l - i)

        return result


class HkSpace:
    """The space of even, self-adjoint k-body operators on fermionic Fock space
    of specific dimension"""

    def __init__(self, n, k):
        self.fock = FockSpace(n)
        self.n = n
        self.k = k

    @property
    def dimension(self):
        return HCkSpace(self.n, self.k).dimension

    def basis(self, name='default'):
        if name == 'default':
            result = []
            orbits = FiniteSet(*range(1, self.n + 1))
            Is = orbits.powerset()

            for K in Is:
                for I in (orbits - K).powerset():
                    for J in ((orbits - K) - I).powerset():
                        total_length = len(I) + len(J) + 2 * len(K)
                        if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                            result.append(self.fock.NCdC(K, I, J))

            return result
        if name == 'orthogonal':
            orbits = FiniteSet(*range(1, self.n + 1))
            result = []

            for K in orbits.powerset():
                for I in (orbits - K).powerset():
                    for J in ((orbits - K) - I).powerset():

                        total_length = len(I) + len(J) + 2 * len(K)

                        if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                            if tuple(I) < tuple(J):
                                basis_element = 0

                                for L in K.powerset():
                                    basis_element += (sympify(-2)) ** (len(L)) * S.ImaginaryUnit * (
                                        self.fock.NCdC(L, I, J) - self.fock.NCdC(L, J, I))
                                result.append(basis_element)
                            if tuple(I) <= tuple(J):
                                basis_element = 0
                                for L in K.powerset():
                                    basis_element += (sympify(-2)) ** (len(L)) * (
                                        self.fock.NCdC(L, I, J) + self.fock.NCdC(L, J, I))
                                result.append(basis_element)

            return result
