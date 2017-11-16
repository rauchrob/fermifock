import itertools
from sympy import Matrix, sqrt, expand, Add, Mul, Integer, sympify, Range, binomial, pprint, Expr, FiniteSet, S

from fermifock.operators import N_Operator
from fermifock.spaces import FockSpace
import copy


def gramMatrix(basis, scalarproduct):
    n = len(basis)
    return Matrix(n,n, lambda i,j: scalarproduct(basis[i],basis[j]))

def gramSchmidt(basis, scalarproduct):
    result = copy.copy(basis)
    for b_i in result:
        i = result.index(b_i)
        if i == 0:
            result[0] = b_i / sqrt(scalarproduct(b_i, b_i))
            continue

        for n_i in result[:i]:
            b_i -= (scalarproduct(b_i, n_i) / scalarProduct(n_i, n_i)) * n_i

        result[i] = expand(b_i / sqrt(scalarproduct(b_i, b_i)))
    return result


class NCdC_Operator(Expr):
    def _sympystr(self, *args, **kwargs):
        return "[%s;%s|%s]" % (self.args[0], self.args[1], self.args[2])

    def state(self):
        return self.args[0:3]

    def sp(self, other):
        if isinstance(other, NCdC_Operator):
            K, A, B = self.args[0:3]
            L, C, D = other.args[0:3]

            if (A == C) and (B == D):
                return sympify(2 ** (n - len(K.union(L).union(A).union(B))))
            else:
                return sympify(0)


NCdC = NCdC_Operator


class HCkDefaultBasis(object):
    def __init__(self, n, k):
        self.n = n
        self.k = k
        self._basis = list(iter(self))

    def __getitem__(self, item):
        return self._basis[item]

    def __len__(self):
        return len(self._basis)

    def __iter__(self):
        return iter(self._generator())

    def index(self, arg):
        return self._basis.index(arg)

    def _generator(self):
        orbits = FiniteSet(*range(1, self.n + 1))
        Is = orbits.powerset()

        for K in Is:
            for I in (orbits - K).powerset():
                for J in ((orbits - K) - I).powerset():
                    total_length = len(I) + len(J) + 2 * len(K)
                    if (total_length % 2 == 0) and (total_length <= 2 * self.k):
                        yield NCdC(K, I, J)

    def _isPairwiseDisjoint(self, *args):
        for i, j in itertools.product(range(len(args)), repeat=2):
            if (i != j) and not args[i].is_disjoint(args[j]):
                return False
        return True


def scalarProduct(A, B):
    if A.func == Add:
        return Add(*[scalarProduct(arg, B) for arg in A.args])

    if B.func == Add:
        return Add(*[scalarProduct(A, arg) for arg in B.args])

    if A.func == Mul:
        number_args = [arg for arg in A.args if arg.is_number]
        non_number_args = [arg for arg in A.args if not arg.is_number]
        return Mul(Mul(*number_args), scalarProduct(Mul(*non_number_args), B))

    if B.func == Mul:
        number_args = [arg for arg in B.args if arg.is_number]
        non_number_args = [arg for arg in B.args if not arg.is_number]
        return Mul(Mul(*number_args), scalarProduct(A, Mul(*non_number_args)))

    if A.func == N_Operator and B.func == N_Operator:
        n = A.fockspace.one_particle_dimension
        return Integer(2) ** (n - len(A.state.union(B.state)))

    if A.func == NCdC_Operator and B.func == NCdC_Operator:
        return A.sp(B)


def ONBofNk(n, k):
    F = FockSpace(n)
    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I))

    return gramSchmidt([F.N(I) for I in Is], scalarProduct)


def ONB_guess(n, k):
    F = FockSpace(n)

    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I))
    return [expand(2 ** (-(sympify(n) / 2)) * F.N_tilde2(I)) for I in Is]


def verify_ONB_guess(n, k):
    return gramMatrix(ONB_guess(n, k), scalarProduct).is_Identity


def verify_ONB_guesses():
    for n in Range(1, 7):
        for k in Range(0, n + 1):
            dim = Add(*[binomial(n, l) for l in Range(0, k + 1)])
            if dim > 100:
                continue

            if verify_ONB_guess(n, k):
                print("n=%d, k=%d (dim N_k=%d): OK" % (n, k, dim))
            else:
                print("n=%d, k=%d (dim N_k=%d): FAILED" % (n, k, dim))
                pprint(G)


def ONBofHCk(n, k):
    basis = sorted(list(HCkDefaultBasis(n, k)), key=lambda A: A.args)
    return gramSchmidt(basis, scalarProduct)


def ONBofHCk_guess(n, k):
    orbits = FiniteSet(*range(1, n + 1))
    result = []

    for K in orbits.powerset():
        for I in (orbits - K).powerset():
            for J in ((orbits - K) - I).powerset():
                basis_element = 0
                total_length = len(I) + len(J) + 2 * len(K)
                if (total_length % 2 == 0) and (total_length <= 2 * k):
                    for L in K.powerset():
                        basis_element += (sympify(-2)) ** (len(L)) * NCdC(L, I, J)

                    result.append(basis_element)

    return result


def verify_ONBofHCk_guess(n, k):
    return gramMatrix(ONBofHCk_guess(n, 2), scalarProduct).is_diagonal()


def HCk_dimension_guess(n, k):
    result = 0
    for l in range(k + 1):
        for i in range(2 * l + 1):
            result += binomial(n, i) * binomial(n, 2 * l - i)

    return result


def HCk_dimension_guess2(n, k):
    result = 0
    for l in range(k + 1):
        for k in range(l + 1):
            for i in range(2 * (l - k) + 1):
                result += binomial(n, k) * binomial(n - k, i) * binomial(n - k - i, 2 * (l - k) - i)

    return result


def ONBofHk_guess(n, k):
    orbits = FiniteSet(*range(1, n + 1))
    result = []

    for K in orbits.powerset():
        for I in (orbits - K).powerset():
            for J in ((orbits - K) - I).powerset():

                total_length = len(I) + len(J) + 2 * len(K)

                if (total_length % 2 == 0) and (total_length <= 2 * k):
                    if tuple(I) < tuple(J):
                        basis_element = 0

                        for L in K.powerset():
                            basis_element += (sympify(-2)) ** (len(L)) * S.ImaginaryUnit * (
                            NCdC(L, I, J) - NCdC(L, J, I))
                        result.append(basis_element)
                    if tuple(I) <= tuple(J):
                        basis_element = 0
                        for L in K.powerset():
                            basis_element += (sympify(-2)) ** (len(L)) * (NCdC(L, I, J) + NCdC(L, J, I))
                        result.append(basis_element)

    return result


print "Checking if our guessed OGBs for Hk are in fact OGBs and have the expected dimension:
for _n in range(7):
    for k in range(0, _n + 1):
        n = _n
        basis = ONBofHk_guess(n, k)
        dimension = len(basis)
        G = gramMatrix(basis, scalarProduct)

        print "n=%d, k=%d: dim_IR H(n,k) = %d" % (n, k, dimension)
        assert (G.is_diagonal())
        assert (G.det() != 0)
        assert (dimension == HCk_dimension_guess(n, k))

exit(0)

# Check if our guessed OGB for HCk is in fact an OGB and has the expected dimensions
# TODO: add multiprocessing capabilities?!
for i in range(0, 6):
    for k in range(0, i + 1):
        n = i

        basis = ONBofHCk_guess(n, k)
        dimension = len(basis)
        print "dim H_C(%d,%d) = %d" % (n, k, dimension)
        G = gramMatrix(basis, scalarProduct)
        assert (G.is_diagonal())
        assert (G.det() != 0)
        assert (dimension == HCk_dimension_guess(n, k))
        assert (dimension == HCk_dimension_guess2(n, k))

        if k == 1:
            assert (dimension == 2 * (n ** 2) - n + 1)

# n = 3
# basis = sorted(list(HCkDefaultBasis(n, 2)), key=lambda A: A.args)
# orthogonalized_basis = gramSchmidt(basis, scalarProduct)
#
# assert(len(basis) == len(orthogonalized_basis))
# assert(gramMatrix(orthogonalized_basis, scalarProduct).is_diagonal())
#
# for i in range(len(basis)):
#     #assert(len(orthogonalized_basis[i].args) == 2 ** len(basis[i].state()[0]))
#     print basis[i], '  |->   ', expand(orthogonalized_basis[i])
