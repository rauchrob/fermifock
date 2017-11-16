import copy

from sympy import Matrix, sqrt, expand, Add, Mul, Range, binomial, pprint

from fermifock.operators import N_Operator, NCdC_Operator
from fermifock.spaces import FockSpace


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
            b_i -= (scalarproduct(b_i, n_i) / scalarproduct(n_i, n_i)) * n_i

        result[i] = expand(b_i / sqrt(scalarproduct(b_i, b_i)))
    return result

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
        return A.sp(B)

    if A.func == NCdC_Operator and B.func == NCdC_Operator:
        return A.sp(B)


#################################################
# TODO: move remaining parts into tests/theorems

def ONBofNk(n, k):
    F = FockSpace(n)
    Is = sorted([I for I in F.basisIndicies() if len(I) <= k], key=lambda I: (len(I), I))

    return gramSchmidt([F.N(I) for I in Is], scalarProduct)


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


from fermifock.spaces import HkSpace, HCkSpace

print "Checking if our guessed OGBs for Hk are in fact OGBs and have the expected dimension"
for n in range(4):
    for k in range(0, n + 1):
        hkspace = HkSpace(n, k)
        basis = hkspace.basis('orthogonal')
        actual_dimension = len(basis)
        G = gramMatrix(basis, scalarProduct)

        print "n=%d, k=%d: dim_IR H(n,k) = %d" % (n, k, actual_dimension)
        assert (G.is_diagonal())
        assert (G.det() != 0)
        assert (actual_dimension == hkspace.dimension)

print "Check if our guessed OGB for HCk is in fact an OGB and has the expected dimensions"
for n in range(0, 6):
    for k in range(0, n + 1):
        hckspace = HCkSpace(n, k)
        basis = hckspace.basis('orthogonal')
        actual_dimension = len(basis)
        print "dim H_C(%d,%d) = %d" % (n, k, actual_dimension)
        G = gramMatrix(basis, scalarProduct)
        assert (G.is_diagonal())
        assert (G.det() != 0)
        assert (actual_dimension == hckspace.dimension)

        if k == 1:
            assert (actual_dimension == 2 * (n ** 2) - n + 1)