from sympy import Matrix, sqrt, expand, Add, Mul, Integer, sympify, Range, binomial, pprint

from fermifock.operators import N_Operator


def gramMatrix(basis, scalarproduct):
    n = len(basis)
    return Matrix(n,n, lambda i,j: scalarproduct(basis[i],basis[j]))

def gramSchmidt(basis, scalarproduct):
    for b_i in basis:
        i = basis.index(b_i)
        if i == 0:
            basis[0] = b_i / sqrt(scalarproduct(b_i, b_i))
            continue
         
        for n_i in basis[:i]:
            b_i -= scalarproduct(b_i,n_i)*n_i

        basis[i]=expand(b_i/sqrt(scalarproduct(b_i,b_i)))
    return basis


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
