import copy

from sympy import Matrix, sqrt, expand, Add, Mul

from fermifock.operators import N_Operator, NCdC_Operator


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