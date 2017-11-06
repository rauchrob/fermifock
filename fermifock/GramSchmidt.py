from sympy import Matrix, sqrt, expand

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
