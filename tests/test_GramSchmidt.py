from unittest import TestCase

from fermifock.GramSchmidt import gramMatrix, gramSchmidt
from sympy import Matrix

class GramSchmidtWithOrthonormalBasis(TestCase):
    def setUp(self):
        self.basis = [1,2,3,4]
        self.scalarProduct = lambda i,j: 1 if i == j else 0

    def test_GramMatrixIsTrivial(self):
        gram = gramMatrix(self.basis, self.scalarProduct)
        self.assertEqual(gram, Matrix.diag(1,1,1,1))

    def test_GramSchmidtReturnsInputBasis(self):
        self.assertEqual(gramSchmidt(self.basis, self.scalarProduct), self.basis)
