import unittest

from sympy import FiniteSet

from fermifock.orthonormalization import gramMatrix, scalarProduct
from fermifock.spaces import FockSpace, HkSpace, HCkSpace, NkSpace


class TestNkSpaceDefaultBasis(unittest.TestCase):
    def setUp(self):
        self.space = NkSpace(3, 2)
        self.basis = self.space.basis("default")
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_actual_dimension_matches_expected_dimension(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_offdiagonal(self):
        self.assertFalse(self.grammatrix.is_diagonal())

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())


class TestNkSpaceOrthogonalBasis(unittest.TestCase):
    def setUp(self):
        self.space = NkSpace(3, 2)
        self.basis = self.space.basis("orthogonal")
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_actual_dimension_matches_expected_dimension(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_diagonal(self):
        self.assertTrue(self.grammatrix.is_diagonal())


class TestHCkSpaceDefaultBasis(unittest.TestCase):
    def setUp(self):
        self.space = HCkSpace(3, 2)
        self.basis = self.space.basis("default")
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_actual_dimension_matches_expected_dimension(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_offdiagonal(self):
        self.assertFalse(self.grammatrix.is_diagonal())

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())


class TestHCkSpaceOrthogonalBasis(unittest.TestCase):
    def setUp(self):
        self.space = HCkSpace(3, 2)
        self.basis = self.space.basis('orthogonal')
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_basis_has_correct_length(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_diagonal(self):
        self.assertTrue(self.grammatrix.is_diagonal())


class TestHkSpaceDefaultBasis(unittest.TestCase):
    def setUp(self):
        self.space = HkSpace(3, 2)
        self.basis = self.space.basis("default")
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_actual_dimension_matches_expected_dimension(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_offdiagonal(self):
        self.assertFalse(self.grammatrix.is_diagonal())

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())


class TestHkSpaceOrthogonalBasis(unittest.TestCase):
    def setUp(self):
        self.space = HkSpace(3, 2)
        self.basis = self.space.basis('orthogonal')
        self.grammatrix = gramMatrix(self.basis, scalarProduct)

    def test_basis_has_correct_length(self):
        self.assertEqual(self.space.dimension, len(self.basis))

    def test_gram_matrix_is_symmetric(self):
        self.assertTrue(self.grammatrix.is_symmetric())

    def test_gram_matrix_is_invertible(self):
        self.assertNotEqual(self.grammatrix.det(), 0)

    def test_gram_matrix_is_diagonal(self):
        self.assertTrue(self.grammatrix.is_diagonal())


class TestNCdC(unittest.TestCase):
    def test_scalarproduct(self):
        A = FockSpace(5).NCdC(FiniteSet(1), FiniteSet(2), FiniteSet(3))

        self.assertEqual(scalarProduct(A, A), 4)
