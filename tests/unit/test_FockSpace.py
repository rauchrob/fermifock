from unittest import TestCase

from sympy import FiniteSet

from fermifock.orthonormalization import scalarProduct
from fermifock.spaces import FockSpace


class TestFockSpace(TestCase):
    def setUp(self):
        self.fock = FockSpace(6)

    def test_is_instantiable(self):
        self.assertIsInstance(self.fock, FockSpace)

    def test_scalarProduct(self):
        N = self.fock.N(FiniteSet(1, 2, 3))
        self.assertEqual(scalarProduct(N, N), 8)
