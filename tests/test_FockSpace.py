from unittest import TestCase

from fermifock.FockSpace import FockSpace, scalarProduct
from sympy import FiniteSet

class TestFermionicFockSpace(TestCase):
    def setUp(self):
        self.fock = FockSpace(6)

    def test_is_instantiable(self):
        self.assertIsInstance(self.fock, FockSpace)

    def test_scalarProduct(self):
        N = self.fock.N(FiniteSet(1,2,3))
        self.assertEqual(scalarProduct(N, N), 8)
