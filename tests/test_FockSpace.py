from unittest import TestCase

from fermifock.FockSpace import FockSpace, scalarProduct
from sympy import FiniteSet

class TestFermionicFockSpace(TestCase):
    def test_is_instantiable(self):
        self.assertIsInstance(FockSpace(3), FockSpace)

    def test_scalarProduct(self):
        fock = FockSpace(6)
        K = FiniteSet(1,2,3)

        self.assertEqual(scalarProduct(fock.N(K), fock.N(K)), 8)
