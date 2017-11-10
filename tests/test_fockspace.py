import unittest

from fermifock import *


class TestFermionicFockBra(unittest.TestCase):
    def test_instantiability(self):
        bra = FermionicFockBra(1, 2)

        self.assertIsInstance(bra, FermionicFockBra)
        self.assertIsInstance(bra.label, FiniteSet)

    def test_accessing_label(self):
        ket = FermionicFockBra(1)
        self.assertEqual(ket.label, FiniteSet(1))

    def test_multiplication(self):
        bra1 = FermionicFockBra(1)
        bra2 = FermionicFockBra(2)
        self.assertEqual(bra1 * bra2, FermionicFockBra(1, 2, 3))
        self.assertEqual(bra1 * bra1, 0)
