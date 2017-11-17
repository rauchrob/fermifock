import unittest

from sympy import FiniteSet

from fermifock.operators import NCdC
from fermifock.spaces import FockSpace, HkSpace


class TestHkSpace(unittest.TestCase):
    def test_instantiable(self):
        self.assertIsInstance(HkSpace(2, 3), HkSpace)

class TestNCdC(unittest.TestCase):
    def test_scalarproduct(self):
        fock = FockSpace(5)
        A = NCdC(FiniteSet(1), FiniteSet(2), FiniteSet(3), fock)

        self.assertEqual(A.sp(A), 4)
