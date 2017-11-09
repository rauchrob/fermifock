import unittest

from sympy import FiniteSet

from fermifock.operators import NCdC
from fermifock.orthonormalization import gramMatrix
from fermifock.spaces import FockSpace, HkSpace, HCkSpace


class TestHkSpace(unittest.TestCase):
    def test_instantiable(self):
        self.assertIsInstance(HkSpace(2, 3), HkSpace)

    @unittest.skip('Fix dimension formula for Hk')
    def test_dimension(self):
        self.assertEqual(HkSpace(4, 1).dimension(), 2 * (4 ** 2) + 4 + 2)


class TestHCkSpace(unittest.TestCase):
    @unittest.skip('Fix dimension formula for HCk')
    def test_dimension(self):
        self.assertEqual(HCkSpace(4, 1).dimension(), 2 * (4 ** 2) + 4 + 2)

    def test_GramMatrixInvertible(self):
        basis = HCkSpace(3, 1).basis()
        G = gramMatrix(list(basis), lambda a, b: a.sp(b))

        self.assertNotEqual(G.det(), 0)


class TestNCdC(unittest.TestCase):
    def test_scalarproduct(self):
        fock = FockSpace(5)
        A = NCdC(fock, FiniteSet(1), FiniteSet(2), FiniteSet(3))

        self.assertEqual(A.sp(A), 4)
