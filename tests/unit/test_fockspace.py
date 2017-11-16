from unittest import TestCase

from sympy import Matrix
from fermifock import *


class TestBra(TestCase):
    def test_instantiability(self):
        bra = Bra(1, 2)
        self.assertIsInstance(bra, Bra)

    def test_dual(self):
        self.assertEqual(Bra(1, 2).dual, Ket(1, 2))

    def test_representation(self):
        self.assertEqual(represent(Bra(1, 2), basis=FermionicFockBasis(2)), Matrix([[0, 0, 0, 1]]))


class TestKet(TestCase):
    def setUp(self):
        self.ket1 = Ket(1)
        self.ket2 = Ket(2, 3)
        self.ket3 = Ket(3, 2, 1)

    def test_instantiating_without_args_should_lead_to_empty_label(self):
        self.assertEqual(Ket().label, ())

    def test_multiplication(self):
        self.assertEqual(self.ket1 * self.ket2, Ket(1, 2, 3))
        self.assertEqual(self.ket2 * self.ket1, Ket(2, 3, 1))

    def test_innerproduct(self):
        self.assertEqual((self.ket1.dual * self.ket1).doit(), 1)
        self.assertEqual((self.ket1.dual * self.ket2).doit(), 0)

    def test_anticommutativity(self):
        self.assertEqual(self.ket1 * self.ket1, 0)
        self.assertEqual(self.ket3, -Ket(1, 2, 3))

    def test_representation(self):
        self.assertEqual(represent(self.ket1, basis=FermionicFockBasis(2)), Matrix([0, 1, 0, 0]))


class TestCd(TestCase):
    def test_antisymmetry(self):
        self.assertEqual(Cd(1, 1), 0)
        self.assertEqual(Cd(2, 1), -Cd(1, 2))

    def test_qapply(self):
        self.assertEqual(qapply(Cd() * Ket()), Ket())
        self.assertEqual(qapply(Cd(1) * Ket()), Ket(1))
        self.assertEqual(qapply(Cd(1) * Ket(1)), 0)
        self.assertEqual(qapply(Cd(1, 3, 6) * Ket(2, 4, 10)), -Ket(2, 1, 3, 4, 6, 10))

    def test_represent(self):
        B = FermionicFockBasis(2)
        self.assertEqual(represent(Cd(), basis=B), Matrix.diag(1, 1, 1, 1))
        self.assertEqual(represent(Cd(1), basis=B), Matrix([[0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0]]))

    def test_trace(self):
        self.assertEqual(trace(Cd(), FermionicFockBasis(2)), 4)
        self.assertEqual(trace(Cd(1, 2, 7)), 0)


class TestC(TestCase):
    def test_antisymmetry(self):
        self.assertEqual(C(1, 1), 0)
        self.assertEqual(C(2, 1), -C(1, 2))

    def test_qapply(self):
        self.assertEqual(qapply(C() * Vacuum), Vacuum)
        self.assertEqual(qapply(C(1) * Ket(1)), Ket())
        self.assertEqual(qapply(C(1) * Ket(2)), 0)
        self.assertEqual(qapply(C(1, 2) * Ket(1, 2, 4, 10)), Ket(4, 10))

    def test_represent(self):
        B = FermionicFockBasis(2)
        self.assertEqual(represent(C(), basis=B), Matrix.diag(1, 1, 1, 1))
        self.assertEqual(represent(C(1), basis=B), represent(Cd(1), basis=B).transpose())

    def test_trace(self):
        self.assertEqual(trace(C(), FermionicFockBasis(2)), 4)
        self.assertEqual(trace(C(1, 2, 7)), 0)


class TestN(TestCase):
    def test_single_particle_state(self):
        self.assertIsInstance(N(1), N)
        self.assertEqual(qapply(N(1) * Ket(2)), 0)
        self.assertEqual(qapply(N(1) * Ket(1)), Ket(1))

    def test_label_deduplication(self):
        self.assertEqual(N(1, 1, 2), N(1, 2))

    def test_qapply(self):
        self.assertEqual(qapply(N(1, 2) * Ket(1, 2, 3)), Ket(1, 2, 3))
        self.assertEqual(qapply(N(1, 3) * Ket(1, 2)), 0)

    def test_trace(self):
        B = FermionicFockBasis(2)
        self.assertEqual(trace(N(1, 2), basis=B), 1)
        self.assertEqual(trace(3 * N(1, 2), basis=B), 3)
        self.assertEqual(trace(N(1, 2) + N(1, 3), basis=B), 2)
        self.assertEqual(hs_scalarproduct(N(), N(2), basis=B), 2)
        self.assertEqual(hs_scalarproduct(N(1), N(1), basis=B), 2)
        self.assertEqual(hs_scalarproduct(N(1), N(2), basis=B), 1)

    def test_representation1(self):
        rep = represent(N(1), basis=FermionicFockBasis(2))
        self.assertEqual(rep, Matrix([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]))

    def test_representation2(self):
        rep = represent(N(1, 2), basis=FermionicFockBasis(2))
        self.assertEqual(rep, Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]))