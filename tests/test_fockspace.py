from unittest import TestCase

from sympy.physics.quantum import represent

from fermifock import *


class TestFermionicFockBra(TestCase):
    def test_instantiability(self):
        bra = FermionicFockBra(1, 2)
        self.assertIsInstance(bra, FermionicFockBra)

    def test_dual(self):
        self.assertEqual(FermionicFockBra(1, 2).dual, FermionicFockKet(1, 2))


class TestFermionicFockKet(TestCase):
    def setUp(self):
        self.ket1 = FermionicFockKet(1)
        self.ket2 = FermionicFockKet(2, 3)
        self.ket3 = FermionicFockKet(2, 1)

    def test_instantiating_without_args_should_lead_to_empty_label(self):
        self.assertEqual(FermionicFockKet().label, ())

    def test_multiplication(self):
        self.assertEqual(self.ket1 * self.ket2, FermionicFockKet(1, 2, 3))
        self.assertEqual(self.ket2 * self.ket1, FermionicFockKet(2, 3, 1))

    def test_innerproduct(self):
        self.assertEqual((self.ket1.dual * self.ket1).doit(), 1)
        self.assertEqual((self.ket1.dual * self.ket2).doit(), 0)

    def test_anticommutativity(self):
        self.assertEqual(self.ket1 * self.ket1, 0)
        self.assertEqual(self.ket3, -FermionicFockKet(1, 2))


class TestN(TestCase):
    def test_single_particle_state(self):
        self.assertIsInstance(N(1), N)
        self.assertEqual(qapply(N(1) * FermionicFockKet(2)), 0)
        self.assertEqual(qapply(N(1) * FermionicFockKet(1)), FermionicFockKet(1))

    def test_label_deduplication(self):
        self.assertEqual(N(1, 1, 2), N(1, 2))

    def test_qapply(self):
        self.assertEqual(qapply(N(1, 2) * FermionicFockKet(1, 2, 3)), FermionicFockKet(1, 2, 3))
        self.assertEqual(qapply(N(1, 3) * FermionicFockKet(1, 2)), 0)

    def test_trace(self):
        self.assertEqual(trace(N(1, 2)), 4)
        self.assertEqual(trace(3 * N(1, 2)), 12)
        self.assertEqual(trace(N(1, 2) + N(1, 3)), 8)
        self.assertRaises(NotImplementedError, trace, N(1, 2) * N(2, 3))

    # @skip('not yet implemented => TODO')
    def test_representation(self):
        rep = represent(N(1, 2), one_particle_hilbertspace_dimension=2)
        self.assertEqual(rep, Matrix([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]))
        rep = represent(N(1), one_particle_hilbertspace_dimension=2)
        self.assertEqual(rep, Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]))