from unittest import TestCase

import fermifock.FockSpace

class TestFermionicFockSpace(TestCase):
    def test_is_instantiable(self):
        self.assertIsInstance(fermifock.FockSpace(3), fermifock.FockSpace)
