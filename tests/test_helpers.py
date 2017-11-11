from unittest import TestCase

from fermifock.helpers import deduplicate_tuple, has_duplicates


class TestHelpers(TestCase):
    def test_tuple_deduplication(self):
        self.assertEqual(deduplicate_tuple((1, 3, 2)), (1, 3, 2))
        self.assertEqual(deduplicate_tuple((1, 3, 2, 3)), (1, 3, 2))

    def test_has_duplicates(self):
        self.assertEqual(has_duplicates((1, 3, 2)), False)
        self.assertEqual(has_duplicates((1, 3, 2, 3)), True)
