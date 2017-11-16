from unittest import TestCase

from fermifock.helpers import deduplicate_tuple, has_duplicates, signed_sort


class TestHelpers(TestCase):
    def test_tuple_deduplication(self):
        self.assertEqual(deduplicate_tuple((1, 3, 2)), (1, 3, 2))
        self.assertEqual(deduplicate_tuple((1, 3, 2, 3)), (1, 3, 2))

    def test_has_duplicates(self):
        self.assertEqual(has_duplicates((1, 3, 2)), False)
        self.assertEqual(has_duplicates((1, 3, 2, 3)), True)

    def test_signed_sort(self):
        self.assertEqual(signed_sort((1, 2, 4, 5)), (1, (1, 2, 4, 5)))
        self.assertEqual(signed_sort((1, 2, 4, 3)), (-1, (1, 2, 3, 4)))
        self.assertEqual(signed_sort((4, 3, 2, 1)), (1, (1, 2, 3, 4)))
        self.assertEqual(signed_sort((3, 2, 1)), (-1, (1, 2, 3)))
