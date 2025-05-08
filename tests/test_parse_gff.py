"""Tests all code in src.parse_gff"""
import unittest

from src.parse_gff import Feature
from src.parse_gff import FeatureCategory

class TestFeatureClass(unittest.TestCase):

    def test_Feature(self):
        """The class 'Feature' returns an object of type feature."""
        # self.assertTrue(Feature.__doc__)
        self.assertEqual(str(type(Feature("1", "contig1", FeatureCategory.Gene, 1, 100, "+", ".", parent_id=None, child_id= ["child1", "child2"]))), "<class 'src.parse_gff.Feature'>")
        self.assertIsInstance(Feature("1", "contig1", FeatureCategory.Gene, 1, 100, "+", ".", parent_id=None, child_id= ["child1", "child2"]), Feature)


    # def test_is_zero_responds_correctly_to_ints(self):
    #     """The function 'is_zero' responds correctly to integers."""
    #     self.assertTrue(is_zero(0))
    #     self.assertFalse(is_zero(1))

    # def test_is_zero_raises_an_exception_upon_non_ints(self):
    #     """The function 'is_zero' raises an exception upon non-ints."""
    #     self.assertRaises(TypeError, is_zero, {1, 2})
    #     self.assertRaises(TypeError, is_zero, "I am a string")
        