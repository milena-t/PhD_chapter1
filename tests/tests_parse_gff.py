"""Tests all code in src.weather.Milena_utils."""
import unittest

from src.parse_gff import Feature

class TestFeatureClass(unittest.TestCase):

    """Class to test the code in src.weather.Milena_utils."""

    def test_Feature_has_documentation(self):
        """The class 'Feature' has documentation."""
        self.assertTrue(Feature.__doc__)

    # def test_is_zero_has_filled_documentation(self):
    #     """The function 'is_zero' has documentation."""
    #     self.assertIsNotNone(is_zero.__doc__)


    # def test_is_zero_responds_correctly_to_ints(self):
    #     """The function 'is_zero' responds correctly to integers."""
    #     self.assertTrue(is_zero(0))
    #     self.assertFalse(is_zero(1))

    # def test_is_zero_raises_an_exception_upon_non_ints(self):
    #     """The function 'is_zero' raises an exception upon non-ints."""
    #     self.assertRaises(TypeError, is_zero, {1, 2})
    #     self.assertRaises(TypeError, is_zero, "I am a string")
        