import unittest
import matplotlib.pyplot as plt
import scanpy as sc
import scez
from scez import __version__


class TestScezConfig(unittest.TestCase):
    def test_version(self):
        self.assertEqual(scez.__version__, __version__)

    def test_scanpy_settings(self):
        self.assertEqual(sc.settings.verbosity, 1)

if __name__ == '__main__':
    unittest.main()