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

    def test_figure_params(self):
        self.assertEqual(sc.settings.set_figure_params, {'dpi': 100, 'dpi_save': 300, 'frameon': False, 'figsize': (5, 5), 'facecolor': 'white'})

    def test_max_open_warning(self):
        self.assertEqual(plt.rcParams['figure.max_open_warning'], 0)

if __name__ == '__main__':
    unittest.main()