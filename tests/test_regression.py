"""Regression test: the numba solver reproduces the original (2023) pure-Python solver.

The reference was produced by the original code for a plunging flat plate (k = 0.1,
10 panels, 100 time steps). Run with:  python -m unittest discover tests
"""
import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class RegressionTest(unittest.TestCase):
    def test_plunging_plate_matches_original(self):
        cwd = os.getcwd()
        os.chdir(ROOT)
        try:
            from Lib.MAIN import solveVortexPanels
            from Lib.Subsolvers import LEmotion, constructBladePanels
            n, k, U = 100, 0.1, -10.0
            f = k * abs(U) / np.pi
            blades, te, wake = constructBladePanels(5 * np.pi / 180, 10, n)
            _, dT, _, _, Uu, V, Theta, _ = LEmotion(n, n / (200 * f), U, 0.0, 0.0, 0.0, 0.0, 0.05, f)
            sol = solveVortexPanels(dT, Uu, V, blades, Theta, te, wake, saveToFile=False)
        finally:
            os.chdir(cwd)
        ref = np.load(ROOT / "tests" / "data" / "reference_plunging_100steps.npz")
        np.testing.assert_allclose(sol, ref["sol"], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose([w.Xv for w in wake if w is not None], ref["wake"], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose([b.Gamma for b in blades], ref["bladeG"], rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
