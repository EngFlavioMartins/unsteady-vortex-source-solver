"""Plunging flat plate with the unsteady vortex panel method (MAIN.ipynb, reduced frequency k = 0.1).

10 panels, angle of attack 5 degrees, plunge amplitude 0.05 chord, 445 time steps. Prints the
run time and saves the bound-circulation history and the final wake to outputs/plunging_plate.png.
"""
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
from Lib.MAIN import solveVortexPanels  # noqa: E402
from Lib.Subsolvers import LEmotion, constructBladePanels  # noqa: E402

k, U, n = 0.1, -10.0, 445
f = k * abs(U) / np.pi
blades, te, wake = constructBladePanels(5 * np.pi / 180, 10, n)
Time, dT, X, Y, Uu, V, Theta, Om = LEmotion(n, n / (200 * f), U, 0.0, 0.0, 0.0, 0.0, 0.05, f)

t0 = time.time()
circulation = solveVortexPanels(dT, Uu, V, blades, Theta, te, wake, saveToFile=False)
print(f"Solved {n} time steps in {time.time() - t0:.2f} s")

out = ROOT / "outputs"
out.mkdir(exist_ok=True)
wx = np.array([w.Xv for w in wake if w is not None])
wg = np.array([w.Gamma for w in wake if w is not None])
le = blades[0].Xa
fig, axs = plt.subplots(2, 1, figsize=(8, 5), constrained_layout=True)
axs[0].plot(Time[:-1], circulation, lw=1.2)
axs[0].set(xlabel="time [s]", ylabel=r"bound circulation $\Gamma$", title=f"Plunging plate, k = {k}")
lim = np.max(np.abs(wg))
axs[1].scatter(wx[:, 0] - le[0], wx[:, 1] - le[1], c=wg, s=4, cmap="RdBu_r", vmin=-lim, vmax=lim)
axs[1].plot([0, blades[-1].Xb[0] - le[0]], [0, blades[-1].Xb[1] - le[1]], "k-", lw=2)
axs[1].set(xlabel="$x/c$", ylabel="$y/c$", title="Wake vortices (colour: circulation)", aspect="equal")
fig.savefig(out / "plunging_plate.png", dpi=150)
print(f"Figure: {out / 'plunging_plate.png'}")
