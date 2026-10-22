"""Animation for the README: the vortex street behind a plunging flat plate (k = 1).

Reduced frequency k = 1, plunge amplitude 0.2 chord, zero mean angle of attack, 10 panels:
the wake rolls up into a reverse von Karman street (thrust-producing plunge).

Runs the solver step by step (the same calls as Lib.MAIN.solveVortexPanels) and draws each
wake vortex as a Gaussian blob of its circulation, in a frame that follows the plate along the stream (so the plate heaves), on black:
teal and red are the two senses of rotation and zero fades into the black.
Writes docs/assets/animation.webp, and the last frame with its colour scale to outputs/render/
for compose_plate.py (the still figure).  Run:  python scripts/render/make_animation.py
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT), str(HERE)]
os.chdir(ROOT)
from gif import save_gif  # noqa: E402
from onda_style import colour  # noqa: E402
from Lib.Classes import wakeVortex  # noqa: E402
from Lib.MAIN import solveLinearSystem  # noqa: E402
from Lib.Subsolvers import (LEmotion, constructBladePanels, initializeInductionMatrix, moveBlade,  # noqa: E402
                            updateInductionMatrix, vortexRollUp, wakeInducedVelocity)

K, U, H_PLUNGE, NSTEPS, DT, EVERY = 1.0, -10.0, 0.2, 520, 0.005, 3
W, H = 1260, 330                       # pixels; view 7.6 x 2.0 chords behind the leading edge
X0, X1, Y0, Y1 = -0.6, 7.0, -1.0, 1.0
SIGMA = 0.07                           # blob radius [chords], about 1.4 x the vortex spacing U dt

f = K * abs(U) / np.pi
blades, te, wake = constructBladePanels(0.0, 10, NSTEPS)
Time, dT, X, Y, Uu, V, Theta, Om = LEmotion(NSTEPS, (NSTEPS - 1) * DT, U, 0.0, 0.0, 0.0, 0.0, H_PLUNGE, f)
A = initializeInductionMatrix(blades, te)

lut = np.array([colour("signed_dark", t) for t in np.linspace(0, 1, 256)])
xs, ys = np.linspace(X0, X1, W), np.linspace(Y1, Y0, H)
gx, gy = np.meshgrid(xs, ys)
states = []
for t, dt in enumerate(dT):
    moveBlade(t, dt, Uu, V, blades, Theta, te)
    updateInductionMatrix(blades, te, A)
    solveLinearSystem(A, wakeInducedVelocity(blades, wake), blades, te)
    vortexRollUp(blades, te, wake, dt)
    wake[t] = wakeVortex(te.Xv)
    wake[t].Gamma = te.Gamma
    if t % EVERY == 0 and t > 0:
        le = np.array([blades[0].Xa[0], 0.0])      # follow the plate along the stream only, so it heaves
        states.append((np.array([w.Xv - le for w in wake if w is not None]),
                       np.array([w.Gamma for w in wake if w is not None]),
                       blades[0].Xa - le, blades[-1].Xb - le))


def field(pts, gam):
    keep = (pts[:, 0] > X0 - 0.4) & (pts[:, 0] < X1 + 0.4) & (np.abs(pts[:, 1]) < 1.5)
    out = np.zeros((H, W))
    for (px, py), g in zip(pts[keep], gam[keep]):
        out += g * np.exp(-((gx - px) ** 2 + (gy - py) ** 2) / (2 * SIGMA ** 2))
    return out


# one colour scale for the whole animation, from the developed street (last third)
scale = np.percentile(np.abs(np.concatenate([field(*s[:2]).ravel() for s in states[-len(states) // 3::8]])), 99.5)
frames = []
for pts, gam, lead, tip in states:
    img = lut[np.clip(((np.tanh(field(pts, gam) / scale) + 1) / 2 * 255).astype(int), 0, 255)]
    for s_ in np.linspace(0, 1, 400):                       # the plate: a light-grey bar
        cx, cy = lead[0] + s_ * (tip[0] - lead[0]), lead[1] + s_ * (tip[1] - lead[1])
        j, i = int((cx - X0) / (X1 - X0) * (W - 1)), int((Y1 - cy) / (Y1 - Y0) * (H - 1))
        img[max(i - 3, 0):i + 4, max(j - 1, 0):j + 2] = 0.72
    frames.append((img * 255).astype(np.uint8))
save_gif(frames, ROOT / "docs" / "assets" / "animation.webp", fps=20, width=840)

# the still: last frame, and the colour scale in vorticity units (a blob of circulation G has
# peak vorticity G / (2 pi SIGMA^2); omega c / |U| with the chord c = 1)
render = ROOT / "outputs" / "render"
render.mkdir(parents=True, exist_ok=True)
Image.fromarray(frames[-1]).save(render / "plate-raw.png")
(render / "plate-scale.json").write_text(json.dumps({"omega_scale": scale / (2 * np.pi * SIGMA**2 * abs(U))}))
