"""Add the LaTeX colour bar to the plunging-plate render.

Usage: python compose_plate.py <raw png> <scale json> <output stem>
The bar follows the tanh compression of make_animation.py, so its ticks sit at their true
vorticity values (omega c / U).
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlay import compose  # noqa: E402

raw, scale_json, stem = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
scale = json.loads(scale_json.read_text())["omega_scale"]
tick = float(f"{1.5 * scale:.1g}")                     # one round value near the saturation
frac = lambda w: 0.5 + 0.5 * math.tanh(w / scale)      # noqa: E731
label = f"{tick:g}"
compose([(raw, None)], stem, bar=dict(name="signed_dark", vmin=-1, vmax=1, title=r"$\omega c/U$",
                                     ticks=[(frac(-tick), "$-$" + label), (frac(0), "0"), (frac(tick), label)]))
