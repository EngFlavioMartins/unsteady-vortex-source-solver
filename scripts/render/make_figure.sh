#!/bin/sh
# Make the README animation and figure (docs/assets/animation.webp, outputs/render/).
# Needs Python with the packages in the README and LaTeX.
set -e
cd "$(dirname "$0")/../.."
R=scripts/render
O=outputs/render
python $R/make_animation.py
python $R/compose_plate.py $O/plate-raw.png $O/plate-scale.json $O/plunging-plate
