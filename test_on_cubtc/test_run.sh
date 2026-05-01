#!/bin/bash
set -euo pipefail

PYTHON="/Users/omert/miniconda3/bin/python"
SCRIPT="/Users/omert/Desktop/MOF2frag/fragmenter.py"
DIR="/Users/omert/Desktop/MOF2frag/test_on_cubtc/"

for cif in "$DIR"/*.cif; do
  [ -e "$cif" ] || continue
  base="${cif%.cif}"

  echo "=== Running normal: $(basename "$cif") ==="
  "$PYTHON" "$SCRIPT" "$cif" --radius 6.0 --output "${base}_frag.xyz"

  echo "=== Running minimize: $(basename "$cif") ==="
  "$PYTHON" "$SCRIPT" "$cif" --radius 6.0 --minimize --output "${base}_frag_min.xyz"
done

echo "Done."

