#!/usr/bin/env bash
set -euo pipefail

# Simple reproducible verification runner.
# Usage:
#   bash scripts/run_verification.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python -V
pip -V

echo "[1/3] Installing package (editable)"
pip install --no-build-isolation -e .

echo "[2/3] Running unit tests"
pytest -q

echo "[3/3] Capturing git and dependency snapshot"
git rev-parse HEAD
pip freeze | sort > verification/requirements-lock.txt

echo "Verification run complete."
