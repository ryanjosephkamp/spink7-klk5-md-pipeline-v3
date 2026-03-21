#!/usr/bin/env bash
# Run the full test suite under Python optimization mode (-O).
# This strips all assert statements from bytecode, exposing any
# validation logic that incorrectly relies on assert.
#
# Prerequisite: L-20 must be implemented first (replace assert-based
# validation with proper exception classes). If L-20 is not yet
# implemented, some tests WILL fail — this is expected and diagnostic.
#
# Usage:
#   bash scripts/test_optimized.sh
#   bash scripts/test_optimized.sh tests/test_integration.py  # specific file
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

echo "========================================="
echo "Running pytest under python -O"
echo "Python: $(python --version)"
echo "Optimization level: -O (assert stripped)"
echo "========================================="

python -O -m pytest "${@:-tests/}" -v --tb=short

echo "========================================="
echo "All tests passed under -O"
echo "========================================="
