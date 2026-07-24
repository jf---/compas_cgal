#!/usr/bin/env bash
# Compile + run the arrangement-redesign exploration prototypes.
# CGAL is header-only but links GMP+MPFR; use the cgal-dev env (CGAL 5.6.1).
# Override CGAL_ENV if your env lives elsewhere.
set -euo pipefail
CGAL_ENV="${CGAL_ENV:-/Users/jelle/mambaforge/envs/cgal-dev}"
cd "$(dirname "$0")"
for src in min pl_bench stock_arr gps_zone harvest; do
    echo "=== $src ==="
    c++ -std=c++20 -O2 -I"$CGAL_ENV/include" "$src.cpp" -o "$src" \
        -L"$CGAL_ENV/lib" -lgmp -lmpfr
    DYLD_LIBRARY_PATH="$CGAL_ENV/lib" "./$src"
    echo
done
