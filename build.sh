#!/usr/bin/env bash
set -e
mkdir -p build
g++ -O3 -march=native -DNDEBUG src/fiber_fem_nedelec.cpp -o build/fiber_fem -I/usr/include/eigen3
echo "OK -> build/fiber_fem"
