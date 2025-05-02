#!/bin/bash
echo '========================================================================='
echo '======================== on the GPU:'
echo 'First analysis/factorize takes longer, to "warmup" the GPU, which is the'
echo 'time to allocate and pin GPU/CPU memory pools.  That only needs to be'
echo 'done once for the entire application, however.'
echo '========================================================================='
set -e
cd ..
make dl_simple
cd Demo
CHOLMOD_USE_GPU=1 CHOLMOD_MATRIXOPS=ON \../build/cholmod_dl_simple \
 -A ../../../RXMesh-AMGX/output/sphere3/A_sym.mtx \
 -B ../../../RXMesh-AMGX/output/sphere3/B.mtx \
 -EX ../../../RXMesh-AMGX/output/sphere3/X.mtx

