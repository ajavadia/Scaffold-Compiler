#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

#usage
#./par-sweep-concatenated.sh <path/to/DIR>

# Benchmark to sweep
bench="../LPFS_Scheds/010k_square_root"

files=($(ls $bench/*lpfs | sed "s/^\(.*\)\..*$/\1/"))

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)

# concatenated codes parameters
caps=(10 100 "inf")
windows=(10 100 "inf")
directions=("--forward" "--backward")
technologies=("ion" "superconductor")

parallel="parallel --delay .2 -j 8 --joblog logs/cc.log --resume"

$parallel "$ROOT/simd_router/router {1} --p {2} --cap {3} --window {4} {5} --tech {6} --usage --ages >/dev/null 2>/dev/null" ::: ${files[@]} ::: ${error_rates[@]} ::: ${caps[@]} ::: ${windows[@]} ::: ${directions[@]} ::: ${technologies[@]}
