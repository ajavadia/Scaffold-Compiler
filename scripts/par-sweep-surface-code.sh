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

# surface codes parameters
attempt_th_yx=(4 8)
attempt_th_drop=(15 20)
layout=("" "--opt")
technologies=("ion" "sup")

parallel="parallel --delay .2 -j 8 --joblog logs/sc.log --resume"

$parallel "$ROOT/braidflash/braidflash {1} --p {2} --yx {3} --drop {4} {5} --tech {6} >/dev/null" ::: ${files[@]} ::: ${error_rates[@]} ::: ${attempt_th_yx[@]} ::: ${attempt_th_drop[@]} ::: "${layout[@]}" ::: ${technologies[@]}
