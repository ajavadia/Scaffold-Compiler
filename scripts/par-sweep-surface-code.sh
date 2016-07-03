#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)

# surface codes
attempt_th_yx=(4 8)
attempt_th_drop=(15 20)
layout=("" "--drop")

/usr/bin/time parallel -j32 $ROOT/braidflash/braidflash {1} --p {2} --yx {3} --drop {4} {5} ::: $(ls $1/*lpfs | sed "s/^\(.*\)\..*$/\1/") ::: ${error_rates[@]} ::: ${attempt_th_yx[@]} ::: ${attempt_th_drop[@]} ::: "${layout[@]}"
