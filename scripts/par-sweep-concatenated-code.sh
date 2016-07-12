#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

#usage
#./par-sweep-concatenated.sh <path/to/DIR>

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)

# concatenated codes parameters
caps=(80 100)
windows=(40 50)
directions=("forward" "backward" "backforth")

/usr/bin/time parallel -j32 $ROOT/simd_router/router {1} --p {2} --cap {3} --window {4} --direction {5} ::: $(ls $1/*lpfs | sed "s/^\(.*\)\..*$/\1/") ::: ${error_rates[@]} ::: ${caps[@]} ::: ${windows[@]} ::: ${directions[@]}
