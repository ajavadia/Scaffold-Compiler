#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)

# concatenated codes parameters
caps=(100)
windows=(50)
directions=("forward" "backward" "backforth")

# sweep problem sizes
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  simd_simulation=$bench_dir/simd_simulation
  
  # sweep error rates
  for p in "${error_rates[@]}"
  do
    echo "================== p = ${p} =================="
    echo "------------- Concatenated Code --------------"        
    # simulate simd_router with concatenated code (different directions, caps, windows)
    for direction in "${directions[@]}"
    do  
      for window in "${windows[@]}"
      do     
        for cap in "${caps[@]}"
        do
          echo "(p = ${p}, CapSize = ${cap}, WindowSize = ${window}, Direction = ${direction})"
          $ROOT/simd_router/router ${bench} --p ${p} --window ${window} --cap ${cap} --${direction}
        done
      done
    done
  done
done
