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
directions=("backforth")

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
      echo "Direction = ${direction}"
      for window in "${windows[@]}"
      do
        echo "WindowSize = ${window}"        
        for cap in "${caps[@]}"
        do
          echo "CapSize = ${cap}"
          $ROOT/simd_router/router ${bench} --p ${p} --window ${window} --cap ${cap} --${direction}
          echo "CapSize = inf"                            
          $ROOT/simd_router/router ${bench} --p ${p} --window ${window} --${direction}          
        done
      done
      echo "WindowSize = inf"        
      for cap in "${caps[@]}"
      do
        echo "CapSize = ${cap}"                  
        $ROOT/simd_router/router ${bench} --p ${p} --cap ${cap} --${direction}
        echo "CapSize = inf"                            
        $ROOT/simd_router/router ${bench} --p ${p} --${direction}          
      done      
    done
  done
done
