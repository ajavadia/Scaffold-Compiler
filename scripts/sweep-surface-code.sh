#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# change the desired thresholds below. that's the only change needed.
error_rates=(9)

# surface codes
attempt_th_yx=(4)
attempt_th_drop=(15)

# sweep problem sizes
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  braid_simulation=$bench_dir/braid_simulation
  
  # sweep error rates
  for p in "${error_rates[@]}"
  do
    echo "================== p = ${p} =================="
    # simulate braidflash with surface code (different opt, attempt_th_yx, attempt_th_drop)
    echo "--------------- Surface Code -----------------"   
    for yx in "${attempt_th_yx[@]}"
    do  
      for drop in "${attempt_th_drop[@]}"
      do
        echo "(p = ${p}, attempt_th_yx = ${yx}, attempt_th_drop = ${drop})"        
        cd $ROOT/braidflash
        ./braidflash ${bench} --p ${p} --yx ${yx} --drop ${drop} > /dev/null
        ./braidflash ${bench} --p ${p} --yx ${yx} --drop ${drop} --opt > /dev/null
        cd -
      done
    done
  done
done
