#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)
# concatenated codes
caps=(5 20 100)
windows=(10 50)
directions=("forward" "backward" "backforth")
# surface codes
attempt_th_yx=(8 20)
attempt_th_drop=(4 15)

# sweep problem sizes
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  simd_simulation=$bench_dir/simd_simulation
  braid_simulation=$bench_dir/braid_simulation
  
  # sweep error rates
  for p in "${error_rates[@]}"
  do
    echo "================== p = ${p} =================="    

    echo "------------- Concatenated Code --------------"       
    list=(${simd_simulation}/${bench_name}.p.${p}.*)
    # find best concatenated kq
    best_kq=$(grep "physical kq:" "$list" | awk '{print $3}')
    best_f=""
    for f in ${simd_simulation}/${bench_name}.p.${p}.*
    do    
      kq=$(grep "physical kq:" $f | awk '{print $3}')
      KQ=$(grep "logical KQ:" $f | awk '{print $3}')
      if ((kq < best_kq)); then
        best_kq=$kq
        best_f=$f
      fi
    done
    echo "best_cc_f: " $best_f
    echo "best_cc_kq: " $best_kq
    echo "KQ: " $KQ

    echo "--------------- Surface Code -----------------"   
    list=(${braid_simulation}/${bench_name}.p.${p}.*)    
    # find best surface kq    
    best_kq=$(grep "physical kq:" "$list" | awk '{print $3}')
    best_f=""
    for f in ${braid_simulation}/${bench_name}.p.${p}.*
    do    
      kq=$(grep "physical kq:" $f | awk '{print $3}')
      KQ=$(grep "logical KQ:" $f | awk '{print $3}')
      if ((kq < best_kq)); then
        best_kq=$kq
        best_f=$f
      fi
    done
    echo "best_sc_f: " $best_f
    echo "best_sc_kq: " $best_kq
    echo "KQ: " $KQ
  done
done
