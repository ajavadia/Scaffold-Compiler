#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# change the desired thresholds below. that's the only change needed.
error_rates=(5 6 7 8 9)
technologies=("ion" "sup")

# sweep problem sizes
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  simd_simulation=$bench_dir/simd_simulation
  braid_simulation=$bench_dir/braid_simulation
 
  for tech in "${technologies[@]}" 
  do
    # sweep error rates
    for p in "${error_rates[@]}"
    do
      echo "================== p = ${p} =================="    

      echo "------------- Concatenated Code --------------"       
      list=(${simd_simulation}/${bench_name}.p.${p}.*.${tech}*.kq)
      # find best concatenated kq
      best_cc_kq=$(grep "physical kq:" "$list" | awk '{print $3}')
      best_cc_f="$list"
      for f in ${simd_simulation}/${bench_name}.p.${p}.*.${tech}*.kq
      do   
        kq=$(grep "physical kq:" $f | awk '{print $3}')
        KQ=$(grep "logical KQ:" $f | awk '{print $3}')
        concat_level=$(grep "concatenation level:" $f | awk '{print $3}')
        if ((kq < best_cc_kq)); then        
          best_cc_kq=$kq
          best_cc_f=$f
        fi
      done
      echo "best_cc_f: " $best_cc_f
      echo "best_cc_kq: " $best_cc_kq
      echo "KQ: " $KQ

      echo "--------------- Surface Code -----------------"   
      list=(${braid_simulation}/${bench_name}.p.${p}.*.${tech}*.kq)    
      # find best surface kq    
      best_sc_kq=$(grep "physical kq:" "$list" | awk '{print $3}')
      best_sc_f="$list"
      for f in ${braid_simulation}/${bench_name}.p.${p}.*.${tech}*.kq
      do 
        kq=$(grep "physical kq:" $f | awk '{print $3}')
        KQ=$(grep "logical KQ:" $f | awk '{print $3}')
        code_distance=$(grep "code distance:" $f | awk '{print $3}')      
        if ((kq < best_sc_kq)); then
          best_sc_kq=$kq
          best_sc_f=$f
        fi
      done
      echo "best_sc_f: " $best_sc_f
      echo "best_sc_kq: " $best_sc_kq
      echo "KQ: " $KQ

      echo "---------------- Comparison ------------------"   
      cc_to_sc=$(bc <<< "$best_cc_kq/$best_sc_kq")
      if ((cc_to_sc < 1)); then
        best_for_p="CONCATENATED"
      else
        best_for_p="SURFACE"
      fi
      echo "(concatenation_level code_distance):" "($concat_level $code_distance)"
      echo $best_for_p "($(bc <<< "scale=4;$best_cc_kq/$best_sc_kq"))"
    done
  done
done

