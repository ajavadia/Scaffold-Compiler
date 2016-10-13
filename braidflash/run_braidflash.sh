#!/bin/bash

# usage: $ bash run_braidflash.sh {benchmark1} {benchmark2} 

physical_errors=(5)
technologies=("sup")
priority_policies=(0 1 2 3 4 5 6)

#compile
make

#run
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)

  for tech in "${technologies[@]}"   
  do
    for p in "${physical_errors[@]}"
    do 
      for pri in "${priority_policies[@]}"
      do
      ./braidflash ${bench} --p ${p} --tech ${tech} --pri ${pri}
      ./braidflash ${bench} --p ${p} --tech ${tech} --pri ${pri} --opt
      done
    done
  done
done
