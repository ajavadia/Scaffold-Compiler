#!/bin/bash

physical_errors=(3 4 5 6 7 8 9)

#compile
rm -f braidflash && make

#run
for b in $(cat file_list)
do
  echo $b  
  for p in "${physical_errors[@]}"
  do
    echo $p
    ./braidflash ${b} --p ${p} > ${b}.p${p}.br  
    ./braidflash ${b} --opt --p ${p} > ${b}.p${p}.opt.br 
    ./braidflash ${b} --cnot --p ${p} > ${b}.p${p}.cnot.br  
    ./braidflash ${b} --cnot --opt --p ${p} > ${b}.p${p}.cnot.opt.br    
  done
done
