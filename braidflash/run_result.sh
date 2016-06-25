#!/bin/sh

#run
for b in $(cat file_list)
do
  f=$b
  echo $f

  # Excel Sheet 'BraidFlash'
  echo "TotalSUCCESS" "TotalCONFLICT" "UniqueCONFLICT" "SerialCLOCK" "ParallelCLOCK" "CriticalCLOCK"  
  echo "-- CNOT only"
  Ser1=$(grep "SerialCLOCK" $f.cnot.br | awk '{sum+=$2} END {print sum}')  
  Par1=$(grep "ParallelCLOCK" $f.p3.cnot.br | awk '{sum+=$2} END {print sum}')  
  Cri1=$(grep "CriticalCLOCK" $f.cnot.br | awk '{sum+=$2} END {print sum}')    
  tots1=$(grep "total_success" $f.cnot.br | awk '{sum+=$2} END {print sum}')  
  totc1=$(grep "total_conflict" $f.cnot.br | awk '{sum+=$2} END {print sum}') 
  unic1=$(grep "unique_conflict" $f.cnot.br | awk '{sum+=$2} END {print sum}')   

  Ser2=$(grep "SerialCLOCK" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')  
  Par2=$(grep "ParallelCLOCK" $f.p3.cnot.opt.br | awk '{sum+=$2} END {print sum}')  
  Cri2=$(grep "CriticalCLOCK" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')    
  tots2=$(grep "total_success" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')  
  totc2=$(grep "total_conflict" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}') 
  unic2=$(grep "unique_conflict" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')  
  echo $tots1 $totc1 $unic1 $Ser1 $Par1 $Cri1
  echo $tots2 $totc2 $unic2 $Ser2 $Par2 $Cri2

  echo "-- CNOT+H"
  Ser3=$(grep "SerialCLOCK" $f.br | awk '{sum+=$2} END {print sum}')  
  Par3=$(grep "ParallelCLOCK" $f.p3.br | awk '{sum+=$2} END {print sum}')  
  Cri3=$(grep "CriticalCLOCK" $f.br | awk '{sum+=$2} END {print sum}')    
  tots3=$(grep "total_success" $f.br | awk '{sum+=$2} END {print sum}')  
  totc3=$(grep "total_conflict" $f.br | awk '{sum+=$2} END {print sum}') 
  unic3=$(grep "unique_conflict" $f.br | awk '{sum+=$2} END {print sum}')   

  Ser4=$(grep "SerialCLOCK" $f.opt.br | awk '{sum+=$2} END {print sum}')  
  Par4=$(grep "ParallelCLOCK" $f.p3.opt.br | awk '{sum+=$2} END {print sum}')  
  Cri4=$(grep "CriticalCLOCK" $f.opt.br | awk '{sum+=$2} END {print sum}')    
  tots4=$(grep "total_success" $f.opt.br | awk '{sum+=$2} END {print sum}')  
  totc4=$(grep "total_conflict" $f.opt.br | awk '{sum+=$2} END {print sum}') 
  unic4=$(grep "unique_conflict" $f.opt.br | awk '{sum+=$2} END {print sum}')  
  echo $tots3 $totc3 $unic3 $Ser3 $Par3 $Cri3
  echo $tots4 $totc4 $unic4 $Ser4 $Par4 $Cri4

  # Excel Sheet 'DroppedGates'
  echo "TotalDroppedGates" "UniqueDroppedGates" 
  echo "-- CNOT only"
  totd1=$(grep "total_dropped_gates" $f.cnot.br | awk '{sum+=$2} END {print sum}')  
  unid1=$(grep "unique_dropped_gates" $f.cnot.br | awk '{sum+=$2} END {print sum}')    
  totd2=$(grep "total_dropped_gates" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')  
  unid2=$(grep "unique_dropped_gates" $f.cnot.opt.br | awk '{sum+=$2} END {print sum}')      
  echo $totd1 $unid1
  echo $totd2 $unid2
  
  echo "-- CNOT+H"
  totd3=$(grep "total_dropped_gates" $f.br | awk '{sum+=$2} END {print sum}')  
  unid3=$(grep "unique_dropped_gates" $f.br | awk '{sum+=$2} END {print sum}')    
  totd4=$(grep "total_dropped_gates" $f.opt.br | awk '{sum+=$2} END {print sum}')  
  unid4=$(grep "unique_dropped_gates" $f.opt.br | awk '{sum+=$2} END {print sum}')      
  echo $totd3 $unid3
  echo $totd4 $unid4

  # Excel Sheet 'ConflictedAttempts'
  echo "Attempts"  
  echo "-- CNOT+H"
  for a in {0..40}
  do
    att_sum=$(grep "attempt\t$a" $f.br | awk '{sum+=$3} END {print sum}') 
    att_sum_opt=$(grep "attempt\t$a" $f.opt.br | awk '{sum+=$3} END {print sum}') 
    echo $a $att_sum $att_sum_opt
  done  

  # Excel Sheet 'ManhattanCost'
  echo "ManhattanCost" "ManhattanCost.Opt"
  mcost1=$(grep "mcost:" $f.opt.br | awk '{print $2}')    
  mcost2=$(grep "mcost_opt:" $f.opt.br | awk '{print $2}')      
  echo $mcost1 $mcost2

  # Excel Sheet 'Area'
  echo "CodeDistance" "MaxLogicalQbits" "MaxPhysicalQbits"
  d=$(grep "code_distance(d):" $f.br | awk '{print $2}')
  Q=$(grep "num_logical_qubits:" $f.br | awk '{print $2}')
  q=$(grep "num_physical_qubits:" $f.br | awk '{print $2}')
  echo $d $Q $q
  # Extra Info 
  echo "EventCount.CNOT" "EventCount.CNOT+H"
  ecount1=$(grep "event_count:" $f.cnot.opt.br | awk '{print $2}')     
  ecount2=$(grep "event_count:" $f.opt.br | awk '{print $2}')        
  echo $ecount1 $ecount2
  
done
