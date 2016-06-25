#!/bin/bash

DIR=$(dirname $0)
ROOT=$DIR/..
OPT=$ROOT/build/Release+Asserts/bin/opt
SCAF=$ROOT/build/Release+Asserts/lib/Scaffold.so

# Capacity of each SIMD region
D=(1024)
# Number of SIMD regions
K=(4)
# Module flattening threshold
# note: thresholds must be picked from the set in scripts/flattening_thresh.py
THRESHOLDS=(001k)
# Full schedule? otherwise only generates metrics (faster)
FULL_SCHED=true

# Create directory to put all byproduct and output files in
for f in $*; do
  b=$(basename $f .scaffold)  
  echo "[gen-lpfs.sh] $b: Creating output directory ..."
  mkdir -p "$b"
  #mv ./*${b}* ${b} 2>/dev/null
done

# Generate .ll file if not done already
for f in $*; do
  b=$(basename $f .scaffold)
  echo "[gen-lpfs.sh] $b: Compiling ..."
  if [ ! -e ${b}/${b}.ll ]; then
    # Generate compiled files
    $ROOT/scaffold.sh -r $f
    mv ${b}11.ll ${b}11.ll.keep_me
    # clean intermediary compilation files (comment out for speed)
    $ROOT/scaffold.sh -c $f
    # Keep the final output for the compilation
    mv ${b}11.ll.keep_me ${b}/${b}.ll
  fi
done

# Module flattening pass with different thresholds
for f in $*; do
  b=$(basename $f .scaffold)
  echo "[gen-lpfs.sh] $b: Flattening ..."
  echo "[gen-lpfs.sh] Computing module gate counts ..."  
  $OPT -S -load $SCAF -ResourceCount2 ${b}/${b}.ll > /dev/null 2> ${b}.out  
  python $DIR/flattening_thresh.py ${b}  
  for th in ${THRESHOLDS[@]}; do      
    if [ ! -e ${b}/${b}.flat${th}.ll ]; then
      echo "[gen-lpfs.sh] Flattening modules smaller than Threshold = $th ..."    
      mv ${b}.flat${th}.txt flat_info.txt
      $OPT -S -load $SCAF -FlattenModule -dce -internalize -globaldce ${b}/${b}.ll -o ${b}/${b}.flat${th}.ll
    fi
  done
  rm -f *flat*.txt ${b}.out     
done

# Perform resource estimation 
for f in $*; do
  b=$(basename $f .scaffold)
  echo "[gen-lpfs.sh] $b: Resource count ..."
  for th in ${THRESHOLDS[@]}; do      
    if [ -n ${b}/${b}.flat${th}.resources ]; then
      echo "[gen-lpfs.sh] Resource count for Threshold = $th flattening ..."
      $OPT -S -load $SCAF -ResourceCount ${b}/${b}.flat${th}.ll > /dev/null 2> ${b}/${b}.flat${th}.resources
    fi
  done
done

# For different K and D values specified above, generate MultiSIMD schedules
for f in $*; do
  b=$(basename $f .scaffold)
  for d in ${D[@]}; do
    for k in ${K[@]}; do
      echo "[gen-lpfs.sh] $b: Generating SIMD K=$k D=$d leaves ..."
      for th in ${THRESHOLDS[@]}; do
        if [ ! -e ${b}/${b}.flat${th}.simd.${k}.${d}.leaves.local ]; then
          echo "[gen-lpfs.sh] GenSIMD for Threshold = $th flattening ..."
          $OPT -load $SCAF -GenLPFSSchedule -simd-kconstraint-lpfs $k -simd-dconstraint-lpfs $d -simd_l 1 -local_mem 1 ${b}/${b}.flat${th}.ll > /dev/null 2> ${b}/${b}.flat${th}.simd.${k}.${d}.leaves.local
        fi
      done
    done
  done
done

# Take into account the penalty of ballistic communication
for f in $*; do
  b=$(basename $f .scaffold)
  cd ${b}
  echo "[gen-lpfs.sh] $b: Adding communication latencies ..."
  ../${DIR}/comm_aware.pl ${b}*.leaves.local
  cd ..
done

# Obtain coarse-grain schedules by co-scheduling modules
for f in $*; do
  b=$(basename $f .scaffold)
  cd ${b}
  for th in ${THRESHOLDS[@]}; do      
    for c in comm_aware_schedule.txt.${b}.flat${th}_*; do
      k=$(perl -e '$ARGV[0] =~ /_K(\d+)/; print $1' $c)
      d=$(perl -e '$ARGV[0] =~ /_D(\d+)/; print $1' $c)
      x=$(perl -e '$ARGV[0] =~ /.*_(.+)/; print $1' $c)
      th=$(perl -e '$ARGV[0] =~ /.flat(\d+[a-zA-Z])/; print $1' $c)    
      echo "[gen-lpfs.sh] $b: Coarse-grain schedule ..."
      mv $c comm_aware_schedule.txt
      if [ ! -e ${b}.flat${th}.simd.${k}.${d}.${x}.time ]; then
        ../$OPT -load ../$SCAF -GenCGSIMDSchedule -simd-kconstraint-cg $k -simd-dconstraint-cg $d ${b}.flat${th}.ll > /dev/null 2> ${b}.flat${th}.simd.${k}.${d}.${x}.time
      fi
    done
  done
  rm -f comm_aware_schedule.txt
  cd ..
done
