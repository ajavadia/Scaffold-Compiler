#!/bin/bash

# usage: $ bash plot_all.sh {benchmark1} {benchmark2} 

#compile
#rm -f router && make

# change the desired thresholds below. that's the only change needed.
error_rates=(5)
caps=("inf")
windows=(10 100 inf)
directions=("forward" "backforth")

#----------- Different Smoothing Directions
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  output_simulation=$bench_dir/simd_simulation/$bench_name
  output_plot=$bench_dir/simd_simulation/simd_plot/$bench_name
  rm -rf $bench_dir/simd_simulation/simd_plot && mkdir -p $bench_dir/simd_simulation/simd_plot

  for p in "${error_rates[@]}"
  do
    for direction in "${directions[@]}"
    do  
      echo ${direction}
      echo -n "" > ${output_plot}.${direction}.usage.data
      echo -n "" > ${output_plot}.${direction}.ages.data
      echo "WindowSize avg peak" > ${output_plot}.${direction}.ages.data  
      for window in "${windows[@]}"
      do
        echo "WindowSize = $window"
        ./router ${bench} --p ${p} --window ${window} --${direction} --usage --ages --storage
        echo "Cycle WindowSize=$window" | cat - ${output_simulation}.p.$p.cap.inf.window.${window}.${direction}.usage > /tmp/out && mv /tmp/out ${output_plot}.p.$p.cap.inf.window.${window}.${direction}.usage #usage
        echo "$window" | paste - ${output_simulation}.p.$p.cap.inf.window.${window}.${direction}.ages > tmp && mv tmp ${output_plot}.p.$p.cap.inf.window.${window}.${direction}.ages                    #ages
        cat ${output_plot}.${direction}.ages.data ${output_plot}.p.$p.cap.inf.window.${window}.${direction}.ages > temp && mv temp ${output_plot}.${direction}.ages.data                     #ages
        if [ "$window" == "${windows[0]}" ];
        then
          cat ${output_plot}.p.$p.cap.inf.window.${window}.${direction}.usage | paste ${output_plot}.${direction}.usage.data - > temp && mv temp ${output_plot}.${direction}.usage.data
        else
          cat ${output_plot}.p.$p.cap.inf.window.${window}.${direction}.usage | awk '{print $2}' - | paste ${output_plot}.${direction}.usage.data - > temp && mv temp ${output_plot}.${direction}.usage.data
        fi
      done

      #----------- usage: create gnuplot script and plot
      # a string that will basically say plot the first column agains columns 2,3,...
      plot_string="plot "
      for i in `seq 1 ${#windows[@]}`
      do
        plot_string+="\"${output_plot}.${direction}.usage.data\" u 1:$((i+1)) w lines lt $((i+1)) title col, "
      done
      plot_string+="\"${output_plot}.${direction}.usage.data\" u 1:$((i+2)) w lines lt $((i+2)) title col"

      # put it all into a .gp script
      echo "
      set output \"${output_plot}.${direction}.usage.ps\"

      set terminal postscript landscape color

      set size 1,.75
      set key outside right top

      set title 'Qubit Usage over Time (${direction} smoothing)'
      set xlabel 'cycles'
      set ylabel 'number of existing qubits'

      set style line 1 lc rgb \"blue\"
      set boxwidth 0.5
      " > plot_usage.gp
      echo ${plot_string} >> plot_usage.gp
      # run the .gp script with gnuplot, saving the plot
      gnuplot plot_usage.gp


      #----------- ages: create gnuplot script and plot
      # a string that will just plot the avg and the peak
      plot_string="plot \"${output_plot}.${direction}.ages.data\" u 2:xtic(1) w boxes lt 2 title col, \"${output_plot}.${direction}.ages.data\" u 3:xtic(1) w boxes lt 3 title col"

      # put it all into a .gp script
      echo "
      set output \"${output_plot}.${direction}.ages.ps\"

      set terminal postscript landscape color

      set size 1,.75
      set key outside right top

      set title 'Lifespan of Qubits (${direction} smoothing)'
      set xlabel 'smoothing window size'
      set ylabel 'qubit age at time of death (cycles)'

      set style data histogram
      set auto x
      set boxwidth 0.5
      " > plot_ages.gp
      echo ${plot_string} >> plot_ages.gp
      # run the .gp script with gnuplot, saving the plot
      gnuplot plot_ages.gp

      rm -f plot_*.gp ${output_plot}.*.window* 

      echo "usage plot written to: " ${output_plot}.${direction}.usage.ps
      echo "ages plot written to: " ${output_plot}.${direction}.ages.ps   
      ps2pdf ${output_plot}.${direction}.usage.ps ${output_plot}.${direction}.usage.pdf
      ps2pdf ${output_plot}.${direction}.ages.ps ${output_plot}.${direction}.ages.pdf
      #open ${output_plot}.${direction}.usage.ps
      #open ${output_plot}.${direction}.ages.ps

    done
  done
done
