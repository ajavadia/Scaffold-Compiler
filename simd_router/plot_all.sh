#!/bin/bash

# usage: $ bash plot_all.sh {benchmark1} {benchmark2} 

#compile
rm -f router && make

# change the desired thresholds below. that's the only change needed.
error_rates=(6 9)
caps=(1 5 15)
windows=(20)
caps_for_windows=(2)
windows_for_caps=(100)
directions=("forward" "backforth")

#----------- Different Smoothing Directions
for bench in $*; do
  bench_dir=$(dirname $bench)
  bench_name=$(basename $bench)
  output_simulation=$bench_dir/output_simulation/$bench_name
  output_plot=$bench_dir/output_plot/$bench_name
  rm -rf $bench_dir/output_simulation && mkdir -p $bench_dir/output_simulation
  rm -rf $bench_dir/output_plot && mkdir -p $bench_dir/output_plot
  
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
        echo "Cycle WindowSize=$window" | cat - ${output_simulation}.cap.inf.window.${window}.${direction}.usage > /tmp/out && mv /tmp/out ${output_plot}.cap.inf.window.${window}.${direction}.usage #usage
        echo "$window" | paste - ${output_simulation}.cap.inf.window.${window}.${direction}.ages > tmp && mv tmp ${output_plot}.cap.inf.window.${window}.${direction}.ages                    #ages
        cat ${output_plot}.${direction}.ages.data ${output_plot}.cap.inf.window.${window}.${direction}.ages > temp && mv temp ${output_plot}.${direction}.ages.data                     #ages
        if [ $window -eq ${windows[0]} ];
        then
          cat ${output_plot}.cap.inf.window.${window}.${direction}.usage | paste ${output_plot}.${direction}.usage.data - > temp && mv temp ${output_plot}.${direction}.usage.data
        else
          cat ${output_plot}.cap.inf.window.${window}.${direction}.usage | awk '{print $2}' - | paste ${output_plot}.${direction}.usage.data - > temp && mv temp ${output_plot}.${direction}.usage.data
        fi
      done

      echo "WindowSize = inf"
      ./router ${bench} --p ${p} --${direction} --usage --ages --storage
      echo "cycle WindowSize=inf" | cat - ${output_simulation}.cap.inf.window.inf.${direction}.usage > /tmp/out && mv /tmp/out ${output_plot}.cap.inf.window.inf.${direction}.usage
      cat ${output_plot}.cap.inf.window.inf.${direction}.usage | awk '{print $2}' - | paste ${output_plot}.${direction}.usage.data - > temp && mv temp ${output_plot}.${direction}.usage.data  #usage
      echo "inf" | paste - ${output_simulation}.cap.inf.window.inf.${direction}.ages > tmp && mv tmp ${output_plot}.cap.inf.window.inf.${direction}.ages             #ages
      cat ${output_plot}.${direction}.ages.data ${output_plot}.cap.inf.window.inf.${direction}.ages > temp && mv temp ${output_plot}.${direction}.ages.data                              #ages

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

      set terminal postscript landscape dashed color

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

      set terminal postscript landscape dashed color

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
      
      rm plot_*.gp ${output_plot}.*.window* 

      #open ${output_plot}.${direction}.usage.ps
      #open ${output_plot}.${direction}.ages.ps

  done
done