#!/bin/bash

file=evolution.dat
temp=temp.dat
out=plot.png

for i in `cat evolution.dat | awk '{ print $1 }' | uniq`
do
  cat ${file} | grep "^${i} " > ${temp}
  echo "\
    reset; \
    set samples 1000, 1000; \
    set xrange [-10: 10]; \
    set yrange [-10: 10]; \
    r(x, y) = sqrt(x * x + y * y); \
    f(x, y) = cos(0.25 * r(x, y)) + exp(1.0); \
    set t png; \
    set o \"${out}\"; \
    set title \"${i}\"; \
    set pm3d; \
    set cbtics 0.5; \
    unset border; \
    unset xtics; \
    unset ytics; \
    unset ztics; \
    set view 60, 350; \
    splot f(x, y) notitle, \
          \"${temp}\" u 2:3:(f(\$2, \$3)) pt 7 notitle; \
    set o" | gnuplot
  mv ${out} `printf '%03d.png' ${i}`
done

exit
