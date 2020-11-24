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
    f(x) = sin(2 * x) * exp(-0.05 * x * x) + pi; \
    set t png; \
    set o \"${out}\"; \
    set title \"${i}\"; \
    plot \"${temp}\" u 2:(f(\$2)) notitle, f(x) notitle; \
    set o" | gnuplot
  mv ${out} `printf '%03d.png' ${i}`
done

gst-launch-1.0 -e multifilesrc location="%03d.png" \
  index=0 caps="image/png,framerate=(fraction)10/1" \
  ! pngdec ! videoconvert ! theoraenc ! oggmux \
  ! filesink location=result.ogv

rm ${temp} *.png

exit

