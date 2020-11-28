#!/bin/bash

file=evolution.dat
temp=chessboard.tex

echo "\documentclass{minimal}" > ${temp}
echo "\usepackage[utf8]{inputenc}" >> ${temp}
echo "\usepackage{skak}" >> ${temp}
echo "\begin{document}" >> ${temp}
echo "\newgame" >> ${temp}
echo -n "\fenboard{" >> ${temp}
echo -n `cat ${file} \
    | awk -F' ' '$NF=="0"' \
    | awk '{print $(NF-6), $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1)}' \
    | head -n 1` >> ${temp}
echo "}" >> ${temp}
echo "\showboard" >> ${temp}
echo "\end{document}" >> ${temp}

pdflatex ${temp}
pdfcrop ${temp%.*}.pdf result.pdf
rm ${temp%.*}*

exit
