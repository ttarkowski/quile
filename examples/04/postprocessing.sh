#!/bin/bash

file=solution.dat
temp=chessboard.tex

echo "\documentclass{minimal}" > ${temp}
echo "\usepackage[utf8]{inputenc}" >> ${temp}
echo "\usepackage{skak}" >> ${temp}
echo "\begin{document}" >> ${temp}
echo "\newgame" >> ${temp}
echo -n "\fenboard{" >> ${temp}
echo -n `cat ${file}` >> ${temp}
echo "}" >> ${temp}
echo "\showboard" >> ${temp}
echo "\end{document}" >> ${temp}

pdflatex ${temp}
pdfcrop ${temp%.*}.pdf result.pdf
rm ${temp%.*}*

exit
