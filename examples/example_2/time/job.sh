#!/bin/bash

for number in 4 8 16 32 64 128 256 512
do
    g++ -Wall -Wextra -pedantic -O3 -std=c++20 -pthread -DNDEBUG -I../../../ \
        example_2.cc -DNUMBER=${number} -o example_2
    for i in `seq 1 10`
    do
	echo -n "${number} "
	(time ./example_2 ) |& grep real | awk '{print $2}'
    done
done | tee results.dat
rm example_2 solution.dat
