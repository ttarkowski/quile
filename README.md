# Quilë
Quilë---C++20 header-only genetic algorithms library released under the terms
of the MIT License.

## Etymology
> qya. *quilë* -- color

The part *chromo-* of the word *chromosome*, which is central term in genetics,
comes from the Greek χρῶμα and means *color*. Similarly, the name of this
library origins from fictional language Neo-Quenya and means *color* as well.

## Features
- Gene- and genotype-level constraints
- Representation: floating-point, integer, binary, permutation
- Mutation operators: Gaussian, swap, random reset
- Recombination operators: arithmetic, one-point crossover, cut-and-crossfill
- Selection probabilities calculations method: fitness proportional selection
- First generation: random, user input
- Population selection schemes: roulette wheel selection, stochastic universal
  sampling (SUS), generational
- Fitness function calculations: database of calculated values, support for
  fitness functions with incalculable results
- Termination conditions: fixed number of iterations, based on maximum fitness
  improvement, fitness function treshold
- Multiple thread calculations
- Set of test functions for benchmark analysis of floating-point optimization

## Dependencies
- Modern compiler with at least partial support of C++20 (library is known to
  work with GCC 10.2.0)

*Please note that code in examples/ directory might require some additional
software.*

## Tutorial and documentation
Please see [Wiki](https://github.com/ttarkowski/quile/wiki).

## Examples
The library is platform independent, but some postprocessing of example's data
requires software specific for GNU/Linux. Hence, in order to follow the
instructions below, you will need GNU/Linux operating system.

1. Directory change (please replace *nn* with example number)  
    *cd examples/nn*
2. Code compilation  
    *g++ -Wall -Wextra -pedantic -O3 -std=c++20 -fconcepts -pthread -I../../
    example.cc*
3. Program execution  
    *./a.out*
4. Postprocessing (please see note below)  
    *chmod u+x ./postprocessing.sh*  
    *./postprocessing.sh*
5. Result view -- open the file called *result.ogv* or *result.pdf*

Note: Depending on the example you will need some additional free software
packages listed below.

List of examples:

- 01: **maximum of 1D function**  
Representation: *floating-point*  
Requires: *gnuplot, GStreamer*  
Result: *OGV movie showing genetic drift*
- 02: **maximum of 2D function**  
Representation: *floating-point*  
Requires: *gnuplot, GStreamer*  
Result: *OGV movie showing genetic drift*
- 03: **the eight queens puzzle**  
Representation: *integer*  
Requires: *TeXLive with skak package and pdfcrop program*  
Result: *PDF file showing chessboard with solution*
- 04: **the eight queens puzzle**  
Representation: *permutation*  
Requires: *TeXLive with skak package and pdfcrop program*  
Result: *PDF file showing chessboard with solution*

## Benchmark analysis
Example benchmark analysis program is located at examples/benchmark/ directory.
It can be compiled and executed with following commands:  
    *g++ -Wall -Wextra -pedantic -O3 -std=c++20 -fconcepts -pthread -I../../
    example.cc*  
    *./a.out*  
Please note that contrary to other examples there is no postprocessing procedure
nor special requirements.

## Evenstar
Evenstar is an example program from computational materials science domain and
is devoted to crystal structure prediction of boron nanowires. It is located at
examples/evenstar/ directory and requires additional software and data:  

- Boost Graph Library (BGL) - compilation dependency
- Quantum ESPRESSO - runtime dependency
- B.pbe-n-kjpaw_psl.1.0.0.UPF pseudopotential file - runtime dependency

The program can be compiled with example following comand:  
    *g++ -Wall -Wextra -pedantic -O3 -std=c++20 -fconcepts -pthread -I../../
    main.cc src/\*.cc -DFLAT=true -DCELL_ATOMS=10*  
Please note that -DFLAT parameter, describing whether nanowire is flat, can be
either true or false and -DCELL_ATOMS, describing number of atoms in nanowire
unit cell, can be integer value greater than 0. The program can be executed with
    *./a.out*
command if only pseudopotential file described above is located in program's
directory. Please note that typical program run is very long compared to another
examples and requires a lot of computer resources to finish successfully.

## Known bugs
- (Code readability) The code is autoformatted with specialized tool, but -
unfortunately - not all new C++ features are supported with it and thus some
code formatting bugs are present in main library file and/or examples.
