# Quilë
Quilë---C++20 header-only genetic algorithms library released under the terms
of the MIT License.

## Etymology
> qya. *quilë* -- color

The part *chromo-* of the word *chromosome*, which is central term in genetics,
comes from the Greek χρῶμα and means *color*. Similarly, the name of this
library origins from fictional language Quenya and means *color* as well.

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

## Dependencies
- Modern compiler with at least partial support of C++20 (library is known to
  work with GCC 10.2.0)

*Please note that code in examples/ directory might require some additional
software.*

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
