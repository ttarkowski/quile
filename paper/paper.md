---
title: 'Quilë: C++ genetic algorithms scientific library'
tags:
  - C++
  - genetic algorithm
  - header-only library
  - scientific library
authors:
  - name: Tomasz Tarkowski
    orcid: 0000-0002-6071-9389
    affiliation: 1
affiliations:
  - name: Chair of Complex Systems Modelling,
          Institute of Theoretical Physics,
          Faculty of Physics,
          University of Warsaw,
          Pasteura 5,
          PL-02093 Warszawa,
          Poland
    index: 1
date: 18 October 2022
bibliography: paper.bib
---

# Summary

This work discusses a general-purpose genetic algorithms [@Holland1975] scientific
header-only library named Quilë. The software is written in C++20 and has been
released under the terms of the MIT license. It is available at
[https://github.com/ttarkowski/quile/](https://github.com/ttarkowski/quile/).
The name of the library come from the fictional language Neo-Quenya and means
"color" (cf. origin of the word *chromosome*).

Genetic algorithms, or more broadly, *evolutionary computations*, is a field
of computer science devoted to population-based, trial-and-error methods of
problem solving. Evolutionary computation was invented by Alan Turing by
noting that the principles of biological evolution and genetics can be applied
to optimization problems [@Turing1948; @Turing1950]. One of the first
evolutionary computations performed on a computer was done by Nils A.
Barricelli [@Barricelli1962; @Galloway2011]. Evolutionary computations have
found wide applications in many disciplines [@Katoch2021; @Ghaheri2015;
@Goudos2016; @Kudjo2017; @Lee2018; @Drachal2021].

The genetic algorithm is a conceptually simple procedure. The aim is to find the
most optimal solution of a given optimization problem. First, one begins with
some sampling of the space of potential solutions; this can be done randomly or
*ad hoc*. This procedure forms the first population of the evolutionary process,
i.e., its first iteration. Each candidate solution from the population is
evaluated in terms of its "fitness". A "fitter" candidate solution has a greater
probability of becoming a parent and entering the next population of the evolution's
subsequent iteration. However, less "fit" individuals still can propagate,
albeit with lower probability. This trait of evolutionary computations helps
preventing premature optimization to the wrong local optimum. Iterations of the
evolution continue as long as the termination condition is not met. One example of
a termination condition is reaching the fitness function *plateau*.

# Overview

*Nature of problem:* Floating-point, integer, binary, or permutation
single-objective constrained finite-dimensional optimization problems of
arbitrary nature, i.e., maximization of $f\colon G \rightarrow \mathbb{R}$, where
$G \subseteq \prod_{i=0}^{c-1} X_i$, where $X_i$ is equal to set of logical
values $\mathbb{B}$ (i.e. false and true) or is bounded subset of set of real
numbers $\mathbb{R}$ or integer numbers $\mathbb{Z}$.

*Solution method:* Single-objective constrained genetic algorithm with pure
floating-point, integer, binary, and permutation representation with set of
exchangeable components (e.g., variation operators, selection probability
functions, selection mechanisms, termination conditions).

*Unique features:* Header-only and easy-to-deploy genetic algorithms C++20
library tailored for problems with computationally expensive fitness functions.

*Functionality:*

- Floating-point, integer, binary, and permutation genotype representations.
- Automatic compile-time representation/variation compatibility checks.
- Mutation operators: Gaussian, self-adaptive, swap, random-reset, and
  bit-flipping.
- Recombination operators: arithmetic, single arithmetic, one-point crossover,
  and cut-and-crossfill.
- Canonical composition of mutation and recombination is available.
- Variations can be applied stochastically.
- Selection mechanisms: stochastic universal sampling, roulette wheel algorithm,
  and generational selection.
- Selection probability functions: fitness proportionate selection with
  windowing procedure, linear and exponential pressure ranking selection.
- Termination conditions: reaching the fitness function *plateau*, reaching the
  maximum number of iterations, reaching the given fitness function value,
  reaching some user defined threshold.
- Conjunction and disjunction of termination conditions are supported.
- Possibility of addition of new variation operators, selection mechanisms,
  selection probability functions, and termination conditions at client-side
  code.

*Limitations:*

- Client-side program compilation time is comparatively long due to
  *header-only* nature of the library and the use of templates.
- Set of variation operators is limited. In case of floating-point
  representation alone, the popular BLX recombination [@Eshelman1993] is not
  implemented.
- Very popular mechanisms of selection to the next generation, e.g.,
  $(\mu + \lambda )$-selection and $(\mu , \lambda )$-selection, are not
  implemented.
- Recombination with number of parents different than 2 is not supported.

*Note*: For information about API documentation, tutorial, installation
instructions, test suite, code statistics, reporting problems with the library,
support inquiries, and feature requests, please see the *README* file in the software
archive.

*Scholarly publications and research projects using the software:*

- research project in computational materials science (cf.
  [Acknowledgements](#acknowledgements) section)
- PhD thesis [@Tarkowski2022b]
- 1 report and 2 articles [@Tarkowski2022a; @Tarkowski2022; @Tarkowski2023]

# Statement of need

Scientific use of C++ genetic algorithms libraries seems to be dominated by four
software packages: GAlib [@Wall], Evolving Objects [@Keijzer2002], OpenBEAGLE
[@Gagne2006], and Evolutionary Computation Framework (ECF) [@Jakobovic]. Evolving
Objects, OpenBEAGLE, and ECF are written in the C++98 standard of the language, while
GAlib is written in a pre-standard version of it. Those libraries, which are
written in C++98, while being comprehensive and feature-rich, also tend to be
relatively hard to use for the novice user inclined toward scientific computation.
The reason is the comparatively complex installation process or complexity of
the library itself.

The Quilë library tries to fill the niche of easy-to-use, high-performance
genetic algorithms scientific libraries by implementing the features using
the modern C++20 standard of the language. The software provides an easy starting point for
researchers and academic teachers who need genetic algorithms and use C++ for
their work. The library is available as a header-only (one file) implementation
that can be installed by simply copying its source code. The user can run
the accompanying examples in matter of mere minutes. On the other hand, the library
also intentionally strives to be minimal, so only a limited set of use cases is
covered.

The library is implemented in generic (template metaprogramming) and partly in
functional programming style with elements of concurrency. Modern elements of
C++ are used, e.g., *concepts* [@Sutton2017], alongside established
constructs, e.g., substitution failure is not an error (SFINAE) [@Vandevoorde2002]. 
In order to compile programs developed with the library, a C++ compiler supporting 
the C++20 standard of the language is needed.

The library employs a database of already-computed fitness function values. The
software is therefore well suited for optimization tasks with fitness functions
calculated from simulation codes like *ab initio* (condensed matter physics) or
finite element method calculations.
This feature, however, does not exclude the use of inexpensive fitness functions.
For the sake of convenience the database itself is available through the
intermediary objects, which are responsible for database cohesion and lifetime
management. The intermediary object type is implemented with the use of the smart
pointer *std::shared_ptr* [@Alexandrescu2001]. Calculations of fitness function values
not yet available in the database  are computed concurrently in order to speed up
the whole evolutionary process; thread pool design pattern is used to optimally
balance the load on CPU cores [@Williams2019].

Example usage of the library features is outlined in the *tutorial/tutorial.txt*
file in the software archive. Example results obtained with the use of the
library are shown in \autoref{fig1} and in \autoref{fig2}.

![Genetic drift over function $f(x, y) = \cos \frac{r(x, y)}{4} + e$, where
$r(x, y) = \sqrt{x^2 + y^2}$, on domain $[-10, 10]^2$. Three genotype
generations were presented: the first one (#0), in half of the evolution (#6)
and at the end (#12). Each generation has the same size. This figure can be
reproduced by running *examples/example_1* example from the software archive.
\label{fig1}](fig1.pdf)

![The $n$ queens puzzle: time of genetic program execution $\Delta t$ and
problem complexity $\xi$. Calculations were performed on low-spec CPU (*thermal
design power* 6 W) for sample size of 10 per measurement point. Number of $n$
queens puzzle solutions can be taken from [@oeisA000170]. Inset: Sample solution
*4Q3/Q7/7Q/5Q2/2Q5/6Q1/1Q6/3Q4 w - - 0 0* for eight queens puzzle reached after
creation of 56 unique genotypes in evolution (on average 123.61 with standard
deviation 111.69 for sample size 100). This figure can be reproduced by running
*examples/example_2/time* example from the software archive.
\label{fig2}](fig2.pdf)

# Acknowledgements

This work is a result of projects funded by the National Science Centre of
Poland (Twardowskiego 16, PL-30312 Kraków, Poland,
[http://www.ncn.gov.pl/](http://www.ncn.gov.pl/)) under grant number
[UMO-2016/23/B/ST3/03575](https://projekty.ncn.gov.pl/en/index.php?projekt_id=359039).

# References
