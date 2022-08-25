#include "src/nanotube.h"
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <quile/quile.h>

using namespace quile;
using namespace mithril;

namespace {

const std::size_t n_phi = N_PHI;
const std::size_t n_z = N_Z;
const std::size_t c = 2 * n_phi * n_z;
using G = genotype<g_binary<2 * n_phi * n_z>>;

template<binary_chromosome G>
G
read_genotype(std::istream& is)
{
  G g{};
  is >> g;
  return g;
}

const double a = .75; // scaling factor

double
y_coord(std::size_t i)
{
  return (i / n_z) * a * std::sqrt(3.) / 2.;
}

double
x_coord(std::size_t i)
{
  return (i % n_z) * a + (((i / n_z) % 2 == 0) ? 0. : a / 2.);
}

auto
neighbors(std::size_t i)
{
  static const hex_lattice_ord h{ n_phi, n_z };
  return h.neighbors(i);
}

template<binary_chromosome G>
void
latex_input(std::ostream& os, const G& g)
{
  os << "\\documentclass[tikz]{standalone}\n"
     << "\n"
     << "\\pgfdeclarelayer{background}\n"
     << "\\pgfsetlayers{background, main}\n"
     << "\n"
     << "\\begin{document}\n"
     << "  \\begin{tikzpicture}\n";
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (g.value(i)) {
      os << "    \\node[ball color=gray, circle] (" << i << ") at ("
         << std::fixed << std::setprecision(9) << x_coord(i) << ", "
         << std::fixed << std::setprecision(9) << y_coord(i) << ") {};\n";
    }
  }
  os << "    \\begin{pgfonlayer}{background}\n";
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (g.value(i)) {
      for (auto j : neighbors(i)) {
        if (g.value(j)) {
          os << "      \\draw[gray, line width=1mm](" << i << ")--(" << j
             << ");\n";
        }
      }
    }
  }
  os << "    \\end{pgfonlayer}\n"
     << "  \\end{tikzpicture}\n"
     << "\\end{document}\n";
}

}

int
main()
{
  latex_input(std::cout, read_genotype<G>(std::cin));
}
