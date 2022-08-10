#include "src/nanotube.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iterator>
#include <quile/quile.h>
#include <string>
#include <vector>

using namespace quile;
using namespace mithril;

namespace {

const std::size_t n_phi = N_PHI;
const std::size_t n_z = N_Z;
const std::size_t c = 2 * n_phi * n_z;
using G = genotype<g_binary<2 * n_phi * n_z>>;

const double decomposition_values[7] = { 0.,     1.7803, 5.1787, 5.6504,
                                         6.2522, 6.5718, 6.5116 };

template<binary_chromosome G>
double
energy(const G& g)
{
  return energy_from_model<G, n_phi, n_z>(g, decomposition_values);
}

const double energy_prec = 1e-4;

[[maybe_unused]] bool
check_precision(double e0, double e1, double prec)
{
  return std::fabs(e0 - e1) <= prec;
}

template<binary_chromosome G>
using abstract_classes = std::vector<G>;

template<chromosome G>
G
deserialize_line(const std::string line)
{
  std::istringstream istream{ line };
  std::size_t i;
  istream >> i;
  G g{};
  istream >> g;
  fitness f{};
  istream >> f;
  assert(check_precision(f, energy(g), energy_prec));
  return g;
}

template<chromosome G>
abstract_classes<G>
deserialize(std::istream& is)
{
  abstract_classes<G> res{};
  for (std::string line{}; std::getline(is, line);) {
    const G g{ deserialize_line<G>(line) };
    const G abstract{ find_min_element_of_abstract_class<G, n_phi, n_z>(g) };
    assert(check_precision(energy(g), energy(abstract), energy_prec));
    if (std::find(std::begin(res), std::end(res), abstract) == res.end()) {
      res.push_back(abstract);
    }
  }
  return res;
}

}

int
main()
{
  std::ifstream ifile{ "evolution.dat" };
  auto res{ deserialize<G>(ifile) };
  std::ranges::stable_sort(res, std::ranges::greater{}, &energy<G>);
  std::ofstream ofile{ "result.dat" };
  for (auto g : res) {
    ofile << g << " : " << number_of_atoms(g) << ' ';
    std::ranges::copy(decomposition<G, n_phi, n_z>(g),
                      std::ostream_iterator<std::size_t>(ofile, " "));
    ofile << d_n_phi<G, n_phi, n_z>(g) << ' ' << d_n_z<G, n_phi, n_z>(g) << ' '
          << std::scientific << std::setprecision(9) << energy(g) << '\n';
  }
}
