// Crystal structure prediction of boron nanotubes
// - representation: binary
// - variation type: bit-flipping mutation, one-point crossover recombination
// - parents/surivor selection: stochastic universal sampling (SUS)
// - termination condition: based on maximum fitness improvement

#include "src/nanotube.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <quile/quile.h>

using namespace quile;
using namespace mithril;

namespace {

const std::size_t n_phi = N_PHI;
const std::size_t n_z = N_Z;
const std::size_t c = 2 * n_phi * n_z;
using G = genotype<g_binary<2 * n_phi * n_z>>;

template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
bool
nanotube_condition(const G& g)
{
  return atoms_connected_in_unit_cell<G, n_phi, n_z>(g) &&
         adjacency_at_unit_cell_boundary_along_nanotube<G, n_phi, n_z>(g) &&
         adjacency_at_unit_cell_boundary_at_circumference<G, n_phi, n_z>(g);
}

const double decomposition_values[7] = { 0.,     1.7803, 5.1787, 5.6504,
                                         6.2522, 6.5718, 6.5116 };

}

int
main()
{
  const auto ff = []<binary_chromosome G>(const G& g) -> fitness {
    return energy_from_model<G, n_phi, n_z>(g, decomposition_values);
  };

  const fitness_db<G> fd{ ff, nanotube_condition<G, n_phi, n_z>, 1 };
  const ranking_selection<G> rs{ fd, linear_ranking_selection(2.) };

  const auto p0 = random_population<nanotube_condition<G, n_phi, n_z>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ rs };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ rs });

  const std::size_t generation_sz{ 2 * c };
  const std::size_t parents_sz{ c };
  const fitness dE{ 1e-3 };
  const auto tc = max_fitness_improvement_termination_2<G>(fd, 100, dE);

  const mutation_fn<G> m{ bit_flipping<G>(1. / G::size()) };
  const recombination_fn<G> r{ one_point_xover<G> };
  const variation<G> v{ stochastic_mutation<G>(m, .5),
                        stochastic_recombination<G>(r, .5) };

  std::ofstream file{ "evolution.dat" };
  for (std::size_t i = 0; const auto& x : evolution<G>(
                            v, p0, p1, p2, tc, generation_sz, parents_sz)) {
    for (const auto& xx : x) {
      file << i << ' ' << xx << ' ' << std::scientific << std::setprecision(9)
           << fd(xx) << '\n';
    }
    ++i;
  }
}
