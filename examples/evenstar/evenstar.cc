// Crystal structure prediction of boron nanowires
// - representation: floating-point
// - variation type: Gaussian mutation (fixed sigma), arithmetic recombination
// - parents/surivor selection: stochastic universal sampling (SUS)
// - termination condition: based on maximum fitness improvement

#include "../common/pwx.h"
#include "../common/system.h"
#include "../common/ts_unordered_map.h"
#include "src/nanowire.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <quile/quile.h>
#include <string>

using namespace quile;
using namespace evenstar;

namespace {

template<typename G>
thread_safe_unordered_map<G, std::string> file_db{};

using type = double;
const bool flat = FLAT;
const std::size_t cell_atoms = CELL_ATOMS;
const range<type> bond_range{ 1.54, 2.10 }; // Angstrom
const pwx_atom atom{ "B", 10.811, "B.pbe-n-kjpaw_psl.1.0.0.UPF" };
const auto d = construct_domain<type, cell_atoms, flat>(bond_range);
using G = genotype<g_floating_point<type, std::size(d), &d>>;

template<typename G>
requires floating_point_chromosome<G>
void
input_file(const std::string& filename, const G& g)
{
  const int k_points = 8;
  file_db<G>.insert_or_modify(g, filename);
  std::ofstream file{ filename };
  const auto [p, h] = geometry<G>(g, atom.symbol, flat);
  const auto max_x = std::ranges::max_element(p, {}, &pwx_position::x)->x;
  const auto max_y = std::ranges::max_element(p, {}, &pwx_position::y)->y;
  const typename G::gene_t free_space = 10.;
  file << pwx_control("scf", filename)
       << pwx_system(0, number_of_atoms<G>(flat), 1, 0., 1.e-2, 6.e+1)
       << pwx_electrons(100, 7.e-1)
       << pwx_cell_parameters_diag(max_x + free_space, max_y + free_space, h)
       << pwx_atomic_species({ atom }) << pwx_atomic_positions(p)
       << pwx_k_points(1, 1, k_points, 0, 0, 1);
  assert(number_of_atoms<G>(flat) == p.size());
}

template<quile::floating_point_chromosome G>
bool
nanowire_condition(const G& g)
{
  const auto ps = geometry_pbc<G>(g, atom.symbol, flat);
  return atoms_not_too_close(ps, bond_range.min()) &&
         all_atoms_connected(ps, bond_range.max());
}

}

int
main()
{
  const auto ff = []<floating_point_chromosome G>(const G& g) -> fitness {
    const std::string input_filename{ pwx_unique_filename() };
    input_file<G>(input_filename, g);
    const auto [o, e] = execute("/bin/bash calc.sh " + input_filename);
    return o == "Calculations failed.\n" ? incalculable : -std::stod(o);
  };

  const fitness_db<G> fd{ ff, nanowire_condition<G> };
  const fitness_proportional_selection<G> fps{ fd };

  const auto p0 = random_population<nanowire_condition<G>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ fps };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ fps });

  const std::size_t generation_sz{ 1000 };
  const std::size_t parents_sz{ 42 };
  const auto tc_1 = max_fitness_improvement_termination<G>(fd, 10, 0.05);
  const auto tc_2 = max_iterations_termination<G>(1000);
  const auto tc = fn_or(tc_1, tc_2);

  const type sigma{ .02 };
  const variation<G> v{ Gaussian_mutation<G>(sigma, 1.),
                        arithmetic_recombination<G> };

  std::ofstream file{ "evolution.dat" };
  for (std::size_t i = 0; const auto& x : evolution<G>(
                            v, p0, p1, p2, tc, generation_sz, parents_sz)) {
    for (const auto& xx : x) {
      file << i << ' ' << xx << ' ' << std::scientific << std::setprecision(9)
           << fd(xx) << ' ' << file_db<G>.at(xx) << '\n';
    }
    ++i;
  }
}
