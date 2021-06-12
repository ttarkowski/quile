// Crystal structure prediction of boron nanowires
// - representation: binary
// - variation type: random-reset mutation, no recombination
// - parents/surivor selection: stochastic universal sampling (SUS)
// - termination condition: based on maximum fitness improvement

#include "../common/pwx.h"
#include "../common/system.h"
#include "../common/ts_unordered_map.h"
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <quile/quile.h>
#include <string>

using namespace quile;

namespace {

template<typename G>
thread_safe_unordered_map<G, std::string> file_db{};

const double bond_length = 1.675; // Angstrom
const std::size_t m_size = M_SIZE;
const std::size_t n_size = N_SIZE;
const pwx_atom atom{ "B", 10.811, "B.pbe-n-kjpaw_psl.1.0.0.UPF" };
using G = genotype<g_binary<m_size * n_size>>;

template<binary_chromosome G>
pwx_positions
geometry(const G& g)
{
  pwx_positions res{};
  for (std::size_t i = 0; auto x : g) {
    if (x) {
      const double x = ((i % m_size) + (i / m_size) * .5) * bond_length;
      const double y = (i / m_size) * bond_length * std::sqrt(3.) / 2.;
      res.push_back({ atom.symbol, x, y, 0. });
    }
    ++i;
  }
  return res;
}

template<binary_chromosome G>
std::size_t
number_of_atoms(const G& g)
{
  // TODO: std::ranges::count_if(g, true) does not work.
  return std::ranges::count_if(g, [](bool b) -> bool { return b; });
}

template<binary_chromosome G>
void
input_file(const std::string& filename, const G& g)
{
  const int k_points = 8;
  const double dz = 15.;
  file_db<G>.insert_or_modify(g, filename);
  std::ofstream file{ filename };
  file << pwx_control("scf", filename)
       << pwx_system(0, number_of_atoms<G>(g), 1, 0., 1.e-2, 6.e+1)
       << pwx_electrons(100, 7.e-1)
       << pwx_cell_parameters_triangle_60deg(bond_length, m_size, n_size, dz)
       << pwx_atomic_species({ atom }) << pwx_atomic_positions(geometry<G>(g))
       << pwx_k_points(k_points, k_points, 1, 1, 1, 0);
}

}

int
main()
{
  const auto ff = []<binary_chromosome G>(const G& g) -> fitness {
    const std::string input_filename{ pwx_unique_filename() };
    input_file<G>(input_filename, g);
    const auto [o, e] = execute("/bin/bash calc.sh " + input_filename);
    return o == "Calculations failed.\n" ? incalculable : -std::stod(o);
  };

  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  const auto p0 = random_population<constraints_satisfied<G>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ fps };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ fps });

  const std::size_t generation_sz{ 100 };
  const std::size_t parents_sz{ 24 };
  const auto tc_1 = max_fitness_improvement_termination<G>(fd, 10, 0.05);
  const auto tc_2 = max_iterations_termination<G>(1000);
  const auto tc = fn_or(tc_1, tc_2);

  const variation<G> v{ random_reset<G>(1. / (m_size * n_size)) };

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
