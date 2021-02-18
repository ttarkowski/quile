// Evolutionary search for a maximum of given function over domain
// - function: f(x, y) = cos(0.25 * r(x, y)) + e
// - domain: [-10, +10] x [-10, +10]
// - representation: floating-point
// - variation type: Gaussian mutation (fixed sigma), no recombination
// - parents/surivor selection: stochastic universal sampling (SUS)
// - termination condition: based on maximum fitness improvement

#include <cmath>
#include <cstddef>
#include <fstream>
#include <numbers>
#include <quile/quile.h>

using namespace quile;

using type = double;

// Function
const auto f = [](type x, type y) -> fitness {
  const auto r = [](type x, type y) -> type {
    return std::sqrt(x * x + y * y);
  };
  return std::cos(0.25 * r(x, y)) + std::numbers::e_v<fitness>;
};

// Domain
const domain<type, 2> d{ range{ -10., +10. }, range{ -10., +10. } };

int
main()
{
  using G = genotype<g_floating_point<type, 2, &d>>;
  const fitness_function<G> ff = [](const G& g) {
    return f(g.value(0), g.value(1));
  };
  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  // First generation creator
  const auto p0 = random_population<constraints_satisfied<G>, G>;
  // Parents selection
  const auto p1 = stochastic_universal_sampling<G>{ fps };
  // Survivor selection
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ fps });

  const std::size_t generation_sz{ 1000 };
  const std::size_t parents_sz{ 42 };
  const auto tc = max_fitness_improvement_termination<G>(fd, 10, 0.05);

  const type sigma{ .2 };
  const variation<G> v{ Gaussian_mutation<G>(sigma, 1.) };

  std::ofstream file{ "evolution.dat" };
  for (std::size_t i = 0; const auto& x : evolution<G>(
                            v, p0, p1, p2, tc, generation_sz, parents_sz)) {
    for (const auto& xx : x) {
      file << i << ' ' << xx << '\n';
    }
    ++i;
  }
}
