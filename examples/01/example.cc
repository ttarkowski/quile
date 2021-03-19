// Evolutionary search for a maximum of given function over domain
// - function: f(x) = sin(2 * x) * exp(-0.05 * x^2) + pi
// - domain: [-10, +10]
// - representation: floating-point
// - variation type: no mutation, arithmetic recombination
// - parents/surivor selection: roulette wheel selection
// - termination condition: fixed number of iterations

#include <cmath>
#include <cstddef>
#include <fstream>
#include <quile/quile.h>

using namespace quile;

using type = double;
const std::size_t dim = 1;

// Function
fitness
f(type x)
{
  return std::sin(2 * x) * std::exp(-0.05 * x * x) + pi<fitness>;
}

// Domain
const domain<type, dim> d{ range{ -10., +10. } };

int
main()
{
  using G = genotype<g_floating_point<type, dim, &d>>;
  const fitness_function<G> ff = [](const G& g) { return f(g.value(0)); };
  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  // First generation creator
  const auto p0 = random_population<constraints_satisfied<G>, G>;
  // Parents selection
  const auto p1 = roulette_wheel_selection<G>{ fps };
  // Survivor selection
  const auto p2 = adapter<G>(roulette_wheel_selection<G>{ fps });

  const std::size_t generation_sz{ 1000 };
  const std::size_t parents_sz{ 42 };
  const auto tc = max_iterations_termination<G>(100);

  const variation<G> v{ arithmetic_recombination<G> };

  std::ofstream file{ "evolution.dat" };
  for (std::size_t i = 0; const auto& x : evolution<G>(
                            v, p0, p1, p2, tc, generation_sz, parents_sz)) {
    for (const auto& xx : x) {
      file << i << ' ' << xx << '\n';
    }
    ++i;
  }
}
