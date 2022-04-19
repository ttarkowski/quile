#include <cmath>
#include <cstddef>
#include <fstream>
#include <quile/quile.h>

using namespace quile;

using type = double;
const std::size_t dim = 2;

fitness
f(type x, type y)
{
  return -(x * x + y * y);
}

int
main()
{
  static const domain<type, dim> d0{ range{ -35., +35. }, range{ -35., +35. } };
  static const domain<type, 2 * dim> d{ self_adaptive_variation_domain(d0,
                                                                       .001) };
  using G = genotype<g_floating_point<type, 2 * dim, &d>>;
  const fitness_function<G> ff = [](const G& g) {
    return f(g.value(0), g.value(1));
  };
  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  const auto p0 = random_population<constraints_satisfied<G>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ fps };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ fps });

  const std::size_t generation_sz{ 1000 };
  const std::size_t parents_sz{ 42 };
  const auto tc = max_fitness_improvement_termination<G>(fd, 10, 0.05);

  const variation<G> v{ self_adaptive_mutation<G>(.002, .002),
                        arithmetic_recombination<G> };

  std::ofstream file{ "evolution.dat" };
  print(file, evolution<G>(v, p0, p1, p2, tc, generation_sz, parents_sz));
}
