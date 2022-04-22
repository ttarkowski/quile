#include <cassert>
#include <cmath>
#include <numeric>
#include <quile/quile.h>

using namespace quile;

template<chromosome G>
fitness
fitness_fn(const G& g)
{
  return std::accumulate(g.data().begin(),
                         g.data().end(),
                         fitness{ 0. },
                         [](fitness acc, auto x) { return acc + std::abs(x); });
}

int
main()
{
  const std::size_t n = 32;
  static constexpr const auto d{ uniform_domain<int, n>(0, 9) };
  using G = genotype<quile::g_integer<int, n, &d>>;

  const fitness_db<G> fd{ fitness_fn<G>, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  const populate_0_fn<G> generator =
    random_population<constraints_satisfied<G>, G>;

  const populate_1_fn<G> parents_selection =
    stochastic_universal_sampling<G>{ fps };

  const populate_2_fn<G> selection_to_the_next_generation =
    adapter<G>(stochastic_universal_sampling<G>{ fps });

  const population<G> p0 = generator(42);
  const population<G> p1 = generator(42);

  const population<G> q0 = parents_selection(10, p0);
  assert(q0.size() == 10);

  const population<G> q1 = selection_to_the_next_generation(42, p0, p1);
  assert(q1.size() == 42);
}
