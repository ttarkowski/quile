#include <iostream>
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
                         [](fitness acc, auto x) { return acc - x * x; });
}

int
main()
{
  const std::size_t n = 8;
  static constexpr const auto d{ uniform_domain<int, n>(0, 9) };
  using G = genotype<quile::g_integer<int, n, &d>>;
  const fitness_db<G> fd{ fitness_fn<G>, constraints_satisfied<G> };
  const auto p = random_population<constraints_satisfied<G>, G>(3);

  const selection_probabilities_fn<G> sp_fns[] = {
    fitness_proportional_selection<G>{ fd },
    ranking_selection<G>{ fd, linear_ranking_selection(2.) },
    ranking_selection<G>{ fd, exponential_ranking_selection }
  };

  for (auto sp_fn : sp_fns) {
    const selection_probabilities sp = sp_fn(p);
    const selection_probabilities cp = cumulative_probabilities(sp_fn, p);
    for (std::size_t i = 0; i < p.size(); ++i) {
      std::cout << p[i] << ": " << sp[i] << ", " << cp[i] << '\n';
    }
  }
}
