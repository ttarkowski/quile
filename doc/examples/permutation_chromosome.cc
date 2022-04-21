#include <algorithm>
#include <cstddef>
#include <iostream>
#include <quile/quile.h>

template<typename G>
requires quile::permutation_chromosome<G> || quile::uniform_chromosome<G>
auto
n_swap_mutation(std::size_t n)
{
  return [=](const G& g) -> quile::population<G> {
    const auto rnd = []() {
      return quile::random_U<std::size_t>(0, G::size() - 1);
    };
    auto res = g.data();
    for (std::size_t i = 0; i < n; ++i) {
      std::swap(res[rnd()], res[rnd()]);
    }
    return quile::population<G>{ G{ res } };
  };
}

int
main()
{
  const std::size_t n = 32;
  quile::genotype<quile::g_permutation<int, n, 1>> g{};
  std::cout << "before n_swap_mutation: " << g.random_reset() << '\n';
  const auto h = n_swap_mutation<decltype(g)>(3)(g)[0];
  std::cout << "after:                  " << h << '\n';
}
