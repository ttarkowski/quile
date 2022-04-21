#include <iostream>
#include <quile/quile.h>

template<typename G>
requires quile::integer_chromosome<G> quile::population<G>
min_mutation(const G& g)
{
  G res{ g };
  const quile::probability p = .2;
  const auto& c = G::constraints();
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (quile::success(p)) {
      res.value(i, c[i].min());
    }
  }
  return quile::population<G>{ res };
}

int
main()
{
  const std::size_t n = 32;
  static constexpr const auto d{ quile::uniform_domain<int, n>(0, 9) };
  quile::genotype<quile::g_integer<int, n, &d>> g{};
  std::cout << "before min_mutation: " << g.random_reset() << '\n';
  const auto h = min_mutation(g)[0];
  std::cout << "after:               " << h << '\n';
}
