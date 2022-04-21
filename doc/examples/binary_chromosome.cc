#include <iostream>
#include <quile/quile.h>

template<typename G>
requires quile::binary_chromosome<G>
auto
min_mutation(quile::probability p)
{
  return [=](const G& g) -> quile::population<G> {
    G res{ g };
    const auto& c = G::constraints();
    for (std::size_t i = 0; i < G::size(); ++i) {
      if (quile::success(p)) {
        res.value(i, c[i].min());
      }
    }
    return quile::population<G>{ res };
  };
}

int
main()
{
  const std::size_t n = 32;
  quile::genotype<quile::g_binary<n>> g{};
  std::cout << "before min_mutation: " << g.random_reset() << '\n';
  const auto h = min_mutation<decltype(g)>(.2)(g)[0];
  std::cout << "after:               " << h << '\n';
}
