#include <cassert>
#include <iostream>
#include <quile/quile.h>

int
main()
{
  const std::size_t n = 32;
  static constexpr const auto d{ quile::uniform_domain<int, n>(0, 9) };
  using G = quile::genotype<quile::g_integer<int, n, &d>>;
  const quile::populate_0_fn<G> generator =
    quile::random_population<quile::constraints_satisfied<G>, G>;
  const quile::population<G> p = generator(10);
  assert(p.size() == 10);
  for (auto& g : p) {
    std::cout << g << '\n';
  }
}
