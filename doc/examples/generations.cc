#include <cassert>
#include <iostream>
#include <quile/quile.h>

int
main()
{
  const std::size_t n = 4;
  static constexpr const auto d{ quile::uniform_domain<int, n>(0, 9) };
  using G = quile::genotype<quile::g_integer<int, n, &d>>;
  const quile::populate_0_fn<G> generator =
    quile::random_population<quile::constraints_satisfied<G>, G>;
  quile::generations<G> gs{};
  for (int i = 0; i < 10; ++i) {
    gs.push_back(generator(5));
  }
  assert(gs.size() == 10);
  for (int i = 0; auto& p : gs) {
    std::cout << '#' << i++ << ": ";
    for (auto& g : p) {
      std::cout << g << "  ";
    }
    std::cout << '\n';
  }
}
