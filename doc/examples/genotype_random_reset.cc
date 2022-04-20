#include <iostream>
#include <quile/quile.h>

int
main()
{
  static const auto d{ quile::uniform_domain<int, 5>(0, 9) };
  using G = quile::genotype<quile::g_integer<int, 5, &d>>;
  G g{};
  std::cout << g << '\n';
  std::cout << g.random_reset(2) << '\n';
  std::cout << g.random_reset() << '\n';
}
