#include <iostream>
#include <quile/quile.h>

int
main()
{
  static const auto d{ quile::uniform_domain<int, 5>(0, 9) };
  std::cout << quile::genotype<quile::g_integer<int, 5, &d>>::random() << '\n';
}
