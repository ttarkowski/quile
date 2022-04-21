#include <cstddef>
#include <iostream>
#include <quile/quile.h>

int
main()
{
  const std::size_t nb = 45;
  quile::genotype<quile::g_binary<nb>> gb{};
  std::cout << "binary:         " << gb.random_reset() << '\n';

  std::cout << '\n';

  const std::size_t nfp = 4;
  static constexpr const auto dfp{ quile::uniform_domain<double, nfp>(-1.,
                                                                      1.) };
  quile::genotype<quile::g_floating_point<double, nfp, &dfp>> gfp{};
  std::cout << "floating-point: " << gfp.random_reset() << '\n';

  std::cout << '\n';

  const std::size_t ni = 32;
  static constexpr const auto di{ quile::uniform_domain<int, ni>(1, 42) };
  quile::genotype<quile::g_integer<int, ni, &di>> gi{};
  std::cout << "integer:        " << gi.random_reset() << '\n';

  std::cout << '\n';

  const std::size_t np = 33;
  quile::genotype<quile::g_permutation<int, np, 0>> gp;
  std::cout << "permutation:    " << gp.random_reset() << '\n';
}
