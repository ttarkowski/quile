#include <iostream>
#include <quile/quile.h>

int
main()
{
  using namespace quile;
  using G = genotype<g_permutation<int, 7, 0>>;

  const auto g0 = G::random();
  const auto g1 = G::random();
  std::cout << g0 << '\n';
  std::cout << g1 << '\n';

  const auto h0 = unary_identity(g0)[0];
  std::cout << h0 << '\n';

  const auto p = binary_identity(g0, g1);
  std::cout << p[0] << '\n';
  std::cout << p[1] << '\n';
}
