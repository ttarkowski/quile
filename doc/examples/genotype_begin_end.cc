#include <iostream>
#include <quile/quile.h>

using p_type = quile::g_permutation<int, 42, 0>;

int
main()
{
  using G = quile::genotype<p_type>;
  const G g{};
  for (G::const_iterator it = g.begin(); it != g.end(); ++it) {
    std::cout << *it << ' ';
  }
  std::cout << '\n';
  for (auto x : g) {
    std::cout << x << ' ';
  }
  std::cout << '\n';
}
