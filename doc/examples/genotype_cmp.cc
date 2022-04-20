#include <cassert>
#include <quile/quile.h>

int
main()
{
  static const auto d{ quile::uniform_domain<int, 2>(0, 2) };
  using G = quile::genotype<quile::g_integer<int, 2, &d>>;
  const G g0{ { 0, 2 } };
  const G g1{ { 1, 2 } };
  assert(g0 == g0);
  assert(g0 != g1);
  assert(g0 < g1);
  assert(g0 <= g1);
  assert(g1 > g0);
  assert(g1 >= g0);
}
