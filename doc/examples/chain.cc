#include <cassert>
#include <quile/quile.h>

int
main()
{
  const quile::domain<int, 3> d{ quile::range{ 0, 5 },
                                 quile::range{ 1, 6 },
                                 quile::range{ 2, 7 } };
  const quile::chain<int, 3> c0{ 0, 1, 2 };
  const quile::chain<int, 3> c1{ quile::chain_min(d) };
  assert(c0 == c1);
}
