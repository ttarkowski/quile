#include <cassert>
#include <quile/quile.h>

int
main()
{
  long x = -80538738812075974;
  long y = 80435758145817515;
  long z = 12602123297335631;
  long k = 42;
  assert(quile::cube(x) + quile::cube(y) + quile::cube(z) == k);
}
