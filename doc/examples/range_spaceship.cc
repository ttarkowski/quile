#include <cassert>
#include <quile/quile.h>

int
main()
{
  quile::range r0{ 0, 2 };
  quile::range r1{ 1, 2 };
  assert(r0 < r1);
  assert(r0 <= r1);
  assert(r1 > r0);
  assert(r1 >= r0);
}
