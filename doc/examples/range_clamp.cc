#include <cassert>
#include <quile/quile.h>

int
main()
{
  quile::range<int> r{ 0, 42 };
  assert(r.clamp(-10) == 0);
  assert(r.clamp(0) == 0);
  assert(r.clamp(7) == 7);
  assert(r.clamp(42) == 42);
  assert(r.clamp(100) == 42);
}
