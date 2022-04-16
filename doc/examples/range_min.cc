#include <cassert>
#include <quile/quile.h>

int
main()
{
  quile::range r{ 0, 42 };
  assert(r.min() == 0);
}
