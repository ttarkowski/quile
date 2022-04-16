#include <cassert>
#include <quile/quile.h>

int
main()
{
  assert(quile::square(3) + quile::square(4) == quile::square(5));
}
