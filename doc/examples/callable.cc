#include <cassert>
#include <quile/quile.h>

template<quile::callable<bool> B>
bool
negate(B b)
{
  return !b();
}

int
main()
{
  assert(negate([]() { return false; }));
}
