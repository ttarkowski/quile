#include <cassert>
#include <quile/quile.h>

int
main()
{
  const auto f0 = [](int i) { return i == 42; };
  const auto f1 = [](int i) { return i % 2 == 1; };
  const auto f = quile::fn_or(f0, f1);
  assert(f(3));
  assert(f(42));
}
