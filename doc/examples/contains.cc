#include <cassert>
#include <quile/quile.h>

int
main()
{
  quile::domain<double, 2> d{ quile::range{ 0., 1. }, quile::range{ 0., 1. } };
  std::array<double, 2> p{ .5, .5 };
  std::array<double, 2> q{ 2., 3. };
  assert(contains(d, p));
  assert(!contains(d, q));
}
