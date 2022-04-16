#include <iostream>
#include <quile/quile.h>
#include <tuple>

using result = std::tuple<int, int>;

result
draw(quile::probability p, int sz)
{
  int t{ 0 };
  int f{ 0 };
  for (int i = 0; i < sz; ++i) {
    if (quile::success(p)) {
      ++t;
    } else {
      ++f;
    }
  }
  return { t, f };
}

void
print(std::ostream& os, quile::probability p, const result& r)
{
  const auto f = [&](char c, int sz) {
    os << c << ": ";
    for (int i = 0; i < sz; ++i) {
      os << '*';
    }
    os << '\n';
  };
  os << "probability: " << p << '\n';
  f('t', std::get<0>(r));
  f('f', std::get<1>(r));
  os << '\n';
}

int
main()
{
  const int n = 40;
  for (auto p : { 0., .25, .5, .75, 1. }) {
    print(std::cout, p, draw(p, n));
  }
}
