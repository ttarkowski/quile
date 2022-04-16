#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <quile/quile.h>

using result = std::map<int, int>;

result
draw(double m, double sd, int sz)
{
  result res{};
  for (int i = 0; i < sz; ++i) {
    ++res[static_cast<int>(std::floor(quile::random_N(m, sd)))];
  }
  return res;
}

result
fill_gaps(result res)
{
  for (int i = res.begin()->first; i <= res.rbegin()->first; ++i) {
    res[i];
  }
  return res;
}

void
print(std::ostream& os, const result& r)
{
  const auto f = [&](int n, int sz) {
    os << std::setw(3) << n << ": ";
    for (int i = 0; i < sz; ++i) {
      os << '*';
    }
    os << '\n';
  };
  for (auto [n, sz] : r) {
    f(n, sz);
  }
}

int
main()
{
  print(std::cout, fill_gaps(draw(0., 5., 200)));
}
