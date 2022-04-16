#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <quile/quile.h>

using result = std::map<int, int>;

result
draw(int min, int max, int sz)
{
  result res{};
  for (int i = 0; i < sz; ++i) {
    ++res[std::floor(quile::random_U(min, max))];
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
  print(std::cout, fill_gaps(draw(-10, 10, 200)));
}
