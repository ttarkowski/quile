#include <iostream>
#include <quile/quile.h>
#include <random>

int
main()
{
  const int n = 100;
  double res = 0;
  for (int i = 0; i < n; ++i) {
    res += std::bernoulli_distribution{ .5 }(quile::random_engine()) ? +1 : -1;
  }
  std::cout << res / n << '\n';
}
