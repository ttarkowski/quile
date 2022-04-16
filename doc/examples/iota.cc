#include <iostream>
#include <quile/quile.h>

int
main()
{
  for (auto x : quile::iota<int, 7>(0)) {
    std::cout << x << ' ';
  }
  std::cout << '\n';
}
