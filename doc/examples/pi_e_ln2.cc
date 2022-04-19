#include <iostream>
#include <quile/quile.h>

int
main()
{
  std::cout << "pi = " << quile::pi<double> << '\n';
  std::cout << "e = " << quile::e<double> << '\n';
  std::cout << "ln 2 = " << quile::ln2<double> << '\n';
}
