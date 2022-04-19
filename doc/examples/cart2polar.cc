#include <iostream>
#include <quile/quile.h>

int
main()
{
  const double x = 0.;
  const double y = 1.;
  std::cout << "(x, y) = (" << x << ", " << y << ")\n";

  const auto [r, p] = quile::cart2polar(x, y);
  std::cout << "(r, phi) = (" << r << ", " << p << ")\n";

  const auto [x2, y2] = quile::polar2cart(r, p);
  std::cout << "(x, y) = (" << x2 << ", " << y2 << ")\n";
}
