#include <iostream>
#include <quile/quile.h>

int
main()
{
  const double x = 1.;
  const double y = 0.;
  const double z = 0.;
  std::cout << "(x, y, z) = (" << x << ", " << y << ", " << z << ")\n";

  const auto [r, t, p] = quile::cart2spher(x, y, z);
  std::cout << "(r, theta, phi) = (" << r << ", " << t << ", " << p << ")\n";

  const auto [x2, y2, z2] = quile::spher2cart(r, t, p);
  std::cout << "(x, y, z) = (" << x2 << ", " << y2 << ", " << z2 << ")\n";
}
