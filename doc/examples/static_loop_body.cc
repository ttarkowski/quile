#include <iostream>
#include <quile/quile.h>

template<int N>
struct test
{
  constexpr static int value = N;
};

int
main()
{
  quile::static_loop<int, 0, 3>::body([=](auto I) {
    using T = test<I>;
    std::cout << T::value << std::endl;
  });
}
