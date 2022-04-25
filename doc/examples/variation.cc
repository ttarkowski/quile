#include <cassert>
#include <iostream>
#include <quile/quile.h>

int
main()
{
  using namespace quile;
  using G = genotype<g_permutation<int, 7, 0>>;

  const variation<G> v0{};
  const auto g0 = G::random();
  const auto g1 = G::random();
  const population<G> p0 = v0(g0, g1);
  assert(p0[0] == g0 && p0[1] == g1);

  const variation<G> v1{ swap_mutation<G> };
  const auto g2 = G::random();
  const auto g3 = G::random();
  const auto p1 = v1(population<G>{ g0, g1, g2, g3 });
  assert(p1.size() == 4);
  std::cout << p1[0] << '\n';
  std::cout << p1[1] << '\n';
  std::cout << p1[2] << '\n';
  std::cout << p1[3] << '\n';

  const variation<G> v2{ cut_n_crossfill<G> };
  try {
    [[maybe_unused]] const auto p2 = v2(population<G>{ g0, g1, g2 });
  } catch (...) {
    std::cout << "Even number of genotypes is required.\n";
  }

  const variation<G> v3{ swap_mutation<G>, cut_n_crossfill<G> };
  const auto p3 = v3(g0, g1);
  std::cout << p3[0] << '\n';
  std::cout << p3[1] << '\n';
}
