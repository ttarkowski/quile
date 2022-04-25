#include <iostream>
#include <quile/quile.h>

int
main()
{
  using namespace quile;
  const std::size_t n = 32;
  using G = genotype<quile::g_binary<n>>;

  mutation_fn<G> m = stochastic_mutation<G>(swap_mutation<G>, 0.5);
  recombination_fn<G> r = stochastic_recombination<G>(one_point_xover<G>, 0.5);

  const auto g0 = G::random();
  const auto g1 = G::random();
  std::cout << "genotypes: \n" << g0 << '\n' << g1 << '\n';

  const population<G> p0{ m(g0)[0], m(g1)[0] };
  std::cout << "after mutation: \n" << p0[0] << '\n' << p0[1] << '\n';

  const population<G> p1{ r(g0, g1) };
  std::cout << "after recombination: \n" << p1[0] << '\n' << p1[1] << '\n';
}
