#include <iostream>
#include <quile/quile.h>
#include <vector>

using namespace quile;

namespace {

template<chromosome G, typename... Rs>
requires(recombination<Rs, G>&&...) auto random_recombination(Rs&&... rs)
{
  return
    [r = std::vector<recombination_fn<G>>{ rs... }](const G& g0, const G& g1) {
      return r[random_U<std::size_t>(0, r.size() - 1)](g0, g1);
    };
}

}

int
main()
{
  const std::size_t n = 32;
  static constexpr const auto d{ uniform_domain<double, n>(0., 1.) };
  using G = genotype<g_floating_point<double, n, &d>>;

  const auto g0 = G::random();
  const auto g1 = G::random();
  std::cout << "Parents:\n";
  std::cout << g0 << '\n';
  std::cout << g1 << '\n';

  const recombination_fn<G> r = random_recombination<G>(
    single_arithmetic_recombination<G>, one_point_xover<G>);
  const population<G> p = r(g0, g1);
  std::cout << "Offspring:\n";
  std::cout << p[0] << '\n';
  std::cout << p[1] << '\n';
}
