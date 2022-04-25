#include <cassert>
#include <cmath>
#include <cstddef>
#include <quile/quile.h>

int
main()
{
  using namespace quile;
  using type = double;
  const std::size_t dim = 2;
  static const auto d = uniform_domain<type, dim>(-10., +10.);
  using G = genotype<g_floating_point<type, dim, &d>>;
  const fitness_function<G> ff = [](const G& g) {
    return -(std::fabs(g.value(0)) + std::fabs(g.value(1)));
  };
  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  for (int i = 0; i < 2; ++i) {
    const auto g = G::random();
    std::cout << "Fitness function value for " << g << " is " << fd(g) << '.'
              << std::endl;
  }
  const population<G> p{ G::random(), G::random() };
  const fitnesses fs = fd(p);
  for (std::size_t i = 0; i < fs.size(); ++i) {
    std::cout << "fitness function value for " << p[i] << " is " << fs[i] << '.'
              << std::endl;
    assert(fs[i] == fd(p[i]));
  }
  std::cout << "Database size is equal to " << fd.size() << '.' << std::endl;
  std::cout << "Database content:" << std::endl;
  for (auto x : fd) {
    std::cout << ' ' << x.first << " " << x.second << std::endl;
  }
  std::cout << "Database content (once again):" << std::endl;
  for (fitness_db<G>::const_iterator it = fd.begin(); it != fd.end(); ++it) {
    std::cout << ' ' << it->first << " " << it->second << std::endl;
  }
  std::cout << "The best genotype is " << fd.rank_order()[0] << '.'
            << std::endl;
}
