#include <iostream>
#include <quile/quile.h>

const std::size_t n = 32;
using genotype_t = quile::genotype<quile::g_binary<n>>;
using population_t = quile::population<genotype_t>;

static_assert(quile::is_population<population_t>::value);
static_assert(!quile::is_population_v<genotype_t>);

template<typename T>
requires quile::genetic_pool<T>
struct test
{};

int
main()
{
  [[maybe_unused]] test<population_t> t{};
  population_t p{};
  for (int i = 0; i < 8; ++i) {
    p.push_back(genotype_t{});
  }
  for (auto& g : p) {
    g.random_reset();
  }
  for (auto& g : p) {
    std::cout << g << '\n';
  }
}
