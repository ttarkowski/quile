#include <quile/quile.h>

template<typename T, typename G, T t>
requires quile::genotype_constraints<T, G>
struct test
{};

int
main()
{
  using G = quile::genotype<quile::g_binary<42>>;
  [[maybe_unused]] test<decltype(quile::constraints_satisfied<G>),
                        G,
                        quile::constraints_satisfied<G>>
    t{};
}
