#include <quile/quile.h>
#include <vector>

using type = quile::domain<double, 42>;

static_assert(quile::is_domain<type>::value);
static_assert(!quile::is_domain_v<std::vector<double>>);

template<typename T>
requires quile::set_of_departure<T>
struct test
{};

int
main()
{
  [[maybe_unused]] test<type> t{};
}
