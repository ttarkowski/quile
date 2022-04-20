#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const auto d{ quile::uniform_domain<double, 42>(0., 1.) };

using fp_type = quile::g_floating_point<double, 42, &d>;

static_assert(std::is_same_v<fp_type::type, double>);
static_assert(fp_type::size() == 42);
static_assert(fp_type::constraints() == d);

static_assert(quile::is_g_floating_point<fp_type>::value);
static_assert(!quile::is_g_floating_point_v<decltype(d)>);

template<typename T>
requires quile::floating_point_representation<T>
struct test
{};

int
main()
{
  fp_type::chain_t c{ quile::chain_min(d) };
  assert(fp_type::valid(c));
  assert(fp_type::default_chain() == c);
  [[maybe_unused]] test<fp_type> t{};
}
