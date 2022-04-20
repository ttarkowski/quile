#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const auto d{ quile::uniform_domain<int, 100>(0, 42) };

using i_type = quile::g_integer<int, 100, &d>;

static_assert(std::is_same_v<i_type::type, int>);
static_assert(i_type::size() == 100);
static_assert(i_type::constraints() == d);

static_assert(quile::is_g_integer<i_type>::value);
static_assert(!quile::is_g_integer_v<decltype(d)>);

template<typename T>
requires quile::integer_representation<T>
struct test
{};

int
main()
{
  i_type::chain_t c{ quile::chain_min(d) };
  assert(i_type::valid(c));
  assert(i_type::default_chain() == c);
  [[maybe_unused]] test<i_type> t{};
}
