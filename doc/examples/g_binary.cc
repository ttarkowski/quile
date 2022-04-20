#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const quile::domain<bool, 42> d{};

using i_type = quile::g_binary<42>;

static_assert(std::is_same_v<i_type::type, bool>);
static_assert(i_type::size() == 42);
static_assert(i_type::constraints() == d);

static_assert(quile::is_g_binary<i_type>::value);
static_assert(!quile::is_g_binary_v<decltype(d)>);

template<typename T>
requires quile::binary_representation<T>
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
