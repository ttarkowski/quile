#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const quile::domain<bool, 42> d{};

using b_type = quile::g_binary<42>;

static_assert(std::is_same_v<b_type::type, bool>);
static_assert(b_type::size() == 42);
static_assert(b_type::constraints() == d);

static_assert(quile::is_g_binary<b_type>::value);
static_assert(!quile::is_g_binary_v<decltype(d)>);

template<typename T>
requires quile::binary_representation<T>
struct test
{};

int
main()
{
  b_type::chain_t c{ quile::chain_min(d) };
  assert(b_type::valid(c));
  assert(b_type::default_chain() == c);
  [[maybe_unused]] test<b_type> t{};
}
