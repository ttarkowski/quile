#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const auto d{ quile::uniform_domain<int, 42>(0, 41) };

using p_type = quile::g_permutation<int, 42, 0>;

static_assert(std::is_same_v<p_type::type, int>);
static_assert(p_type::size() == 42);
static_assert(p_type::constraints() == d);

static_assert(quile::is_g_permutation<p_type>::value);
static_assert(!quile::is_g_permutation_v<decltype(d)>);

template<typename T>
requires quile::permutation_representation<T>
struct test
{};

int
main()
{
  p_type::chain_t c{ quile::iota<int, 42>(0) };
  assert(p_type::valid(c));
  assert(p_type::default_chain() == c);
  [[maybe_unused]] test<p_type> t{};
}
