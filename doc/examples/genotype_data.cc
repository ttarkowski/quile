#include <exception>
#include <iostream>
#include <quile/quile.h>
#include <type_traits>

using p_type = quile::g_permutation<int, 42, 0>;
using c_type = quile::chain<int, 42>;

static_assert(std::is_same_v<c_type, p_type::chain_t>);

int
main()
{
  using G = quile::genotype<p_type>;
  const G g0{};
  c_type c{ g0.data() };
  for (auto& x : c) {
    x = 0;
  }
  try {
    const G g1{ c };
  } catch (std::exception& e) {
    std::cout << e.what() << '\n';
  }
}
