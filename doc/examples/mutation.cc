#include <iostream>
#include <quile/quile.h>

using namespace quile;

namespace {

const auto id = [](auto x) { return x; };

auto
compose()
{
  return id;
}

template<typename F, typename... Fs>
auto
compose(F&& f, Fs&&... fs)
{
  return [=](auto x) { return f(compose(fs...)(x)); };
}

template<chromosome G, typename... Ms>
requires(mutation<Ms, G>&&...) auto mutation_composition(Ms&&... ms)
{
  return [=](const G& g) -> population<G> {
    return { compose([=](const G& g) -> G { return ms(g)[0]; }...)(g) };
  };
}

template<typename G>
requires integer_chromosome<G>
auto
deterministic_mutation(std::size_t pos, typename G::gene_t val)
{
  return [=](G g) -> population<G> {
    g.value(pos, val);
    return population<G>{ g };
  };
}

}

int
main()
{
  const std::size_t n = 5;
  static constexpr const auto d{ uniform_domain<int, n>(0, 9) };
  using G = genotype<g_integer<int, n, &d>>;

  const G g{};

  mutation_fn<G> m0 = mutation_composition<G>();
  std::cout << m0(g)[0] << '\n';

  mutation_fn<G> m1 = mutation_composition<G>(deterministic_mutation<G>(1, 1));
  std::cout << m1(g)[0] << '\n';

  mutation_fn<G> m2 = mutation_composition<G>(deterministic_mutation<G>(2, 2));
  std::cout << m2(g)[0] << '\n';

  mutation_fn<G> m12 = mutation_composition<G>(deterministic_mutation<G>(1, 1),
                                               deterministic_mutation<G>(2, 2));
  std::cout << m12(g)[0] << '\n';

  mutation_fn<G> m134 =
    mutation_composition<G>(deterministic_mutation<G>(1, 1),
                            deterministic_mutation<G>(3, 3),
                            deterministic_mutation<G>(4, 4));
  std::cout << m134(g)[0] << '\n';
}
