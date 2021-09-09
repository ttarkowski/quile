#include <algorithm>
#include <ios>
#include <iostream>
#include <quile/quile.h>
#include <string>

using namespace quile;
using namespace quile::test_functions;

namespace {

using type = double;
const std::size_t dim = DIM;

const bool ranking = RANKING;
const bool linear_ranking = LINEAR_RANKING;
const type linear_s = LINEAR_S;

const type eps_f = EPS_F;
const type eps_x = EPS_X;
const type frac = FRAC;

const test_function<type, dim> benchmark_function = TEST_FN;
const auto fn_domain = benchmark_function.function_domain();
const auto p_min = benchmark_function.p_min();
using G = genotype<g_floating_point<type, dim, &fn_domain>>;
static_assert(dim == fn_domain.size());

template<typename T, std::size_t N>
T thinnest_dimension(const domain<T, N>& d) {
  const auto r =
    std::ranges::min(d,
                     std::less{},
                     [](const auto& x) { return x.max() - x.min(); });
  return r.max() - r.min();
}

const type sigma = frac * thinnest_dimension(fn_domain);
const std::size_t generation_sz = GEN_SZ;
const std::size_t parents_sz = PAR_SZ;
static_assert(generation_sz >= parents_sz);

constexpr const type mutation_probability = M_PROB;
constexpr const type recombination_probability = R_PROB;
static_assert(mutation_probability || recombination_probability);

template<typename G>
requires chromosome<G>
point<type, dim> phenotype(const G& g)
{
  point<type, dim> p{};
  std::ranges::copy(g.begin(), g.begin() + dim, p.begin());
  return p;
}

template<typename G>
termination_condition_fn<G>
position_threshold_termination(const point<type, dim>& p_min, type eps_x)
{
  const auto f =
    [=](const G& g) { return distance(p_min, phenotype(g)) <= eps_x; };
  return threshold_termination<G, decltype(f)>(f);
}

template<typename G>
auto
get_selection(const fitness_db<G>& fd)
{
  if constexpr (ranking && linear_ranking) {
    return ranking_selection<G>{ fd, linear_ranking_selection(linear_s) };
  } else if constexpr (ranking && !linear_ranking) {
    return ranking_selection<G>{ fd, exponential_ranking_selection };
  } else {
    return fitness_proportional_selection<G>{ fd };
  }
}

} // anonymous namespace

int main()
{
  const std::size_t max_generations = 100000;
  const variation<G> v{
    stochastic_mutation<G>(Gaussian_mutation<G>(sigma, 1. / dim),
			   mutation_probability),
    stochastic_recombination<G>(arithmetic_recombination<G>,
				recombination_probability)
  };
  const auto f = [](const point<type, dim>& p) {
    return -benchmark_function(p);
  };
  const fitness_function<G> ff = [&](const G& g) { return f(phenotype(g)); };
  const fitness_db<G> fd{ ff, constraints_satisfied<G>, 1 };
  const auto sel{ get_selection<G>(fd) };
  const auto p0 = random_population<constraints_satisfied<G>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ sel };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ sel });
  const fitness tr = f(benchmark_function.p_min());
  const auto tc_1 = fn_and(fitness_threshold_termination<G>(fd, tr, eps_f),
                           position_threshold_termination<G>(p_min, eps_x));
  const auto tc_2 = max_iterations_termination<G>(max_generations);
  const auto tc = fn_or(tc_1, tc_2);
  evolution<G>(v, p0, p1, p2, tc, generation_sz, parents_sz, 1);
  const auto g_min = fd.rank_order()[0];
  const type dist_f = std::fabs(fd(g_min) - tr);
  const type dist_x = distance(p_min, phenotype(g_min));
  const auto res = dist_f <= eps_f && dist_x <= eps_x ? fd.size() : 0;
  std::cout << dim << ' '
            << (!ranking? "FPS" : linear_ranking? "lin-RS" : "exp-RS") << ' '
            << std::scientific << linear_s << ' '
            << std::scientific << eps_f << ' '
            << std::scientific << eps_x << ' '
            << std::scientific << frac << ' '
            << std::scientific << mutation_probability << ' '
            << std::scientific << recombination_probability << ' '
            << generation_sz << ' '
            << parents_sz << ' '
            << benchmark_function.name() << ": ";
  if (res) {
    std::cout << res << ' '
              << std::scientific << dist_f << ' '
              << std::scientific << dist_x << '\n';
  } else {
    std::cout << "FAIL\n";
  }
}
