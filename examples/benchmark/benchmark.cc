// Test functions for floating-point optimization
// - representation: floating-point
// - variation type: Gaussian mutation (fixed sigma), arithmetic recombination
// - parents/surivor selection: roulette wheel selection
// - termination condition: fitness function treshold/fixed number of iterations

#include <cmath>
#include <concepts>
#include <functional>
#include <iostream>
#include <iterator>
#include <numbers>
#include <quile/quile.h>
#include <string>

using namespace quile;
using namespace quile::test_functions;

namespace {

using type = double;
const std::size_t dim = 2;

const test_function<double, dim> benchmark_functions[] = {
  // Unimodal nonseparable
  unimodal::nonseparable::Aluffi_Pentini<type>,
  unimodal::nonseparable::Beale<type>,
  unimodal::nonseparable::Brown<type, dim>,
  unimodal::nonseparable::Matyas<type>,
  unimodal::nonseparable::Rosenbrock<type, dim>,
  unimodal::nonseparable::Schwefel_1<type, dim>,
  unimodal::nonseparable::Schwefel_2<type, dim>,
  // Unimodal separable
  unimodal::separable::Easom<type>,
  unimodal::separable::step<type, dim>,
  unimodal::separable::stepint<type, dim>,
  unimodal::separable::sphere<type, dim>,
  // Multimodal nonseparable
  multimodal::nonseparable::Ackley<type, dim>,
  multimodal::nonseparable::Booth<type>,
  multimodal::nonseparable::Bukin_2<type>,
  multimodal::nonseparable::exponential<type, dim>,
  multimodal::nonseparable::Goldstein_Price<type>,
  multimodal::nonseparable::Himmelblau<type>,
  multimodal::nonseparable::Hosaki<type>,
  multimodal::nonseparable::Leon<type>,
  multimodal::nonseparable::Mexican_hat<type>,
  // Multimodal separable
  multimodal::separable::Alpine<type, dim>
};

template<std::size_t I>
const auto fn_domain = benchmark_functions[I].function_domain();

template<std::size_t I>
using benchmark_genotype = genotype<g_floating_point<type, dim, &fn_domain<I>>>;

template<std::size_t I>
std::size_t
benchmark(const variation<benchmark_genotype<I>>& v,
          std::size_t generation_sz,
          fitness eps,
          std::size_t max_generations)
{
  using G = benchmark_genotype<I>;
  const auto f = [](const point<type, dim>& p) {
    return -benchmark_functions[I](p);
  };
  const fitness_function<G> ff = [&](const G& g) {
    return f(point<type, dim>{ g.value(0), g.value(1) });
  };
  const auto constraints = [](const G&) { return true; };
  const fitness_db<G> fd{ ff, constraints };
  const fitness_proportional_selection<G> fps{ fd };
  const auto p0 = random_population<constraints, G>{};
  const auto p1 = roulette_wheel_selection<G>{ fps };
  const auto p2 = adapter<G>(roulette_wheel_selection<G>{ fps });
  const std::size_t parents_sz{ 2 };
  const fitness tr = f(benchmark_functions[I].p_min());
  const auto tc_1 = fitness_treshold_termination<G>(fd, tr, eps);
  const auto tc_2 = max_iterations_termination<G>(max_generations);
  const auto tc = fn_or(tc_1, tc_2);
  evolution<G>(v, p0, p1, p2, tc, generation_sz, parents_sz, 1);
  return std::fabs(fd(fd.rank_order()[0]) - tr) <= eps ? fd.size() : 0;
}

}

int
main()
{
  const std::size_t generation_sz = 100;
  const fitness eps = 1e-2;
  const std::size_t max_generations = 200000;
  const std::size_t n = 5;
  static_loop<std::size_t, 0, std::size(benchmark_functions)>::body(
    [=](auto I) {
      using G = benchmark_genotype<I>;
      const variation<G> v{ Gaussian_mutation<G>(1e-2, 1. / dim),
                            arithmetic_recombination<G> };
      std::cout << benchmark_functions[I].name() << ": " << std::flush;
      for (std::size_t i = 0; i < n; ++i) {
        const auto res = benchmark<I>(v, generation_sz, eps, max_generations);
        if (res) {
          std::cout << res << " " << std::flush;
        } else {
          std::cout << "FAIL " << std::flush;
        }
      }
      std::cout << '\n';
    });
}
