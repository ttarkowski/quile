#include <cassert>
#include <quile/quile.h>
#include <type_traits>

constexpr const auto d{ quile::uniform_domain<double, 42>(0., 1.) };

using representation = quile::g_floating_point<double, 42, &d>;
using fp_genotype = quile::genotype<representation>;

static_assert(std::is_same_v<fp_genotype::gene_t, double>);
static_assert(std::is_same_v<fp_genotype::genotype_t, representation>);
static_assert(fp_genotype::size() == 42);
static_assert(fp_genotype::constraints() == d);
static_assert(fp_genotype::uniform_domain);

static_assert(quile::is_genotype<fp_genotype>::value);
static_assert(!quile::is_genotype_v<decltype(d)>);

template<typename T>
requires quile::chromosome_representation<T>
struct test_0
{};

template<typename T>
requires quile::chromosome<T>
struct test_1
{};

template<typename T>
requires quile::floating_point_chromosome<T>
struct test_2
{};

int
main()
{
  fp_genotype::chain_t c{ quile::chain_min(d) };
  assert(fp_genotype::valid(c));
  assert(fp_genotype{}.data() == c);
  [[maybe_unused]] test_0<representation> t0{};
  [[maybe_unused]] test_1<fp_genotype> t1{};
  [[maybe_unused]] test_2<fp_genotype> t2{};
}
