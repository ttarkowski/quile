/*
 * MIT License
 *
 * Copyright (c) 2020 Tomasz Tarkowski
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef QUILE_H
#define QUILE_H

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace quile {

  /////////////////
  // Thread pool //
  /////////////////

  class thread_pool {
  public:
    explicit thread_pool(std::size_t sz) : free_threads_{sz} {}

    template<typename T>
    std::future<T> async(std::launch policy, const std::function<T()>& f) {
      return std::async(policy, [this, f]() {
                          acquire();
                          if constexpr (std::is_same_v<T, void>) {
                            f();
                            release();
                          } else {
                            const T res = f();
                            release();
                            return res;
                          }
                        });
    }

  private:
    inline void acquire() {
      std::unique_lock<std::mutex> ul{m_};
      cv_.wait(ul, [this]() { return free_threads_ != 0; });
      --free_threads_;
    }
    
    inline void release() {
      std::unique_lock<std::mutex> ul{m_};
      ++free_threads_;
      cv_.notify_one();
    }

  private:
    std::mutex m_{};
    std::condition_variable cv_{};
    std::size_t free_threads_;
  };
  
  //////////////////////
  // Callable concept //
  //////////////////////
  
  template<typename F, typename R, typename... Args>
  concept callable = std::convertible_to<std::invoke_result_t<F, Args...>, R>;

  ///////////
  // Range //
  ///////////
  
  template<typename T>
  class range {
  public:
    constexpr range(T min, T max) : min_{min}, max_{max} {
      if (min > max) {
        throw std::invalid_argument{"range: min greater than max"};
      }
    }

    template<
      typename U = T,
      typename = std::enable_if_t<std::numeric_limits<U>::is_specialized>>
    constexpr range()
      : range{std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max()}
    {}
    
    constexpr range(const range&) = default;
    constexpr range(range&&) = default;
    range& operator=(const range&) = default;
    range& operator=(range&&) = default;
    T min() const { return min_; }
    T max() const { return max_; }

    template<typename U = T,
             typename = std::enable_if_t<!std::is_same_v<U, bool>>>
    T midpoint() const { return std::midpoint(min_, max_); }
    
    template<typename U = T,
             typename = std::enable_if_t<!std::is_same_v<U, bool>>>
    T clamp(T t) const { return std::clamp(t, min_, max_); }
    
    bool contains(T t) const { return t >= min_ && t <= max_; }
    auto operator<=>(const range<T>& r) const = default;

  private:
    T min_;
    T max_;
  };

  template<typename T>
  std::ostream& operator<<(std::ostream& os, const range<T>& r)
  { return (os << '[' << r.min() << ", " << r.max() << ']'); }
  
  template<typename T, std::size_t N>
  std::array<T, N> iota(T t) {
    std::array<T, N> res{};
    std::ranges::transform(std::views::iota(t) | std::views::take(N),
                           std::begin(res),
                           std::identity{});
    return res;
  }
    
  ////////////////////
  // Random numbers //
  ////////////////////

  using probability = double;

  inline std::mt19937& random_engine() {
    // Only one global hidden variable engine (regardless of number of
    // translation units).
    static std::mt19937 engine{std::random_device{}()};
    return engine;
  }
  
  inline bool success(probability success_probability) {
    return std::bernoulli_distribution{ success_probability }(random_engine());
  }

  template<std::floating_point T>
  T random_from_normal_distribution(T mean, T standard_deviation) {
    auto& generator{random_engine()};
    return std::normal_distribution<T>{mean, standard_deviation}(generator);
  }

  template<std::floating_point T>
  T normal(T mean, T standard_deviation) {
    return random_from_normal_distribution<T>(mean, standard_deviation);
  }

  template<typename T>
  T random_from_uniform_distribution(T a, T b) {
    auto& generator{random_engine()};
    if constexpr (std::is_floating_point_v<T>) {
      assert(a < b);
      return // Check whether b - a overflows [N4861, 26.6.8.2.2].
        a > 0. || b <= std::numeric_limits<T>::max() + a
          ? std::uniform_real_distribution<T>{a, b}(generator) // [a, b)
          : random_from_uniform_distribution(false, true)
            ? random_from_uniform_distribution<T>(a, std::midpoint(a, b))
            : random_from_uniform_distribution<T>(std::midpoint(a, b), b);
    } else if constexpr (std::is_same_v<T, bool>) {
      return a == b ? a : std::bernoulli_distribution{0.5}(generator);
    } else { // [a, b]
      assert(a <= b);
      return std::uniform_int_distribution<T>{a, b}(generator);
    }
  }

  template<typename T>
  T uniform(T a, T b) {
    return random_from_uniform_distribution(a, b);
  }
  
  template<typename T>
  T random(const range<T>& r)
  { return uniform<T>(r.min(), r.max()); }

  ////////////
  // Domain //
  ////////////
  
  template<typename T, std::size_t N>
  using domain = std::array<range<T>, N>;
  
  template<typename T, std::size_t N>
  bool contains(const domain<T, N>& d, const std::array<T, N>& p) {
    bool res = true;
    for (std::size_t i = 0; res && i < N; ++i) {
      res &= d[i].contains(p[i]);
    }
    return res;
  }

  template<typename T, std::size_t N>
  constexpr domain<T, N> uniform_domain(const range<T>& r) {
    domain<T, N> res{};
    std::generate_n(std::begin(res), N, [&]() { return r; });
    return res;
  }

  template<typename T, std::size_t N>
  constexpr bool uniform(const domain<T, N>& d) {
    const auto x0 = N == 0? range<T>{} : d[0];
    return std::ranges::all_of(d, [&](const auto& x) { return x0 == x; });
  }
  
  template<typename T, std::size_t N>
  using chain = std::array<T, N>;
  
  template<typename T, std::size_t N>
  chain<T, N> chain_min(const domain<T, N>& d) {
    chain<T, N> res{};
    std::ranges::transform(d, std::begin(res), std::identity{}, &range<T>::min);
    return res;
  }
  
  //////////////
  // Genotype //
  //////////////
  
  template<typename T, std::size_t N, const domain<T, N>* D>
  requires std::floating_point<T>
  struct g_floating_point {
    static_assert(D != nullptr);
    static_assert(N > 0);
    using type = T;
    static constexpr std::size_t size() { return N; }
    static constexpr const domain<type, size()>& constraints() { return *D; }
    using chain_t = chain<type, size()>;
    
    static bool valid(const chain<type, size()>& c)
    { return contains(constraints(), c); }

    static chain_t default_chain() { return chain_min(constraints()); }
  };
  
  template<typename T>
  struct is_g_floating_point : std::false_type {};
  
  template<typename T, std::size_t N, const domain<T, N>* D>
  struct is_g_floating_point<g_floating_point<T, N, D>> : std::true_type {};
  
  template<typename T>
  inline constexpr bool is_g_floating_point_v = is_g_floating_point<T>::value;
  
  template<typename T>
  concept floating_point_representation = is_g_floating_point_v<T>;
  
  template<typename T, std::size_t N, const domain<T, N>* D>
  requires std::integral<T> && (!std::is_same_v<T, bool>)
  struct g_integer {
    static_assert(D != nullptr);
    static_assert(N > 0);
    using type = T;
    static constexpr std::size_t size() { return N; }
    static constexpr const domain<type, size()>& constraints() { return *D; }
    using chain_t = chain<type, size()>;
    
    static bool valid(const chain<type, size()>& c)
    { return contains(constraints(), c); }
    
    static chain_t default_chain() { return chain_min(constraints()); }
  };
  
  template<typename T>
  struct is_g_integer : std::false_type {};
  
  template<typename T, std::size_t N, const domain<T, N>* D>
  struct is_g_integer<g_integer<T, N, D>> : std::true_type {};
  
  template<typename T>
  inline constexpr bool is_g_integer_v = is_g_integer<T>::value;
  
  template<typename T>
  concept integer_representation = is_g_integer_v<T>;
  
  template<std::size_t N>
  struct g_binary {
    static_assert(N > 0);
    using type = bool;
    static constexpr std::size_t size() { return N; }
    static constexpr const domain<type, size()> constraints()
    { return domain<type, size()>{}; }
    
    using chain_t = chain<type, size()>;
    
    static bool valid(const chain<type, size()>& c) { return true; }
    static chain_t default_chain() { return chain_min(constraints()); }
  };
  
  template<typename T>
  struct is_g_binary : std::false_type {};
  
  template<std::size_t N>
  struct is_g_binary<g_binary<N>> : std::true_type {};
  
  template<typename T>
  inline constexpr bool is_g_binary_v = is_g_binary<T>::value;
  
  template<typename T>
  concept binary_representation = is_g_binary_v<T>;
  
  template<typename T, std::size_t N, T M>
  requires std::integral<T> && (!std::is_same_v<T, bool>)
  struct g_permutation {
    using type = T;
    static constexpr std::size_t size() { return N; }
    
    static constexpr const domain<type, size()> constraints()
    { return uniform_domain<type, size()>(range<type>{M, M + N - 1}); }
    
    using chain_t = chain<type, size()>;
    
    static bool valid(const chain<type, size()>& c) {
      const auto i = iota<type, size()>(M);
      return contains(constraints(), c)
        && std::is_permutation(std::begin(c), std::end(c), std::begin(i));
    }
    
    static chain_t default_chain() { return iota<type, size()>(M); }
  };
  
  template<typename T>
  struct is_g_permutation : std::false_type {};
  
  template<typename T, std::size_t N, T M>
  struct is_g_permutation<g_permutation<T, N, M>> : std::true_type {};
  
  template<typename T>
  inline constexpr bool is_g_permutation_v = is_g_permutation<T>::value;
  
  template<typename T>
  concept permutation_representation = is_g_permutation_v<T>;
  
  template<typename T>
  concept chromosome_representation = floating_point_representation<T>
                                            || integer_representation<T>
                                            || binary_representation<T>
                                            || permutation_representation<T>;
  
  template<typename R> requires chromosome_representation<R>
  class genotype {
  public:
    using chain_t = chain<typename R::type, R::size()>;
    using const_iterator = typename chain_t::const_iterator;
    using gene_t = R::type;
    using genotype_t = R;
    static constexpr std::size_t size() { return R::size(); }
    
    static constexpr const domain<gene_t, size()> constraints()
    { return R::constraints(); }
    
    static constexpr bool uniform_domain = uniform(constraints());
    static bool valid(const chain_t& c) { return R::valid(c); }
    
  public:
    genotype() : chain_{R::default_chain()} {}
    
    explicit genotype(const chain_t& c)
      : chain_{c} {
      if (!valid(c)) {
        throw std::invalid_argument{"invalid chain"};
      }
    }
    
    genotype(const genotype&) = default;
    genotype(genotype&&) = default;
    genotype& operator=(const genotype&) = default;
    genotype& operator=(genotype&&) = default;
    
    gene_t value(std::size_t i) const { return chain_[i]; }
    
    template<typename = std::enable_if_t<!permutation_representation<R>>>
    genotype& value(std::size_t i, gene_t v) {
      if (!constraints()[i].contains(v)) {
        throw std::invalid_argument{"bad value"};
      }
      chain_[i] = v;
      return *this;
    }
    
    genotype& random_reset() {
      if constexpr (permutation_representation<R>) {
        std::shuffle(chain_.begin(), chain_.end(), random_engine());
      } else {
        for (std::size_t i = 0; i < size(); ++i) {
          random_reset(i);
        }
      }
      return *this;
    }
    
    template<typename = std::enable_if_t<!permutation_representation<R>>>
    genotype& random_reset(std::size_t i) {
      const auto& c = constraints();
      chain_[i] = uniform<gene_t>(c[i].min(), c[i].max());
      return *this;
    }
    
    auto operator<=>(const genotype& g) const { return chain_ <=> g.chain_; }
    bool operator==(const genotype& g) const { return chain_ == g.chain_; }
    
    const chain_t& data() const { return chain_; }
    const_iterator begin() const { return chain_.begin(); }
    const_iterator end() const { return chain_.end(); }
    
  private:
    chain_t chain_;
  };
  
  template<typename T>
  struct is_genotype : std::false_type {};
  
  template<typename T>
  struct is_genotype<genotype<T>> : std::true_type {};
  
  template<typename T>
  inline constexpr bool is_genotype_v = is_genotype<T>::value;
  
  template<typename G>
  concept chromosome = is_genotype_v<G>;
  
  template<typename G>
  concept floating_point_chromosome = chromosome<G>
    && floating_point_representation<typename G::genotype_t>;

  template<typename G>
  concept integer_chromosome = chromosome<G>
    && integer_representation<typename G::genotype_t>;
  
  template<typename G>
  concept binary_chromosome = chromosome<G>
    && binary_representation<typename G::genotype_t>;
  
  template<typename G>
  concept permutation_chromosome = chromosome<G>
    && permutation_representation<typename G::genotype_t>;

  template<typename G>
  concept uniform_chromosome = chromosome<G> && G::uniform_domain;
  
  template<typename F, typename G>
  concept genotype_constraints = std::predicate<F, G> && chromosome<G>;
  
  template<typename G> requires chromosome<G>
  const auto constraints_satisfied = [](const G&) { return true; };
  
  template<typename G> requires chromosome<G>
  std::ostream& operator<<(std::ostream& os, const G& g) {
    os << '[';
    for (std::size_t i = 0; i < G::size(); ++i) {
      os << g.value(i) << i + 1 < G::size()? ", " : ']';
    }
    return os;
  }
  
} // namespace quile

template<typename G> requires quile::chromosome<G>
struct std::hash<G> {
  std::size_t operator()(const G& g) const noexcept {
    const std::size_t sz{sizeof(std::size_t) * CHAR_BIT};
    std::size_t res{0};
    for (std::size_t i = 0; i < g.size(); ++i) {
      res ^= std::hash<typename G::gene_t>{}(g.value(i)) << i % sz;
    }
    return res;
  }
};
  
namespace quile {

  ////////////////
  // Population //
  ////////////////
  
  template<typename G> requires chromosome<G>
  using population = std::vector<G>;
  
  template<typename T>
  struct is_population : std::false_type {};
  
  template<typename G>
  struct is_population<population<G>> : std::true_type {};
  
  template<typename G>
  inline constexpr bool is_population_v = is_population<G>::value;
  
  template<typename G>
  concept genetic_pool = is_population_v<G>;
  
  // Population generators/selectors
  // - first generation creator
  template<typename G> requires chromosome<G>
  using populate_0_fn = std::function<population<G>(std::size_t)>;
  // - parents selection
  template<typename G> requires chromosome<G>
  using populate_1_fn = std::function<population<G>(std::size_t,
                                                    const population<G>&)>;
  // - survivor selection
  template<typename G> requires chromosome<G>
  using populate_2_fn = std::function<population<G>(std::size_t,
                                                    const population<G>&,
                                                    const population<G>&)>;
  
  template<typename G> requires chromosome<G>
  using generations = std::deque<population<G>>;
  
  //////////////////////////////
  // Mutation & recombination //
  //////////////////////////////
  
  template<typename M, typename G>
  concept mutation = requires (M m, G g)
    {
     { m(g) } -> std::convertible_to<population<G>>;
    } && chromosome<G>;
  
  template<typename G> requires chromosome<G>
  using mutation_fn = std::function<population<G>(const G&)>;
  
  template<typename R, typename G>
  concept recombination = requires (R r, G g)
    {
     { r(g, g) } -> std::convertible_to<population<G>>;
    } && chromosome<G>;
  
  template<typename G> requires chromosome<G>
  using recombination_fn = std::function<population<G>(const G&, const G&)>;
  
  template<typename G> requires chromosome<G>
  population<G> unary_identity(const G& g) { return population<G>{g}; }
  
  template<typename G> requires chromosome<G>
  population<G> binary_identity(const G& g0, const G& g1)
  { return population<G>{g0, g1}; }

  template<typename G> requires chromosome<G>
  class variation {
  public:
    variation(const mutation_fn<G>& m, const recombination_fn<G>& r)
      : m_{m}, r_{r}
    {}

    variation()
      : variation{unary_identity<G>, binary_identity<G>}
    {}

    explicit variation(const mutation_fn<G>& m)
      : variation{m, binary_identity<G>}
    {}

    explicit variation(const recombination_fn<G>& r)
      : variation{unary_identity<G>, r}
    {}

    population<G> operator()(const G& g0, const G& g1) const {
      population<G> res{};
      for (const auto& g : r_(g0, g1)) {
        res.push_back(m_(g).at(0));
      }
      assert(res.size() == 1 || res.size() == 2);
      return res;
    }
    
    population<G> operator()(const population<G>& p) const {
      if (p.size() % 2) {
        throw std::invalid_argument{"wrong population size"};
      }
      population<G> res;
      for (std::size_t i = 0; i < p.size(); i += 2) {
        for (const auto& g : this->operator()(p[i], p[i + 1])) {
          res.push_back(g);
        }
      }
      assert(res.size() == p.size() / 2 || res.size() == p.size());
      return res;
    }
    
  private:
    mutation_fn<G> m_;
    recombination_fn<G> r_;
  };

  template<typename G> requires chromosome<G>
  auto stochastic_mutation(const mutation_fn<G>& m, probability p)
  { return [=](const G& g) { return success(p)? m(g) : population<G>{g}; }; }

  template<typename G> requires chromosome<G>
  auto stochastic_recombination(const recombination_fn<G>& r, probability p) {
    return
      [=](const G& g0, const G& g1) {
        const auto tmp = r(g0, g1);
        return success(p)
          ? tmp
          : tmp.size() == 2
            ? population<G>{g0, g1}
            : success(.5)? population<G>{g0} : population<G>{g1};
      };
  }
  
  ///////////////
  // Evolution //
  ///////////////
  
  template<typename F, typename G>
  concept termination_condition = std::predicate<F,
                                                 std::size_t, generations<G>>;
  
  template<typename G> requires chromosome<G>
  using termination_condition_fn = std::function<bool(std::size_t,
                                                      const generations<G>&)>;
  
  // TODO: Is there any way to reduce number of arguments of this function
  // without increasing solution's complexity?
  template<typename G> requires chromosome<G>
  generations<G> evolution(const variation<G> v,
                           const population<G>& first_generation,
                           const populate_1_fn<G>& p1,
                           const populate_2_fn<G>& p2,
                           const termination_condition_fn<G>& tc,
                           std::size_t parents_sz,
                           std::size_t max_history = 0) {
    generations<G> res{};
    const std::size_t generation_sz = first_generation.size();
    for (std::size_t i = 0; !tc(i, res); ++i) {
      const population<G> p{i == 0
                            ? first_generation
                            : p2(generation_sz,
                                 res.back(),
                                 v(p1(parents_sz, res.back())))};
      res.push_back(p);
      if (max_history && res.size() > max_history) {
        res.pop_front();
      }
    }
    return res;
  }
  
  template<typename G> requires chromosome<G>
  generations<G> evolution(const variation<G>& v,
                           const populate_0_fn<G>& p0,
                           const populate_1_fn<G>& p1,
                           const populate_2_fn<G>& p2,
                           const termination_condition_fn<G>& tc,
                           std::size_t generation_sz,
                           std::size_t parents_sz,
                           std::size_t max_history = 0) {
    return
      evolution<G>(v, p0(generation_sz), p1, p2, tc, parents_sz, max_history);
  }
  
  //////////////////////
  // Fitness function //
  //////////////////////
  
  using fitness = double;
  using fitnesses = std::vector<fitness>;

  template<typename G> requires chromosome<G>
  using fitness_function = std::function<fitness(const G&)>;
  
  const fitness incalculable = -std::numeric_limits<fitness>::infinity();
  
  template<typename G> requires chromosome<G>
  class fitness_db {
  public:
    explicit fitness_db(const fitness_function<G>& f,
                        const genotype_constraints<G> auto& gc,
                        unsigned int thread_sz =
                          std::thread::hardware_concurrency())
      : function_{[=](const G& g) { return gc(g)? f(g) : incalculable; }}
      , thread_sz_{thread_sz}
    {}
    
    fitness_db(const fitness_db&) = default;
    fitness_db& operator=(const fitness_db&) = default;
    
    fitness operator()(const G& g) const {
      const auto it{ fitness_values_->find(g) };
      return it != fitness_values_->end()
        ? it->second
        : ((*fitness_values_)[g] = function_(g));
    }
    
    fitnesses operator()(const population<G>& p) const {
      if (thread_sz_ > 1 && p.size() > 1) {
        multithreaded_calculations(p);
      }
      fitnesses res{};
      std::ranges::transform(p, std::back_inserter(res),
                             [this](const G& g) { return operator()(g); });
      return res;
    }
    
    std::size_t size() const { return fitness_values_->size(); }
    
  private:
    auto uncalculated_fitnesses(const population<G>& p) const {
      std::unordered_set<G> res{};
      std::ranges::copy_if(p, std::inserter(res, std::end(res)),
                           [this](const G& g) {
                             return !fitness_values_->contains(g);
                           });
      return res;
    }
    
    void multithreaded_calculations(const population<G>& p) const {
      using type = std::pair<G, fitness>;
      thread_pool tp{thread_sz_};
      std::vector<std::future<type>> v{};
      for (const auto& x : uncalculated_fitnesses(p)) {
        v.push_back(tp.async<type>(std::launch::async,
                                   [this, x]() {
                                     const fitness xf = this->function_(x);
                                     return type{x, xf};
                                   }));
      }
      for (auto& x : v) {
        fitness_values_->insert(x.get());
      }
    }
    
  private:
    fitness_function<G> function_;
    unsigned int thread_sz_;
    std::shared_ptr<std::unordered_map<G, fitness>> fitness_values_ =
      std::make_shared<std::unordered_map<G, fitness>>();
  };
  
  /////////////////////////////
  // Selection probabilities //
  /////////////////////////////
  
  using selection_probabilities = std::vector<probability>;
  
  template<typename G> requires chromosome<G>
  using selection_probabilities_fn =
    std::function<selection_probabilities(const population<G>&)>;
  
  template<typename G> requires chromosome<G>
  selection_probabilities
  cumulative_probabilities(const selection_probabilities_fn<G>& spf,
                           const population<G>& p) {
    auto res = spf(p);
    std::partial_sum(res.begin(), res.end(), res.begin());
    // Last element should be exactly equal to 1. and another part of algorithm
    // might require this exact identity. Unfortunately, numerical calculations
    // might not be so precise. Let's check if last element is calculated with
    // 1% precision (basic requirement):
    assert(res.back() > .99 && res.back() < 1.01);
    // Then, let's correct the value:
    res.back() = 1.;
    return res;
  }
  
  template<typename C, typename T>
  C select_different_than(const C& c,
                          T t,
                          bool require_nonempty_result) {
    C res{};
    std::ranges::copy_if(c, std::back_inserter(res),
                         [=](auto x) { return x != t; });
    if (require_nonempty_result && res.size() == 0) {
      throw std::runtime_error{"empty result"};
    }
    return res;
  }
  
  inline fitnesses select_calculable(const fitnesses& fs,
                                     bool require_nonempty_result = false) {
    return select_different_than(fs,
                                 incalculable,
                                 require_nonempty_result);
  }
  
  template<typename G>
  class fitness_proportional_selection {
  public:
    explicit fitness_proportional_selection(const fitness_db<G>& ff)
      : ff_{ff}
    {}

    // FPS with windowing with workarounds for:
    // a) population of equally fit genotypes and
    // b) populations containing genotypes which fitnesses cannot be calculated
    // Please note that in b) case, there should be at least one genotype,
    // which fitness can be calculated.
    selection_probabilities operator()(const population<G>& p) const {
      const fitnesses fs{ff_(p)};
      const auto cal = select_calculable(fs, true);
      const fitness min = *std::ranges::min_element(cal);
      const auto n = cal.size();
      const fitness delta = 1. / n;
      const fitness sum =
        std::accumulate(std::begin(cal), std::end(cal), fitness{0.})
        - n * min + 1;
      selection_probabilities res{};
      std::ranges::transform(fs, std::back_inserter(res),
                             [=](fitness f) {
                               return
                                 f == incalculable
                                   ? .0
                                   : (f - min + delta) / sum;
                             });
      return res;
    }
    
  private:
    const fitness_db<G> ff_;
  };
  
  ////////////////////////////
  // Generators & selectors //
  ////////////////////////////
  
  namespace detail {
    
    template<typename G> requires chromosome<G>
    population<G> generate(std::size_t lambda, const std::function<G()>& f) {
      population<G> res{};
      std::generate_n(std::back_inserter(res), lambda, f);
      return res;
    }

  } // namespace detail

  template<typename G> requires chromosome<G>
  populate_2_fn<G> adapter(const populate_1_fn<G>& fn) {
    return [=](std::size_t sz,
               const population<G>& p0,
               const population<G>& p1) {
             population<G> p{p0};
             p.insert(p.end(), p1.begin(), p1.end());
             return fn(sz, p);
           };
  }
  
  template<auto C, typename G>
  requires genotype_constraints<decltype(C), G> && chromosome<G>
  class random_population {
  public:
    random_population() = default;
    
    explicit random_population(unsigned int thread_sz)
      : thread_sz_{thread_sz}
    {}

    population<G> operator()(std::size_t lambda) const {
      thread_pool tp{thread_sz_};
      std::vector<std::future<G>> v{};
      for (std::size_t i = 0; i < lambda; ++i) {
        v.push_back(tp.async<G>(std::launch::async,
                                [g = G{}, this]() mutable -> G {
                                  while (!C(g.random_reset()));
                                  return g;
                                }));
      }
      population<G> res{};
      for (auto& x : v) {
        res.push_back(x.get());
      }
      return res;
    }
    
  private:
    unsigned int thread_sz_ = std::thread::hardware_concurrency();
  };
  
  template<typename G> requires chromosome<G>
  class roulette_wheel_selection {
  public:
    explicit roulette_wheel_selection(const selection_probabilities_fn<G>& spf)
      : spf_{spf}
    {}
    
    population<G> operator()(std::size_t lambda, const population<G>& p) const {
      const auto f = [&, c = cumulative_probabilities(spf_, p)]() -> G {
          return p.at(std::distance(c.begin(),
                                    std::lower_bound(c.begin(),
                                                     c.end(),
                                                     uniform<double>(0., 1.))));
        };
      return detail::generate<G>(lambda, f);
    }
    
  private:
    const selection_probabilities_fn<G> spf_;
  };

  template<typename G> requires chromosome<G>
  class stochastic_universal_sampling {
  public:
    explicit
    stochastic_universal_sampling(const selection_probabilities_fn<G>& spf)
      : spf_{spf}
    {}
    
    population<G> operator()(std::size_t lambda, const population<G>& p) const {
      const auto a = cumulative_probabilities(spf_, p);
      auto r = uniform<double>(0., 1. / lambda);
      
      population<G> res{};
      for (std::size_t i = 0, j = 0; j < lambda; ++i) {
        for (; r <= a.at(i) && j < lambda; r += 1. / lambda, ++j) {
          res.push_back(p.at(i));
        }
      }
      std::shuffle(res.begin(), res.end(), random_engine());
      return res;
    }
    
  private:
    const selection_probabilities_fn<G> spf_;
  };

  template<typename G> requires chromosome<G>
  population<G>
  generational_survivor_selection(std::size_t sz,
                                  const population<G>& generation,
                                  const population<G>& offspring) {
    if (generation.size() != sz || offspring.size() != sz) {
      throw std::invalid_argument{"bad size"};
    }
    return offspring;
  }

  ///////////////////////////
  // Termination condition //
  ///////////////////////////
  
  inline fitness max(const fitnesses& fs) {
    const fitnesses calc{select_calculable(fs, true)};
    return *std::ranges::max_element(calc);
  }
  
  template<typename G> requires chromosome<G>
  fitness max(const population<G>& p, const fitness_db<G>& ff) {
    return max(ff(p));
  }
  
  template<typename G> requires chromosome<G>
  fitnesses max(const generations<G>& gs, const fitness_db<G>& ff) {
    fitnesses res{};
    std::ranges::transform(gs, std::back_inserter(res),
                           [&ff](const population<G>& p) { return max(p, ff); });
    return res;
  }
  
  inline fitness min(const fitnesses& fs) {
    const fitnesses calc{select_calculable(fs, true)};
    return *std::ranges::min_element(calc);
  }
  
  template<typename G> requires chromosome<G>
  fitness min(const population<G>& p, const fitness_db<G>& ff) {
    return min(ff(p));
  }
  
  template<typename G> requires chromosome<G>
  fitnesses min(const generations<G>& gs, const fitness_db<G>& ff) {
    fitnesses res{};
    std::ranges::transform(gs, std::back_inserter(res),
                           [&ff](const population<G>& p) { return min(p, ff); });
    return res;
  }
  
  template<typename G>
  termination_condition_fn<G> max_iterations_termination(std::size_t max) {
    return [=](std::size_t i, const generations<G>&) { return i == max; };
  }
  
  template<typename G>
  termination_condition_fn<G>
  max_fitness_improvement_termination(const fitness_db<G>& ff,
                                      std::size_t n,
                                      double frac) {
    return
      [=]([[maybe_unused]] std::size_t i, const generations<G>& gs) {
        assert(i == gs.size());
        if (gs.size() <= n) {
          return false;
        } else {
          const fitnesses fs{max(gs, ff)};
          const fitness min_last_n =
            *std::min_element(fs.end() - n, fs.end());
          const double x = (max(fs) - min_last_n) / (max(fs) - min(fs));
          return x <= frac;
        }
      };
  }

  template<typename G>
  termination_condition_fn<G>
  fitness_treshold_termination(const fitness_db<G>& fd,
                               fitness tr,
                               fitness eps) {
    return
      [=](std::size_t, const generations<G>& gs) {
        return gs.empty()
          ? false
          : std::ranges::any_of(gs.back(),
                                [=, &fd](const G& g) {
                                  return std::fabs(fd(g) - tr) <= eps;
                                });
      };
  }
  
  /////////////////////////////////////////////////
  // Concrete mutation & recombination operators //
  /////////////////////////////////////////////////
  
  template<typename G> requires floating_point_chromosome<G>
  auto Gaussian_mutation(typename G::gene_t sigma, probability p) {
    return [=](const G& g) -> population<G> {
      G res{};
      const auto c = G::constraints();
      for (std::size_t i = 0; i < G::size(); ++i) {
        if (success(p)) {
          res.value(i, c[i].clamp(g.value(i) + sigma * normal(0., 1.)));
        }
      }
      return population<G>{res};
    };
  }
  
  template<typename G> requires uniform_chromosome<G>
  population<G> swap_mutation(const G& g) {
    const std::size_t n = G::size();
    auto d = g.data();
    std::swap(d[uniform<std::size_t>(0, n - 1)],
              d[uniform<std::size_t>(0, n - 1)]);
    return population<G>{G{d}};
  }
  
  template<typename G> requires chromosome<G>
  auto random_reset(probability p) {
    return
      [=](const G& g) -> population<G> {
      G res{g};
      for (std::size_t i = 0; i < G::size(); ++i) {
        if (success(p)) {
          res.random_reset();
        }
      }
      return population<G>{res};
    };
  }
  
  template<typename G> requires floating_point_chromosome<G>
  population<G> arithmetic_recombination(const G& g0, const G& g1) {
    G res{};
    for (std::size_t i = 0; i < G::size(); ++i) {
      res.value(i, std::midpoint(g0.value(i), g1.value(i)));
    }
    return population<G>{res};
  }

  template<typename G> requires uniform_chromosome<G>
  population<G> one_point_xover(const G& g0, const G& g1) {
    auto d0 = g0.data();
    auto d1 = g1.data();
    const std::size_t n = G::size();
    const auto cp = uniform<std::size_t>(0, n - 1);
    for (std::size_t i = cp; i < n; ++i) {
      std::swap(d0[i], d1[i]);
    }
    return population<G>{G{d0}, G{d1}};
  }

  template<typename G> requires permutation_chromosome<G>
  population<G> cut_n_crossfill(const G& g0, const G& g1) {
    const auto f =
      [cp = uniform<std::size_t>(1, G::size() - 1)](const G& g, auto d) {
        auto it = std::begin(d);
        std::advance(it, cp);
        for (auto x : g) {
          if (std::find(std::begin(d), it, x) == it) {
            *it++ = x;
          }
        }
        assert(it == std::end(d));
        return d;
      };
    return population<G>{G{f(g1, g0.data())}, G{f(g0, g1.data())}};
  }
  
} // namespace quile

#endif // QUILE_H
