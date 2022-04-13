/*
 * MIT License
 *
 * Copyright (c) 2020, 2021 Tomasz Tarkowski
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
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <numeric>
#include <random>
#include <ranges>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

/**
 * \file
 * `quile/quile.h` is one and only file of the library implementation.
 *
 * Please include `quile/quile.h` header using the code `#include
 * <quile/quile.h>` and compiling your program with apropriate include flag,
 * e.g. `-I/home/user/repos/quile` for GCC and Clang compilers.
 *
 * Example compilation command: `g++ -std=c++20 -DNDEBUG -O3 -Wall -Wextra
 * -pedantic -I/home/user/repos/quile -pthread program.cc`.
 *
 * Please remove `-DNDEBUG` to enable assertions.
 *
 * Please add `-DQUILE_ENABLE_LOGGING` to enable logging.
 *
 * Clang compilation command is similar: please change `g++` to `clang++`.
 */

/**
 * \def QUILE_LOG(x)
 * `QUILE_LOG(x)` macro prints diagnostic information.
 *
 * `QUILE_LOG(x)` macro prints diagnostic information to the standard error
 * stream `std::cerr` provided that `QUILE_ENABLE_LOGGING` token is defined.
 * Otherwise it has no effect.
 *
 * Example:
 * \include QUILE_LOG.cc
 *
 * Result:
 * \verbinclude QUILE_LOG.out
 */

#ifdef QUILE_ENABLE_LOGGING
#define QUILE_LOG(x)                                                           \
  do {                                                                         \
    std::cerr << "# Quile log: " << x << '\n';                                 \
  } while (0)
#else
#define QUILE_LOG(x)
#endif

namespace quile {

//////////////
// TMP loop //
//////////////

/**
 * `static_loop` is a `for`-loop replacement with loop index usable at compile
 * time.
 *
 * \tparam T Loop index type.
 * \tparam I First index value (inclusive).
 * \tparam N Last index value (exclusive).
 */
template<std::integral T, T I, T N>
struct static_loop
{
  /**
   * `static_loop::body` performs `f` in a loop iterating from `I` (inclusively)
   * to `N` (exclusively), i.e. it replaces the loop of form `for (T i = I; i <
   * N; ++i) { f(i); }`.
   * @param f  Callable object to be invoked in loop. It must accept argument of
   * type convertible from `T`, being the value of current loop index.
   *
   * Example:
   * \include static_loop.cc
   *
   * Result:
   * \verbinclude static_loop.out
   */
  static void body(auto&&) {}
};

template<std::integral T, T I, T N>
requires(I < N) struct static_loop<T, I, N>
{
  static void body(auto&& f)
  {
    f(std::integral_constant<T, I>{});
    static_loop<T, I + 1, N>::body(std::forward<decltype(f)>(f));
  }
};

//////////////////////////////////
// Functional logical operators //
//////////////////////////////////

/**
 * `fn_and` is an higher-order function returning conjunction of its arguments.
 *
 * @param fs Callable object returning Boolean value.
 *
 * Example:
 * \include fn_and.cc
 *
 * Result: (might be empty)
 * \verbinclude fn_and.out
 *
 * \note `fn_and` can be useful to describe complex genetic algorithm
 * termination conditions.
 */
auto
fn_and(const auto&... fs)
{
  return [=]<typename... Xs>(const Xs&... xs) { return (fs(xs...) && ...); };
}

/**
 * `fn_or` is an higher-order function returning disjunction of its arguments.
 *
 * @param fs Callable object returning Boolean value.
 *
 * Example:
 * \include fn_or.cc
 *
 * Result: (might be empty)
 * \verbinclude fn_or.out
 *
 * \note `fn_or` can be useful to describe complex genetic algorithm
 * termination conditions.
 */
auto
fn_or(const auto&... fs)
{
  return [=]<typename... Xs>(const Xs&... xs) { return (fs(xs...) || ...); };
}

/////////////////
// Thread pool //
/////////////////

/**
 * thread_pool implements modified thread pool design pattern.
 */
class thread_pool
{
public:
  /**
   * `thread_pool` constructor.
   *
   * @param sz Number of threads for concurrent calculations.
   */
  explicit thread_pool(std::size_t sz)
    : free_threads_{ sz }
  {}

  /**
   * `thread_pool::async` asynchronically executes callable object `f`
   * postponing start of `f` until number of concurrently executing threads in
   * pool drops below number `sz`, described in constructor.
   *
   * @param policy Lauch policy (see `std::launch` documentation).
   * @param f Callable object to be concurrently executed.
   *
   * Example:
   * \include thread_pool.cc
   *
   * Result: (might be different due to concurrent execution)
   * \verbinclude thread_pool.out
   */
  template<typename T>
  std::future<T> async(std::launch policy, const std::function<T()>& f)
  {
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
  inline void acquire()
  {
    std::unique_lock<std::mutex> ul{ m_ };
    cv_.wait(ul, [this]() { return free_threads_ != 0; });
    --free_threads_;
  }

  inline void release()
  {
    std::unique_lock<std::mutex> ul{ m_ };
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

/**
 * `callable` specifies that `F` is convertible to the callable object type,
 * which accepts arguments of types `Args...` and returns value of type `R`.
 *
 * \tparam F Type convertible to the callable object type.
 * \tparam R Return type.
 * \tparam Args Argument type.
 *
 * Example:
 * \include callable.cc
 *
 * Result: (might be empty)
 * \verbinclude callable.out
 */
template<typename F, typename R, typename... Args>
concept callable = std::convertible_to<std::invoke_result_t<F, Args...>, R>;

///////////
// Range //
///////////

/**
 * `range` represents subrange (closed interval) of type `T`.
 *
 * \tparam T Base type.
 */
template<typename T>
class range
{
public:
  /**
   * `range` constructor creates object representing closed interval \f$[{\rm
   * min}, {\rm max}]_{\rm T}\f$.
   *
   * @param min Range infimum.
   * @param max Range supremum.
   *
   * \throws std::invalid_argument Exception is raised if `min` is greater than
   * `max`.
   */
  constexpr range(T min, T max)
    : min_{ min }
    , max_{ max }
  {
    if (min > max) {
      throw std::invalid_argument{ "range: min greater than max" };
    }
  }

  /**
   * `range` constructor creates object representing full (closed) range for
   * type `T`.
   *
   * \note This constructor is available for types possesing
   * `std::numeric_limits` specialization.
   */
  template<typename U = T,
           typename = std::enable_if_t<std::numeric_limits<U>::is_specialized>>
  constexpr range()
    : range{ std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max() }
  {}

  constexpr range(const range&) = default;
  constexpr range(range&&) = default;
  range& operator=(const range&) = default;
  range& operator=(range&&) = default;

  /**
   * `range::min` returns range infimum.
   *
   * @return Range infimum, i.e. left endpoint of interval \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ (`min` value) corresponding to the range.
   */
  T min() const { return min_; }

  /**
   * `range::max` returns range supremum.
   *
   * @return Range supremum, i.e. right endpoint of interval \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ (`max` value) corresponding to the range.
   */
  T max() const { return max_; }

  /**
   * `range::midpoint` returns midpoint of interval represented by range.
   *
   * @return Range midpoint, i.e. center of interval \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ (arithmetic mean of `min` and `max` values) corresponding
   * to the range. If `T` is of integer type and the sum of `min` and `max` is
   * odd, the result is rounded towards `min`.
   *
   * \note This method is disabled for ranges of type `bool`.
   */
  template<typename U = T,
           typename = std::enable_if_t<!std::is_same_v<U, bool>>>
  T midpoint() const
  {
    return std::midpoint(min_, max_);
  }

  /**
   * `range::clamp` returns left endpoint if `t` compares less than left
   * endpoint; otherwise returns right endpoint if right endpoint compares less
   * than `t`; otherwise returns `t`.
   *
   * @param t Value to be clamped.
   * @return Value clamped to the range.
   */
  template<typename U = T,
           typename = std::enable_if_t<!std::is_same_v<U, bool>>>
  T clamp(T t) const
  {
    return std::clamp(t, min_, max_);
  }

  /**
   * `range::contains` checks if its argument is contained within interval
   * \f$[{\rm min}, {\rm max}]_{\rm T}\f$ represented by the range.
   *
   * @param t Argument to be checked.
   * @return Boolean value of check: `true` if `t` is in \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ interval and `false` otherwise.
   */
  bool contains(T t) const { return t >= min_ && t <= max_; }

  /**
   * `range::operator<=>` performs default lexicographical comparison with use
   * of left and right endpoints.
   *
   * @param r Range to be compared with `*this`.
   * @return Ordering (cf. `std::strong_ordering`, `std::weak_ordering`,
   * `std::partial_ordering`).
   */
  auto operator<=>(const range<T>& r) const = default;

private:
  T min_;
  T max_;
};

/**
 * `operator<<` prints range to the stream.
 *
 * @param os Stream to use.
 * @param r Range to be printed.
 * @return Reference to the `os` stream.
 */
template<typename T>
std::ostream&
operator<<(std::ostream& os, const range<T>& r)
{
  return (os << '[' << r.min() << ", " << r.max() << ']');
}

/**
 * `iota` returns `std::array` container filled with arithmetic sequence of
 * length `N`, with difference of 1 and starting from value `t`.
 *
 * \tparam T Number type.
 * \tparam N Returned container size.
 * @param t Starting value.
 * @return `std::array<T, N>` with consecutive numbers starting from value `t`.
 */
template<typename T, std::size_t N>
std::array<T, N>
iota(T t)
{
  std::array<T, N> res{};
  std::ranges::transform(std::views::iota(t) | std::views::take(N),
                         std::begin(res),
                         std::identity{});
  return res;
}

////////////////////
// Random numbers //
////////////////////

/**
 * `probability` type represents probability values, i.e. numbers from \f$[0,
 * 1]_{\mathbb{R}}\f$ interval.
 */
using probability = double;

/**
 * `random_engine` returns pseudo-random number generator engine based on
 * Mersenne Twister.
 * @return Reference to static object with Mersenne Twister engine
 * `std::mt19937` initialized with `std::random_device{}()`.
 */
inline std::mt19937&
random_engine()
{
  // Only one global hidden variable engine (regardless of number of
  // translation units).
  static std::mt19937 engine{ std::random_device{}() };
  return engine;
}

/**
 * `success` returns true with probability `success_probability` and false with
 * probability `1 - success_probability`, i.e. it implements Bernoulli
 * distribution \f${\rm B}(1, {\rm success\_{}probability})\f$.
 * @param success_probability Probability of returning `true` value.
 * @return Logic value drawn from \f${\rm B}(1, {\rm
 * success\_{}probability})\f$.
 */
inline bool
success(probability success_probability)
{
  return std::bernoulli_distribution{ success_probability }(random_engine());
}

/**
 * `random_N` returns random number from normal distribution with mean `mean`
 * and standard deviation `standard_deviation`.
 *
 * \tparam T Result type (floating-point).
 * @param mean Mean of normal distribution.
 * @param standard_deviation Standard deviation of normal distribution.
 * @return Number drawn from normal distribution.
 */
template<std::floating_point T>
T
random_N(T mean, T standard_deviation)
{
  auto& generator{ random_engine() };
  return std::normal_distribution<T>{ mean, standard_deviation }(generator);
}

/**
 * `random_U` returns random number from uniform distribution:
 *   - from interval \f$[a, b]_{\mathbb{R}}\f$ for floating-point type `T`
 *   - from set \f$\{a, b\}\f$ for Boolean type `T`
 *   - from interval \f$[a, b]_{\mathbb{Z}}\f$ for integer type `T`
 *
 * \tparam T Return type.
 * @param a Parameter describing aforementioned interval or set.
 * @param b Parameter describing aforementioned interval or set.
 * @return Value drawn from uniform distribution.
 *
 * \note For floating-point types overflow may occur for `std::nextafter(b,
 * std::numeric_limits<T>::max()) - a` (cf. N4861, 26.6.8.2.2).
 */
template<typename T>
T
random_U(T a, T b)
{
  auto& generator{ random_engine() };
  if constexpr (std::is_floating_point_v<T>) {
    assert(a < b);
    return std::uniform_real_distribution<T>{
      a, std::nextafter(b, std::numeric_limits<T>::max())
    }(generator); // [a, b]
  } else if constexpr (std::is_same_v<T, bool>) {
    return a == b ? a : std::bernoulli_distribution{ 0.5 }(generator);
  } else { // [a, b]
    assert(a <= b);
    return std::uniform_int_distribution<T>{ a, b }(generator);
  }
}

/////////////////////
// Some basic math //
/////////////////////

/**
 * `square` returns second power of its argument.
 *
 * \tparam T Argument and return type (floating-point or integer type).
 * \param x Argument to be raised to the second power.
 * \return Argument raised to the second power, i.e. \f$x^2\f$.
 */
template<typename T>
requires std::floating_point<T> || std::integral<T> T
square(T x)
{
  return x * x;
}

/**
 * `cube` returns thirdd power of its argument.
 *
 * \tparam T Argument and return type (floating-point or integer type).
 * \param x Argument to be raised to the third power.
 * \return Argument raised to the third power, i.e. \f$x^3\f$.
 */
template<typename T>
requires std::floating_point<T> || std::integral<T> T
cube(T x)
{
  return x * x * x;
}

template<std::floating_point T>
const T pi = std::numbers::pi_v<T>;

template<std::floating_point T>
const T e = std::numbers::e_v<T>;

template<std::floating_point T>
const T ln2 = std::numbers::ln2_v<T>;

namespace detail {

template<std::floating_point T>
T
angle(T x, T y)
{
  const T t = std::atan2(y, x);
  return t >= T{ 0 } ? t : (2 * std::numbers::pi_v<T> + t);
}

} // namespace detail

template<std::floating_point T>
std::tuple<T, T, T>
cart2spher(T x, T y, T z)
{
  return std::tuple<T, T, T>{
    std::hypot(x, y, z),                // r
    std::acos(z / std::hypot(x, y, z)), // theta
    detail::angle(x, y)                 // phi
  };
}

template<std::floating_point T>
std::tuple<T, T, T>
spher2cart(T r, T theta, T phi)
{
  return std::tuple<T, T, T>{ r * std::sin(theta) * std::cos(phi),
                              r * std::sin(theta) * std::sin(phi),
                              r * std::cos(theta) };
}

template<std::floating_point T>
std::tuple<T, T>
cart2polar(T x, T y)
{
  return std::tuple<T, T>{ std::hypot(x, y), detail::angle(x, y) };
}

template<std::floating_point T>
std::tuple<T, T>
polar2cart(T r, T phi)
{
  return std::tuple<T, T>{ r * std::cos(phi), r * std::sin(phi) };
}

////////////
// Domain //
////////////

template<typename T, std::size_t N>
using domain = std::array<range<T>, N>;

template<typename T>
struct is_domain : std::false_type
{};

template<typename T, std::size_t N>
struct is_domain<domain<T, N>> : std::true_type
{};

template<typename T>
inline constexpr bool is_domain_v = is_domain<T>::value;

template<typename T>
concept set_of_departure = is_domain_v<T>;

template<typename T, std::size_t N>
bool
contains(const domain<T, N>& d, const std::array<T, N>& p)
{
  bool res = true;
  for (std::size_t i = 0; res && i < N; ++i) {
    res &= d[i].contains(p[i]);
  }
  return res;
}

template<typename T, std::size_t N>
constexpr domain<T, N>
uniform_domain(const range<T>& r)
{
  domain<T, N> res{};
  std::generate_n(std::begin(res), N, [&]() { return r; });
  return res;
}

template<typename T, std::size_t N>
constexpr domain<T, N>
uniform_domain(T lo, T hi)
{
  return uniform_domain<T, N>(range<T>{ lo, hi });
}

template<typename T, std::size_t N>
constexpr bool
uniform(const domain<T, N>& d)
{
  const auto x0 = N == 0 ? range<T>{} : d[0];
  return std::ranges::all_of(d, [&](const auto& x) { return x0 == x; });
}

template<typename T, std::size_t N>
constexpr domain<T, 2 * N>
self_adaptive_variation_domain(const domain<T, N>& d, T lo)
{
  const T s = .5;
  domain<T, 2 * N> res{};
  for (std::size_t i = 0; i < N; ++i) {
    res[i] = d[i];
    res[i + N] =
      range{ lo, s * std::max(std::fabs(d[i].min()), std::fabs(d[i].max())) };
  }
  return res;
}

template<typename T, std::size_t N>
using chain = std::array<T, N>;

template<typename T, std::size_t N>
chain<T, N>
chain_min(const domain<T, N>& d)
{
  chain<T, N> res{};
  std::ranges::transform(d, std::begin(res), std::identity{}, &range<T>::min);
  return res;
}

//////////////
// Genotype //
//////////////

template<typename T, std::size_t N, const domain<T, N>* D>
requires std::floating_point<T>
struct g_floating_point
{
  static_assert(D != nullptr);
  static_assert(N > 0);
  using type = T;
  static constexpr std::size_t size() { return N; }
  static constexpr const domain<type, size()>& constraints() { return *D; }
  using chain_t = chain<type, size()>;

  static bool valid(const chain<type, size()>& c)
  {
    return contains(constraints(), c);
  }

  static chain_t default_chain() { return chain_min(constraints()); }
};

template<typename T>
struct is_g_floating_point : std::false_type
{};

template<typename T, std::size_t N, const domain<T, N>* D>
struct is_g_floating_point<g_floating_point<T, N, D>> : std::true_type
{};

template<typename T>
inline constexpr bool is_g_floating_point_v = is_g_floating_point<T>::value;

template<typename T>
concept floating_point_representation = is_g_floating_point_v<T>;

template<typename T, std::size_t N, const domain<T, N>* D>
requires std::integral<T> &&(!std::is_same_v<T, bool>)struct g_integer
{
  static_assert(D != nullptr);
  static_assert(N > 0);
  using type = T;
  static constexpr std::size_t size() { return N; }
  static constexpr const domain<type, size()>& constraints() { return *D; }
  using chain_t = chain<type, size()>;

  static bool valid(const chain<type, size()>& c)
  {
    return contains(constraints(), c);
  }

  static chain_t default_chain() { return chain_min(constraints()); }
};

template<typename T>
struct is_g_integer : std::false_type
{};

template<typename T, std::size_t N, const domain<T, N>* D>
struct is_g_integer<g_integer<T, N, D>> : std::true_type
{};

template<typename T>
inline constexpr bool is_g_integer_v = is_g_integer<T>::value;

template<typename T>
concept integer_representation = is_g_integer_v<T>;

template<std::size_t N>
struct g_binary
{
  static_assert(N > 0);
  using type = bool;
  static constexpr std::size_t size() { return N; }
  static constexpr const domain<type, size()> constraints()
  {
    return domain<type, size()>{};
  }

  using chain_t = chain<type, size()>;

  static bool valid(const chain<type, size()>&) { return true; }
  static chain_t default_chain() { return chain_min(constraints()); }
};

template<typename T>
struct is_g_binary : std::false_type
{};

template<std::size_t N>
struct is_g_binary<g_binary<N>> : std::true_type
{};

template<typename T>
inline constexpr bool is_g_binary_v = is_g_binary<T>::value;

template<typename T>
concept binary_representation = is_g_binary_v<T>;

template<typename T, std::size_t N, T M>
requires std::integral<T> &&(!std::is_same_v<T, bool>)struct g_permutation
{
  using type = T;
  static constexpr std::size_t size() { return N; }

  static constexpr const domain<type, size()> constraints()
  {
    return uniform_domain<type, size()>(M, M + N - 1);
  }

  using chain_t = chain<type, size()>;

  static bool valid(const chain<type, size()>& c)
  {
    const auto i = iota<type, size()>(M);
    return contains(constraints(), c) &&
           std::is_permutation(std::begin(c), std::end(c), std::begin(i));
  }

  static chain_t default_chain() { return iota<type, size()>(M); }
};

template<typename T>
struct is_g_permutation : std::false_type
{};

template<typename T, std::size_t N, T M>
struct is_g_permutation<g_permutation<T, N, M>> : std::true_type
{};

template<typename T>
inline constexpr bool is_g_permutation_v = is_g_permutation<T>::value;

template<typename T>
concept permutation_representation = is_g_permutation_v<T>;

template<typename T>
concept chromosome_representation =
  floating_point_representation<T> || integer_representation<T> ||
  binary_representation<T> || permutation_representation<T>;

template<typename R>
requires chromosome_representation<R>
class genotype
{
public:
  using chain_t = chain<typename R::type, R::size()>;
  using const_iterator = typename chain_t::const_iterator;
  using gene_t = typename R::type;
  using genotype_t = R;
  static constexpr std::size_t size() { return R::size(); }

  static constexpr const domain<gene_t, size()> constraints()
  {
    return R::constraints();
  }

  static constexpr bool uniform_domain = uniform(constraints());
  static bool valid(const chain_t& c) { return R::valid(c); }

public:
  genotype()
    : chain_{ R::default_chain() }
  {}

  explicit genotype(const chain_t& c)
    : chain_{ c }
  {
    if (!valid(c)) {
      throw std::invalid_argument{ "invalid chain" };
    }
  }

  genotype(const genotype&) = default;
  genotype(genotype&&) = default;
  genotype& operator=(const genotype&) = default;
  genotype& operator=(genotype&&) = default;

  gene_t value(std::size_t i) const { return chain_[i]; }

  template<typename = std::enable_if_t<!permutation_representation<R>>>
  genotype& value(std::size_t i, gene_t v)
  {
    if (!constraints()[i].contains(v)) {
      throw std::invalid_argument{ "bad value" };
    }
    chain_[i] = v;
    return *this;
  }

  genotype& random_reset()
  {
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
  genotype& random_reset(std::size_t i)
  {
    const auto& c = constraints();
    chain_[i] = random_U<gene_t>(c[i].min(), c[i].max());
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
struct is_genotype : std::false_type
{};

template<typename T>
struct is_genotype<genotype<T>> : std::true_type
{};

template<typename T>
inline constexpr bool is_genotype_v = is_genotype<T>::value;

template<typename G>
concept chromosome = is_genotype_v<G>;

template<typename G>
concept floating_point_chromosome =
  chromosome<G> && floating_point_representation<typename G::genotype_t>;

template<typename G>
concept integer_chromosome =
  chromosome<G> && integer_representation<typename G::genotype_t>;

template<typename G>
concept binary_chromosome =
  chromosome<G> && binary_representation<typename G::genotype_t>;

template<typename G>
concept permutation_chromosome =
  chromosome<G> && permutation_representation<typename G::genotype_t>;

template<typename G>
concept uniform_chromosome = chromosome<G> && G::uniform_domain;

template<typename F, typename G>
concept genotype_constraints = std::predicate<F, G> && chromosome<G>;

template<typename G>
requires chromosome<G>
const auto constraints_satisfied = [](const G&) { return true; };

template<typename G>
requires chromosome<G> std::ostream&
operator<<(std::ostream& os, const G& g)
{
  for (std::size_t i = 0; i < G::size(); ++i) {
    os << g.value(i) << (i + 1 < G::size() ? " " : "");
  }
  return os;
}

template<typename G>
requires floating_point_chromosome<G> std::ostream&
operator<<(std::ostream& os, const G& g)
{
  const auto prec = std::numeric_limits<typename G::gene_t>::digits10;
  for (std::size_t i = 0; i < G::size(); ++i) {
    os << std::scientific << std::setprecision(prec) << g.value(i)
       << (i + 1 < G::size() ? " " : "");
  }
  return os;
}

} // namespace quile

template<typename G>
requires quile::chromosome<G>
struct std::hash<G>
{
  std::size_t operator()(const G& g) const noexcept
  {
    const std::size_t sz{ sizeof(std::size_t) * CHAR_BIT };
    std::size_t res{ 0 };
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

template<typename G>
requires chromosome<G>
using population = std::vector<G>;

template<typename T>
struct is_population : std::false_type
{};

template<typename G>
struct is_population<population<G>> : std::true_type
{};

template<typename G>
inline constexpr bool is_population_v = is_population<G>::value;

template<typename G>
concept genetic_pool = is_population_v<G>;

// Population generators/selectors
// - first generation creator
template<typename G>
requires chromosome<G>
using populate_0_fn = std::function<population<G>(std::size_t)>;
// - parents selection
template<typename G>
requires chromosome<G>
using populate_1_fn =
  std::function<population<G>(std::size_t, const population<G>&)>;
// - survivor selection
template<typename G>
requires chromosome<G>
using populate_2_fn = std::function<
  population<G>(std::size_t, const population<G>&, const population<G>&)>;

template<typename G>
requires chromosome<G>
using generations = std::deque<population<G>>;

//////////////////////////////
// Mutation & recombination //
//////////////////////////////

template<typename M, typename G>
concept mutation = requires(M m, G g)
{
  {
    m(g)
    } -> std::convertible_to<population<G>>;
}
&&chromosome<G>;

template<typename G>
requires chromosome<G>
using mutation_fn = std::function<population<G>(const G&)>;

template<typename R, typename G>
concept recombination = requires(R r, G g)
{
  {
    r(g, g)
    } -> std::convertible_to<population<G>>;
}
&&chromosome<G>;

template<typename G>
requires chromosome<G>
using recombination_fn = std::function<population<G>(const G&, const G&)>;

template<typename G>
requires chromosome<G> population<G>
unary_identity(const G& g)
{
  return population<G>{ g };
}

template<typename G>
requires chromosome<G> population<G>
binary_identity(const G& g0, const G& g1)
{
  return population<G>{ g0, g1 };
}

template<typename G>
requires chromosome<G>
class variation
{
public:
  variation(const mutation_fn<G>& m, const recombination_fn<G>& r)
    : m_{ m }
    , r_{ r }
  {}

  variation()
    : variation{ unary_identity<G>, binary_identity<G> }
  {}

  explicit variation(const mutation_fn<G>& m)
    : variation{ m, binary_identity<G> }
  {}

  explicit variation(const recombination_fn<G>& r)
    : variation{ unary_identity<G>, r }
  {}

  population<G> operator()(const G& g0, const G& g1) const
  {
    QUILE_LOG("Variation: " << g0 << ", " << g1);
    population<G> res{};
    for (const auto& g : r_(g0, g1)) {
      res.push_back(m_(g).at(0));
    }
    assert(res.size() == 1 || res.size() == 2);
    return res;
  }

  population<G> operator()(const population<G>& p) const
  {
    if (p.size() % 2) {
      throw std::invalid_argument{ "wrong population size" };
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

template<typename G>
requires chromosome<G>
auto
stochastic_mutation(const mutation_fn<G>& m, probability p)
{
  return [=](const G& g) { return success(p) ? m(g) : population<G>{ g }; };
}

template<typename G>
requires chromosome<G>
auto
stochastic_recombination(const recombination_fn<G>& r, probability p)
{
  return [=](const G& g0, const G& g1) {
    const auto tmp = r(g0, g1);
    return success(p)        ? tmp
           : tmp.size() == 2 ? population<G>{ g0, g1 }
           : success(.5)     ? population<G>{ g0 }
                             : population<G>{ g1 };
  };
}

///////////////
// Evolution //
///////////////

template<typename F, typename G>
concept termination_condition = std::predicate<F, std::size_t, generations<G>>;

template<typename G>
requires chromosome<G>
using termination_condition_fn =
  std::function<bool(std::size_t, const generations<G>&)>;

// TODO: Is there any way to reduce number of arguments of this function
// without increasing solution's complexity?
template<typename G>
requires chromosome<G> generations<G>
evolution(const variation<G> v,
          const population<G>& first_generation,
          const populate_1_fn<G>& p1,
          const populate_2_fn<G>& p2,
          const termination_condition_fn<G>& tc,
          std::size_t parents_sz,
          std::size_t max_history = 0)
{
  generations<G> res{};
  const std::size_t generation_sz = first_generation.size();
  for (std::size_t i = 0; !tc(i, res); ++i) {
    QUILE_LOG("Generation #" << i);
    const population<G> p{
      i == 0 ? first_generation
             : p2(generation_sz, res.back(), v(p1(parents_sz, res.back())))
    };
    res.push_back(p);
    if (max_history && res.size() > max_history) {
      res.pop_front();
    }
  }
  return res;
}

template<typename G>
requires chromosome<G> generations<G>
evolution(const variation<G>& v,
          const populate_0_fn<G>& p0,
          const populate_1_fn<G>& p1,
          const populate_2_fn<G>& p2,
          const termination_condition_fn<G>& tc,
          std::size_t generation_sz,
          std::size_t parents_sz,
          std::size_t max_history = 0)
{
  return evolution<G>(
    v, p0(generation_sz), p1, p2, tc, parents_sz, max_history);
}

//////////////////////
// Fitness function //
//////////////////////

using fitness = double;
using fitnesses = std::vector<fitness>;

template<typename G>
requires chromosome<G>
using fitness_function = std::function<fitness(const G&)>;

const fitness incalculable = -std::numeric_limits<fitness>::infinity();

template<typename G>
requires chromosome<G>
class fitness_db
{
private:
  using database = std::unordered_map<G, fitness>;

public:
  using const_iterator = typename database::const_iterator;

public:
  explicit fitness_db(
    const fitness_function<G>& f,
    const genotype_constraints<G> auto& gc,
    unsigned int thread_sz = std::thread::hardware_concurrency())
    : function_{ [=](const G& g) { return gc(g) ? f(g) : incalculable; } }
    , thread_sz_{ thread_sz }
  {}

  fitness_db(const fitness_db&) = default;
  fitness_db& operator=(const fitness_db&) = default;

  fitness operator()(const G& g) const
  {
    const auto it{ fitness_values_->find(g) };
    const bool b = it != fitness_values_->end();
    const fitness res = b ? it->second : ((*fitness_values_)[g] = function_(g));
    QUILE_LOG("Fitness value for ["
              << g << "]: " << res
              << (b ? " (taken from database)" : " (calculated on demand)"));
    return res;
  }

  fitnesses operator()(const population<G>& p) const
  {
    if (thread_sz_ > 1 && p.size() > 1) {
      multithreaded_calculations(p);
    }
    fitnesses res{};
    QUILE_LOG("Fitness values for population of size " << p.size());
    std::ranges::transform(
      p, std::back_inserter(res), [this](const G& g) { return operator()(g); });
    return res;
  }

  std::size_t size() const { return fitness_values_->size(); }

  const_iterator begin() const { return fitness_values_->begin(); }
  const_iterator end() const { return fitness_values_->end(); }

  population<G> rank_order() const
  {
    population<G> res{};
    std::ranges::transform(*this,
                           std::back_inserter(res),
                           std::identity{},
                           &database::value_type::first);
    std::ranges::sort(res, [this](const G& g0, const G& g1) {
      return this->operator()(g0) > this->operator()(g1);
    });
    return res;
  }

private:
  auto uncalculated_fitnesses(const population<G>& p) const
  {
    std::unordered_set<G> res{};
    std::ranges::copy_if(
      p, std::inserter(res, std::end(res)), [this](const G& g) {
        return !fitness_values_->contains(g);
      });
    return res;
  }

  void multithreaded_calculations(const population<G>& p) const
  {
    using type = std::pair<G, fitness>;
    thread_pool tp{ thread_sz_ };
    std::vector<std::future<type>> v{};
    for (const auto& x : uncalculated_fitnesses(p)) {
      QUILE_LOG("Asynchronous fitness value calculations (multithreaded)");
      v.push_back(tp.async<type>(std::launch::async, [this, x]() {
        const fitness xf = this->function_(x);
        return type{ x, xf };
      }));
    }
    for (auto& x : v) {
      const auto p = x.get();
      QUILE_LOG("Fitness value for ["
                << p.first << "]: " << p.second
                << " (calculated asynchronously on demand)");
      fitness_values_->insert(p);
    }
  }

private:
  fitness_function<G> function_;
  unsigned int thread_sz_;
  std::shared_ptr<database> fitness_values_ = std::make_shared<database>();
};

template<typename G>
requires chromosome<G>
void
print(std::ostream& os,
      const generations<G>& gs,
      const fitness_db<G>* fd = nullptr)
{
  for (std ::size_t i = 0; const auto& x : gs) {
    const auto prec = std::numeric_limits<fitness>::digits10;
    for (const auto& xx : x) {
      if (fd) {
        os << i << ' ' << xx << ' ' << std::scientific
           << std::setprecision(prec) << (*fd)(xx) << '\n';
      } else {
        os << i << ' ' << xx << '\n';
      }
    }
    ++i;
  }
}

/////////////////////////////
// Selection probabilities //
/////////////////////////////

using selection_probabilities = std::vector<probability>;

template<typename G>
requires chromosome<G>
using selection_probabilities_fn =
  std::function<selection_probabilities(const population<G>&)>;

template<typename G>
requires chromosome<G> selection_probabilities
cumulative_probabilities(const selection_probabilities_fn<G>& spf,
                         const population<G>& p)
{
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
C
select_different_than(const C& c, T t, bool require_nonempty_result)
{
  C res{};
  std::ranges::copy_if(
    c, std::back_inserter(res), [=](auto x) { return x != t; });
  if (require_nonempty_result && res.size() == 0) {
    throw std::runtime_error{ "empty result" };
  }
  return res;
}

inline fitnesses
select_calculable(const fitnesses& fs, bool require_nonempty_result = false)
{
  return select_different_than(fs, incalculable, require_nonempty_result);
}

template<typename G>
class fitness_proportional_selection
{
public:
  explicit fitness_proportional_selection(const fitness_db<G>& ff)
    : ff_{ ff }
  {}

  // FPS with windowing with workarounds for:
  // a) population of equally fit genotypes and
  // b) populations containing genotypes which fitnesses cannot be calculated
  // Please note that in b) case, there should be at least one genotype,
  // which fitness can be calculated.
  selection_probabilities operator()(const population<G>& p) const
  {
    const fitnesses fs{ ff_(p) };
    const auto cal = select_calculable(fs, true);
    const fitness min = *std::ranges::min_element(cal);
    const auto n = cal.size(); // Value of n is guaranteed to be greater than 0.
    const fitness delta = 1. / n;
    const fitness sum =
      std::accumulate(std::begin(cal), std::end(cal), fitness{ 0. }) - n * min +
      1;
    selection_probabilities res{};
    std::ranges::transform(fs, std::back_inserter(res), [=](fitness f) {
      return f == incalculable ? .0 : (f - min + delta) / sum;
    });
    return res;
  }

private:
  const fitness_db<G> ff_;
};

namespace detail {

template<typename It>
It
advance_cpy(It it, std::size_t n)
{
  std::advance(it, n);
  return it;
}

template<typename It, typename Compare = std::less<>>
std::vector<std::size_t>
rank(It first, It last, Compare comp = {})
{
  std::vector<std::size_t> v(std::distance(first, last), 0);
  std::ranges::generate(v, [i = std::size_t{ 0 }]() mutable { return i++; });
  std::ranges::stable_sort(v, [first, &comp](auto a, auto b) {
    return comp(*advance_cpy(first, a), *advance_cpy(first, b));
  });
  std::vector<std::size_t> res(v.size(), 0);
  for (std::size_t i = 0; auto x : v) {
    res.at(x) = i++;
  }
  return res;
}

template<typename T, typename U>
T
id(U u)
{
  return u;
}

} // namespace detail

inline auto
linear_ranking_selection(double s)
{
  return [s](std::size_t mu, std::size_t j) -> probability {
    assert(mu > 0 && 1. < s && s <= 2. && j < mu);
    const auto id = detail::id<probability, std::size_t>;
    return mu == 1
             ? 1.
             : (2 - s) / id(mu) + 2 * id(j) * (s - 1) / (id(mu) * (id(mu) - 1));
  };
}

inline probability
exponential_ranking_selection(std::size_t mu, std::size_t j)
{
  assert(mu > 0 && j < mu);
  const auto e = std::numbers::e_v<probability>;
  const auto id = detail::id<probability, std::size_t>;
  // Consider difference between -id(mu) and -mu for std::size_t mu.
  return mu == 1 ? 1.
                 : (1. - std::exp(-id(j))) * (1. - e) /
                     (id(mu) * (1. - e) + e - std::exp(1 - id(mu)));
}

template<typename G>
class ranking_selection
{
private:
  using probability_fn = std::function<probability(std::size_t, std::size_t)>;

public:
  ranking_selection(const fitness_db<G>& ff, const probability_fn& pf)
    : ff_{ ff }
    , pf_{ pf }
  {}

  // RS with workarounds for populations containing genotypes which fitnesses
  // cannot be calculated. Please note that there should be at least one
  // genotype, which fitness can be calculated.
  selection_probabilities operator()(const population<G>& p) const
  {
    const fitnesses fs{ ff_(p) };
    auto r = detail::rank(
      std::begin(p), std::end(p), [&](const auto& a, const auto& b) {
        return ff_(a) < ff_(b);
      });
    const auto mu = select_calculable(fs, true).size();
    const auto nq = p.size() - mu;
    selection_probabilities res{};
    std::ranges::transform(r, std::back_inserter(res), [=, this](auto j) {
      return j < nq ? 0. : pf_(mu, j - nq);
    });
    return res;
  }

private:
  const fitness_db<G> ff_;
  const probability_fn pf_;
};

////////////////////////////
// Generators & selectors //
////////////////////////////

namespace detail {

template<typename G>
requires chromosome<G> population<G>
generate(std::size_t lambda, const std::function<G()>& f)
{
  population<G> res{};
  std::generate_n(std::back_inserter(res), lambda, f);
  return res;
}

} // namespace detail

template<typename G>
requires chromosome<G> populate_2_fn<G>
adapter(const populate_1_fn<G>& fn)
{
  return [=](std::size_t sz, const population<G>& p0, const population<G>& p1) {
    population<G> p{ p0 };
    p.insert(p.end(), p1.begin(), p1.end());
    return fn(sz, p);
  };
}

template<auto C, typename G>
requires genotype_constraints<decltype(C), G> && chromosome<G> population<G>
random_population(std::size_t lambda)
{
  population<G> res{};
  auto g = G{};
  for (std::size_t i = 0; i < lambda; ++i) {
    while (!C(g.random_reset()))
      ;
    res.push_back(g);
  }
  return res;
}

template<typename G>
requires chromosome<G>
class roulette_wheel_selection
{
public:
  explicit roulette_wheel_selection(const selection_probabilities_fn<G>& spf)
    : spf_{ spf }
  {}

  population<G> operator()(std::size_t lambda, const population<G>& p) const
  {
    QUILE_LOG("Roulette wheel selection");
    const auto f = [&, c = cumulative_probabilities(spf_, p)]() -> G {
      return p.at(std::distance(
        c.begin(),
        std::lower_bound(c.begin(), c.end(), random_U<double>(0., 1.))));
    };
    return detail::generate<G>(lambda, f);
  }

private:
  const selection_probabilities_fn<G> spf_;
};

template<typename G>
requires chromosome<G>
class stochastic_universal_sampling
{
public:
  explicit stochastic_universal_sampling(
    const selection_probabilities_fn<G>& spf)
    : spf_{ spf }
  {}

  population<G> operator()(std::size_t lambda, const population<G>& p) const
  {
    QUILE_LOG("Stochastic Universal Sampling");
    const auto a = cumulative_probabilities(spf_, p);
    auto r = random_U<double>(0., 1. / lambda);

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

template<typename G>
requires chromosome<G> population<G>
generational_survivor_selection(std::size_t sz,
                                const population<G>& generation,
                                const population<G>& offspring)
{
  if (generation.size() != sz || offspring.size() != sz) {
    throw std::invalid_argument{ "bad size" };
  }
  return offspring;
}

///////////////////////////
// Termination condition //
///////////////////////////

inline fitness
max(const fitnesses& fs)
{
  const fitnesses calc{ select_calculable(fs, true) };
  return *std::ranges::max_element(calc);
}

template<typename G>
requires chromosome<G> fitness
max(const population<G>& p, const fitness_db<G>& ff)
{
  return max(ff(p));
}

template<typename G>
requires chromosome<G> fitnesses
max(const generations<G>& gs, const fitness_db<G>& ff)
{
  fitnesses res{};
  std::ranges::transform(gs,
                         std::back_inserter(res),
                         [&ff](const population<G>& p) { return max(p, ff); });
  return res;
}

inline fitness
min(const fitnesses& fs)
{
  const fitnesses calc{ select_calculable(fs, true) };
  return *std::ranges::min_element(calc);
}

template<typename G>
requires chromosome<G> fitness
min(const population<G>& p, const fitness_db<G>& ff)
{
  return min(ff(p));
}

template<typename G>
requires chromosome<G> fitnesses
min(const generations<G>& gs, const fitness_db<G>& ff)
{
  fitnesses res{};
  std::ranges::transform(gs,
                         std::back_inserter(res),
                         [&ff](const population<G>& p) { return min(p, ff); });
  return res;
}

template<typename G>
termination_condition_fn<G>
max_iterations_termination(std::size_t max)
{
  return [=](std::size_t i, const generations<G>&) { return i == max; };
}

template<typename G>
termination_condition_fn<G>
max_fitness_improvement_termination(const fitness_db<G>& ff,
                                    std::size_t n,
                                    double frac)
{
  return [=]([[maybe_unused]] std::size_t i, const generations<G>& gs) {
    assert(i == gs.size());
    if (gs.size() <= n) {
      return false;
    } else {
      const fitnesses fs{ max(gs, ff) };
      const fitness min_last_n = *std::min_element(fs.end() - n, fs.end());
      const double x = (max(fs) - min_last_n) / (max(fs) - min(fs));
      return x <= frac;
    }
  };
}

template<typename G>
termination_condition_fn<G>
max_fitness_improvement_termination_2(const fitness_db<G>& ff,
                                      std::size_t n,
                                      fitness delta)
{
  return [=]([[maybe_unused]] std::size_t i, const generations<G>& gs) {
    assert(i == gs.size() && delta >= .0);
    if (gs.size() <= n) {
      return false;
    } else {
      const fitnesses fs{ max(gs, ff) };
      const fitness max_0 = *std::max_element(fs.begin(), fs.end() - n);
      const fitness max_1 = *std::max_element(fs.end() - n, fs.end());
      return max_1 <= max_0 + delta;
    }
  };
}

template<typename G, typename F>
requires chromosome<G> && std::predicate<F, G> termination_condition_fn<G>
threshold_termination(const F& thr)
{
  return [=](std::size_t, const generations<G>& gs) {
    return gs.empty() ? false : std::ranges::any_of(gs.back(), thr);
  };
}

template<typename G>
termination_condition_fn<G>
fitness_threshold_termination(const fitness_db<G>& fd, fitness thr, fitness eps)
{
  const auto f = [=](const G& g) { return std::fabs(fd(g) - thr) <= eps; };
  return threshold_termination<G, decltype(f)>(f);
}

/////////////////////////////////////////////////
// Concrete mutation & recombination operators //
/////////////////////////////////////////////////

template<typename G>
requires floating_point_chromosome<G>
auto
Gaussian_mutation(typename G::gene_t sigma, probability p)
{
  return [=](const G& g) -> population<G> {
    G res{};
    const auto c = G::constraints();
    for (std::size_t i = 0; i < G::size(); ++i) {
      if (success(p)) {
        res.value(i, c[i].clamp(g.value(i) + sigma * random_N(0., 1.)));
      }
    }
    return population<G>{ res };
  };
}

template<typename G>
requires floating_point_chromosome<G> &&
  (G::size() % 2 == 0) auto self_adaptive_mutation(typename G::gene_t a0,
                                                   typename G::gene_t a1)
{
  return [=, n = G::size() / 2, c = G::constraints()](const G& g) {
    using type = typename G::gene_t;
    const type p0 = random_N(0., 1.) * a0 / std::sqrt(2 * n);
    const type t1 = a1 / std::sqrt(2 * std::sqrt(n));
    G res{};
    for (std::size_t i = 0; i < n; ++i) {
      const type sigma =
        c[i + n].clamp(g.value(i + n) * std::exp(p0 + t1 * random_N(0., 1.)));
      res.value(i, c[i].clamp(g.value(i) + sigma * random_N(0., 1.)));
      res.value(i + n, sigma);
    }
    return population<G>{ g };
  };
}

template<typename G>
requires uniform_chromosome<G> population<G>
swap_mutation(const G& g)
{
  const std::size_t n = G::size();
  auto d = g.data();
  std::swap(d[random_U<std::size_t>(0, n - 1)],
            d[random_U<std::size_t>(0, n - 1)]);
  return population<G>{ G{ d } };
}

template<typename G>
requires floating_point_chromosome<G> || integer_chromosome<G> ||
  binary_chromosome<G>
auto
random_reset(probability p)
{
  return [=](const G& g) -> population<G> {
    G res{ g };
    for (std::size_t i = 0; i < G::size(); ++i) {
      if (success(p)) {
        res.random_reset(i);
      }
    }
    return population<G>{ res };
  };
}

template<typename G>
requires binary_chromosome<G>
auto
bit_flipping(probability p)
{
  return [=](const G& g) -> population<G> {
    G res{ g };
    for (std::size_t i = 0; i < G::size(); ++i) {
      if (success(p)) {
        res.value(i, !res.value(i));
      }
    }
    return population<G>{ res };
  };
}

template<typename G>
requires floating_point_chromosome<G> population<G>
arithmetic_recombination(const G& g0, const G& g1)
{
  G res{};
  for (std::size_t i = 0; i < G::size(); ++i) {
    res.value(i, std::midpoint(g0.value(i), g1.value(i)));
  }
  return population<G>{ res };
}

template<typename G>
requires floating_point_chromosome<G> population<G>
single_arithmetic_recombination(const G& g0, const G& g1)
{
  G res0{ g0 };
  G res1{ g1 };
  const auto cp = random_U<std::size_t>(0, G::size() - 1);
  const auto mid = std::midpoint(res0.value(cp), res1.value(cp));
  res0.value(cp, mid);
  res1.value(cp, mid);
  return population<G>{ res0, res1 };
}

template<typename G>
requires floating_point_chromosome<G> || integer_chromosome<G> ||
  binary_chromosome<G>
    population<G>
    one_point_xover(const G& g0, const G& g1)
{
  auto d0 = g0.data();
  auto d1 = g1.data();
  const std::size_t n = G::size();
  const auto cp = random_U<std::size_t>(0, n - 1);
  for (std::size_t i = cp; i < n; ++i) {
    std::swap(d0[i], d1[i]);
  }
  return population<G>{ G{ d0 }, G{ d1 } };
}

template<typename G>
requires permutation_chromosome<G> population<G>
cut_n_crossfill(const G& g0, const G& g1)
{
  const auto f = [cp = random_U<std::size_t>(1, G::size() - 1)](const G& g,
                                                                auto d) {
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
  return population<G>{ G{ f(g1, g0.data()) }, G{ f(g0, g1.data()) } };
}

////////////////////////////////////////////////////
// Test functions for floating-point optimization //
////////////////////////////////////////////////////

namespace test_functions
{

  template<std::floating_point T, std::size_t N>
  using point = std::array<T, N>;

  template<std::floating_point T, std::size_t N>
  T distance(const point<T, N>& p0, const point<T, N>& p1)
  {
    T res = .0;
    for (std::size_t i = 0; i < N; ++i) {
      res += square(p0[i] - p1[i]);
    }
    return std::sqrt(res);
  }

  template<std::floating_point T>
  std::tuple<T, T> coordinates(const point<T, 2>& p)
  {
    return std::tuple<T, T>{ p[0], p[1] };
  }

  template<std::floating_point T>
  std::tuple<T, T, T> coordinates(const point<T, 3>& p)
  {
    return std::tuple<T, T, T>{ p[0], p[1], p[2] };
  }

  template<std::floating_point T, std::size_t N>
  point<T, N> uniform_point(T v)
  {
    point<T, N> res{};
    std::ranges::generate(res, [=]() -> T { return v; });
    return res;
  }

  template<std::floating_point T, std::size_t N>
  class test_function
  {
  public:
    using function = std::function<T(const point<T, N>&)>;
    using domain_fn = std::function<domain<T, N>()>;
    using point_fn = std::function<point<T, N>()>;

  public:
    test_function(const std::string& name,
                  const function& fn,
                  const domain_fn& d,
                  const point_fn& p_min)
      : name_{ name }
      , fn_{ fn }
      , d_{ d }
      , p_min_{ p_min }
    {}

    std::string name() const { return name_; }
    T operator()(const point<T, N>& p) const { return fn_(p); }
    domain<T, N> function_domain() const { return d_(); }
    point<T, N> p_min() const { return p_min_(); }

  private:
    std::string name_;
    function fn_;
    domain_fn d_;
    point_fn p_min_;
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> Ackley{
    "Ackley",
    [](const point<T, N>& p) {
      T s0 = 0.;
      T s1 = 0.;
      for (auto x : p) {
        s0 += square(x);
        s1 += std::cos(2 * pi<T> * x);
      }
      return -20. * std::exp(-.02 * std::sqrt(s0) / std::sqrt(N)) -
             std::exp(s1 / N) + 20. + e<T>;
    },
    []() { return uniform_domain<T, N>(-35., 35.); },
    []() { return uniform_point<T, N>(0.); }
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> Alpine{
    "Alpine",
    [](const point<T, N>& p) {
      return std::transform_reduce(
        std::begin(p), std::end(p), T{ 0. }, std::plus<T>{}, [](auto x) {
          return std::fabs(x * std::sin(x) + .1 * x);
        });
    },
    []() { return uniform_domain<T, N>(-10., 10.); },
    []() { return uniform_point<T, N>(0.); }
  };

  template<std::floating_point T>
  const test_function<T, 2> Aluffi_Pentini{
    "Aluffi-Pentini",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return ((.25 * x * x - .5) * x + .1) * x + 0.5 * y * y;
    },
    []() { return uniform_domain<T, 2>(-10., 10.); },
    []() {
      return point<T, 2>{
        [q = 0.1](int k) -> T {
          return 2. * std::sqrt(3.) *
                 std::cos(std::acos(-3. * std::sqrt(3.) * q / 2.) / 3. -
                          2. * std::numbers::pi_v<T> * k / 3.) /
                 3.;
        }(2),
        0.
      };
    }
  };

  template<std::floating_point T>
  const test_function<T, 2> Booth{
    "Booth",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return square(x + 2. * y - 7.) + square(2. * x + y - 5.);
    },
    []() { return uniform_domain<T, 2>(-10., 10.); },
    []() {
      return point<T, 2>{ 1., 3. };
    }
  };

  template<std::floating_point T>
  const test_function<T, 4> Colville{
    "Colville",
    [](const point<T, 4>& p) {
      return 100. * square(p[0] - square(p[1])) + square(1. - p[0]) +
             90. * square(p[3] - p[2] * p[2]) + square(1. - p[2]) +
             10.1 * square(p[1] - 1.) + square(p[3] - 1.) +
             19.8 * (p[1] - 1.) * (p[3] - 1.);
    },
    []() { return uniform_domain<T, 4>(-10., 10.); },
    []() { return uniform_point<T, 4>(1.); }
  };

  template<std::floating_point T>
  const test_function<T, 2> Easom{
    "Easom",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return -std::cos(x) * std::cos(y) *
             std::exp(-square(x - pi<T>) - square(y - pi<T>));
    },
    []() { return uniform_domain<T, 2>(-100., 100); },
    []() { return uniform_point<T, 2>(pi<T>); }
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> exponential{
    "exponential",
    [](const point<T, N>& p) {
      return -std::exp(
        -.5 *
        std::transform_reduce(
          std::begin(p), std::end(p), T{ 0. }, std::plus<T>{}, square<T>));
    },
    []() { return uniform_domain<T, N>(-1., 1.); },
    []() { return uniform_point<T, N>(0.); }
  };

  template<std::floating_point T>
  const test_function<T, 2> Goldstein_Price{
    "Goldstein-Price",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      const auto [x2, y2] = std::tuple<T, T>{ x * x, y * y };
      const auto xy = x * y;
      return (1. + square(x + y + 1.) *
                     (19. - 14. * x + 3. * x2 - 14. * y + 6. * xy + 3. * y2)) *
             (30. + square(2. * x - 3. * y) * (18. - 32. * x + 12. * x2 +
                                               48. * y - 36. * xy + 27. * y2));
    },
    []() { return uniform_domain<T, 2>(-2., 2.); },
    []() {
      return point<T, 2>{ 0., -1. };
    }
  };

  template<std::floating_point T>
  const test_function<T, 2> Hosaki{
    "Hosaki",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return (1. + x * (-8. + x * (7. + x * (-7. / 3. + x / 4.)))) * y * y *
             std::exp(-y);
    },
    []() { return uniform_domain<T, 2>(-10., 10.); },
    []() {
      return point<T, 2>{ 4., 2. };
    }
  };

  template<std::floating_point T>
  const test_function<T, 2> Leon{
    "Leon",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return 100. * square(y - x * x) + square(1. - x);
    },
    []() { return uniform_domain<T, 2>(-1.2, 1.2); },
    []() {
      return point<T, 2>{ 1., 1. };
    }
  };

  template<std::floating_point T>
  const test_function<T, 2> Matyas{
    "Matyas",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      return .26 * (x * x + y * y) - .48 * x * y;
    },
    []() { return uniform_domain<T, 2>(-10., 10.); },
    []() { return uniform_point<T, 2>(0.); }
  };

  template<std::floating_point T>
  const test_function<T, 2> Mexican_hat{
    "Mexican hat",
    [](const point<T, 2>& p) {
      const auto [x, y] = coordinates(p);
      const auto f = [&, x = x, y = y]() {
        return .1 + std::sqrt(square(x - 4.) + square(y - 4.));
      };
      return -20. * std::sin(f()) / f();
    },
    []() { return uniform_domain<T, 2>(-10., 10.); },
    []() { return uniform_point<T, 2>(4.); }
  };

  template<std::floating_point T>
  const test_function<T, 4> Miele_Cantrell{
    "Miele-Cantrell",
    [](const point<T, 4>& p) {
      return std::pow(std::exp(-p[0]) - p[1], 4.) +
             100. * std::pow(p[1] - p[2], 6.) +
             std::pow(std::tan(p[2] - p[3]), 4.) + std::pow(p[0], 8.);
    },
    []() { return uniform_domain<T, 4>(-1., 1); },
    []() {
      return point<T, 4>{ 0., 1., 1., 1. };
    }
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> Rosenbrock{
    "Rosenbrock",
    [](const point<T, N>& p) {
      T res = 0.;
      for (std::size_t i = 0; i < N - 1; ++i) {
        res += 100. * square(p[i + 1] - square(p[i])) + square(p[i] - 1.);
      }
      return res;
    },
    []() { return uniform_domain<T, N>(-30., 30.); },
    []() { return uniform_point<T, N>(1.); }
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> Schwefel{
    "Schwefel",
    [](const point<T, N>& p) {
      T res = 0.;
      for (T sum = 0.; auto x : p) {
        res += square(sum += x);
      }
      return res;
    },
    []() { return uniform_domain<T, N>(-100., 100.); },
    []() { return uniform_point<T, N>(0.); }
  };

  template<std::floating_point T, std::size_t N>
  const test_function<T, N> sphere{
    "sphere",
    [](const point<T, N>& p) {
      return std::transform_reduce(
        std::begin(p), std::end(p), T{ 0. }, std::plus<T>{}, square<T>);
    },
    []() { return uniform_domain<T, N>(0., 10.); },
    []() { return uniform_point<T, N>(0.); }
  };

} // namespace test_functions

} // namespace quile

#endif // QUILE_H
