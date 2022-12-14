/*
 * MIT License
 *
 * Copyright (c) 2020, 2021, 2022 Tomasz Tarkowski
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
#include <array>
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
 * @mainpage Introduction
 *
 * Quilë is a C++20 header-only general purpose genetic algorithms library with
 * no external dependencies supporting floating-point, integer, binary and
 * permutation representations. It is released under the terms of the MIT
 * License. Source code is available at https://github.com/ttarkowski/quile.
 *
 * The name of this library origins from fictional language Neo-Quenya and
 * means \em color.
 *
 * This work is a result of the project funded by National Science Centre Poland
 * (Twardowskiego 16, PL-30312 Kraków, Poland, http://www.ncn.gov.pl/) under the
 * grant number UMO-2016/23/B/ST3/03575.
 *
 * Example compilation command: `g++ -std=c++20 -DNDEBUG -O3 -Wall -Wextra
 * -pedantic -I/home/user/repos/quile -pthread program.cc`. Clang compilation
 * flags are identical. Please remove `-DNDEBUG` to enable assertions. Please
 * add `-DQUILE_ENABLE_LOGGING` to enable logging.
 *
 * Please note that examples from `doc/examples` directory were compiled with
 * following command `g++ -std=c++20 -DQUILE_ENABLE_LOGGING -O3 -Wall -Wextra
 * -pedantic -I../../ -pthread` and their output consisting of `std::cout` and
 * `std::cerr` streams were wrapped to fit 80 colums lines.
 */

/**
 * @file
 * `quile/quile.h` is one and only file of the library implementation.
 *
 * Please include `quile/quile.h` header using the code `#include
 * <quile/quile.h>` and compiling your program with apropriate include flag,
 * e.g. `-I/home/user/repos/quile` for GCC and Clang compilers.
 */

/**
 * @def QUILE_LOG(x)
 * `QUILE_LOG(x)` macro prints diagnostic information.
 *
 * `QUILE_LOG(x)` macro prints diagnostic information to the standard error
 * stream `std::cerr` provided that `QUILE_ENABLE_LOGGING` token is defined.
 * Otherwise it has no effect.
 *
 * Example:
 * @include QUILE_LOG.cc
 *
 * Result:
 * @verbinclude QUILE_LOG.out
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
 * @tparam T Loop index type.
 * @tparam I First index value (inclusive).
 * @tparam N Last index value (exclusive).
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
   * @include static_loop_body.cc
   *
   * Result:
   * @verbinclude static_loop_body.out
   */
  static void body([[maybe_unused]] auto&& f) {}
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
 * @include fn_and.cc
 *
 * Result (might be empty):
 * @verbinclude fn_and.out
 *
 * @note `fn_and` can be useful to describe complex genetic algorithm
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
 * @include fn_or.cc
 *
 * Result (might be empty):
 * @verbinclude fn_or.out
 *
 * @note `fn_or` can be useful to describe complex genetic algorithm
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
  {
  }

  /**
   * `thread_pool::async` asynchronically executes callable object `f`
   * postponing start of `f` until number of concurrently executing threads in
   * pool drops below number `sz`, described in constructor.
   *
   * @param policy Lauch policy (see `std::launch` documentation).
   * @param f Callable object to be concurrently executed.
   *
   * Example:
   * @include thread_pool_async.cc
   *
   * Result (might be different due to concurrent execution):
   * @verbinclude thread_pool_async.out
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
 * @tparam F Type convertible to the callable object type.
 * @tparam R Return type.
 * @tparam Args Argument type.
 *
 * Example:
 * @include callable.cc
 *
 * Result (might be empty):
 * @verbinclude callable.out
 */
template<typename F, typename R, typename... Args>
concept callable = std::convertible_to<std::invoke_result_t<F, Args...>, R>;

///////////
// Range //
///////////

/**
 * `range` represents subrange (closed interval) of type `T`.
 *
 * @tparam T Base type.
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
   * @throws std::invalid_argument Exception is raised if `min` is greater than
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
   * @note This constructor is available for types possesing
   * `std::numeric_limits` specialization.
   */
  template<typename U = T,
           typename = std::enable_if_t<std::numeric_limits<U>::is_specialized>>
  constexpr range()
    : range{ std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max() }
  {
  }

  /**
   * Default copy constructor `range::range`.
   */
  constexpr range(const range&) = default;

  /**
   * Default move constructor `range::range`.
   */
  constexpr range(range&&) = default;

  /**
   * Default assignment operator `range::operator=`.
   */
  range& operator=(const range&) = default;

  /**
   * Default move assignment operator `range::operator=`.
   */
  range& operator=(range&&) = default;

  /**
   * `range::min` returns range infimum.
   *
   * @returns Range infimum, i.e. left endpoint of interval \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ (`min` value) corresponding to the range.
   *
   * Example:
   * @include range_min.cc
   *
   * Result (might be empty):
   * @verbinclude range_min.out
   */
  T min() const { return min_; }

  /**
   * `range::max` returns range supremum.
   *
   * @returns Range supremum, i.e. right endpoint of interval \f$[{\rm min},
   * {\rm max}]_{\rm T}\f$ (`max` value) corresponding to the range.
   *
   * Example:
   * @include range_max.cc
   *
   * Result (might be empty):
   * @verbinclude range_max.out
   */
  T max() const { return max_; }

  /**
   * `range::midpoint` returns midpoint of interval represented by range.
   *
   * @returns Range midpoint, i.e. center of interval \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ (arithmetic mean of `min` and `max` values) corresponding
   * to the range. If `T` is of integer type and the sum of `min` and `max` is
   * odd, the result is rounded towards `min`.
   *
   * @note This method is disabled for ranges of type `bool`.
   *
   * Example:
   * @include range_midpoint.cc
   *
   * Result (might be empty):
   * @verbinclude range_midpoint.out
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
   * @returns Value clamped to the range.
   *
   * Example:
   * @include range_clamp.cc
   *
   * Result (might be empty):
   * @verbinclude range_clamp.out
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
   * @returns Boolean value of check: `true` if `t` is in \f$[{\rm min}, {\rm
   * max}]_{\rm T}\f$ interval and `false` otherwise.
   *
   * Example:
   * @include range_contains.cc
   *
   * Result (might be empty):
   * @verbinclude range_contains.out
   */
  bool contains(T t) const { return t >= min_ && t <= max_; }

  /**
   * `range::operator<=>` performs default lexicographical comparison with use
   * of left and right endpoints.
   *
   * @param r Range to be compared with `*this`.
   * @returns Ordering (cf. `std::strong_ordering`, `std::weak_ordering`,
   * `std::partial_ordering`).
   *
   * Example:
   * @include range_spaceship.cc
   *
   * Result (might be empty):
   * @verbinclude range_spaceship.out
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
 * @returns Reference to the `os` stream.
 *
 * Example:
 * @include range_stream.cc
 *
 * Result:
 * @verbinclude range_stream.out
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
 * @tparam T Number type.
 * @tparam N Returned container size.
 * @param t Starting value.
 * @returns `std::array<T, N>` with consecutive numbers starting from value `t`.
 *
 * Example:
 * @include iota.cc
 *
 * Result:
 * @verbinclude iota.out
 */
template<typename T, std::size_t N>
std::array<T, N>
iota(T t)
{
  std::array<T, N> res{};
  std::generate_n(std::begin(res), N, [t]() mutable { return t++; });
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
 *
 * @returns Reference to static object with Mersenne Twister engine
 * `std::mt19937` initialized with `std::random_device{}()`.
 *
 * Example:
 * @include random_engine.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude random_engine.out
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
 *
 * @param success_probability Probability of returning `true` value.
 * @returns Logic value drawn from \f${\rm B}(1, {\rm
 * success\_{}probability})\f$.
 *
 * Example:
 * @include success.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude success.out
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
 * @tparam T Result type (floating-point).
 * @param mean Mean of normal distribution.
 * @param standard_deviation Standard deviation of normal distribution.
 * @returns Number drawn from normal distribution.
 *
 * Example:
 * @include random_N.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude random_N.out
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
 * @tparam T Return type.
 * @param a Parameter describing aforementioned interval or set.
 * @param b Parameter describing aforementioned interval or set.
 * @returns Value drawn from uniform distribution.
 *
 * @note For floating-point types overflow may occur for `std::nextafter(b,
 * std::numeric_limits<T>::max()) - a` (cf. N4861, 26.6.8.2.2).
 *
 * Example:
 * @include random_U.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude random_U.out
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
 * @tparam T Argument and return type (floating-point or integer type).
 * @param x Argument to be raised to the second power.
 * @returns Argument raised to the second power, i.e. \f$x^2\f$.
 *
 * Example:
 * @include square.cc
 *
 * Result (might be empty):
 * @verbinclude square.out
 */
template<typename T>
requires std::floating_point<T> || std::integral<T> T
square(T x)
{
  return x * x;
}

/**
 * `cube` returns third power of its argument.
 *
 * @tparam T Argument and return type (floating-point or integer type).
 * @param x Argument to be raised to the third power.
 * @returns Argument raised to the third power, i.e. \f$x^3\f$.
 *
 * Example:
 * @include cube.cc
 *
 * Result (might be empty):
 * @verbinclude cube.out
 */
template<typename T>
requires std::floating_point<T> || std::integral<T> T
cube(T x)
{
  return x * x * x;
}

/**
 * `pi` is an approximation of \f$\pi\f$ number.
 *
 * @tparam T Floating-point type.
 *
 * Example:
 * @include pi_e_ln2.cc
 *
 * Result:
 * @verbinclude pi_e_ln2.out
 */
template<std::floating_point T>
const T pi = std::numbers::pi_v<T>;

/**
 * `e` is an approximation of \f$e\f$ number.
 *
 * @tparam T Floating-point type.
 *
 * Example:
 * @include pi_e_ln2.cc
 *
 * Result:
 * @verbinclude pi_e_ln2.out
 */
template<std::floating_point T>
const T e = std::numbers::e_v<T>;

/**
 * `ln2` is an approximation of \f$\ln 2\f$ number.
 *
 * @tparam T Floating-point type.
 *
 * Example:
 * @include pi_e_ln2.cc
 *
 * Result:
 * @verbinclude pi_e_ln2.out
 */
template<std::floating_point T>
const T ln2 = std::numbers::ln2_v<T>;

/**
 * `detail` contains library implementation details and is not inteded for use
 * by library end-user.
 */
namespace detail {

/**
 * `detail::angle` returns angle \f$\phi\f$ in polar coordinate system.
 *
 * @tparam T Argument and return type (floating-point).
 * @param x \f$x\f$ coordinate in Cartesian coordinate system.
 * @param y \f$y\f$ coordinate in Cartesian coordinate system.
 * @returns \f$\phi\f$ coordinate in polar coordinate system corresponding to
 * \f$(x, y)\f$ point.
 */
template<std::floating_point T>
T
angle(T x, T y)
{
  const T t = std::atan2(y, x);
  return t >= T{ 0 } ? t : (2 * std::numbers::pi_v<T> + t);
}

} // namespace detail

/**
 * `cart2spher` changes coordinate system from Cartesian to spherical.
 *
 * @tparam T Argument type and base for return type (floating-point).
 * @param x \f$x\f$ coordinate in Cartesian coordinate system.
 * @param y \f$y\f$ coordinate in Cartesian coordinate system.
 * @param z \f$y\f$ coordinate in Cartesian coordinate system.
 * @returns Tuple consisting of coordinates of \f$(r, \theta , \phi)\f$ point in
 * spherical coordinate system.
 *
 * @note \f$\theta \in [0, \pi]\f$, \f$\phi \in [0, 2\pi )\f$.
 *
 * Example:
 * @include cart2spher.cc
 *
 * Result:
 * @verbinclude cart2spher.out
 */
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

/**
 * `spher2cart` changes coordinate system from spherical to Cartesian.
 *
 * @tparam T Argument type and base for return type (floating-point).
 * @param r \f$r\f$ coordinate in spherical coordinate system.
 * @param theta \f$\theta\f$ coordinate in spherical coordinate system.
 * @param phi \f$\phi\f$ coordinate in spherical coordinate system.
 * @returns Tuple consisting of coordinates of \f$(x, y, z)\f$ point in
 * Cartesian coordinate system.
 *
 * @note \f$\theta \in [0, \pi]\f$, \f$\phi \in [0, 2\pi )\f$.
 *
 * Example:
 * @include cart2spher.cc
 *
 * Result:
 * @verbinclude cart2spher.out
 */
template<std::floating_point T>
std::tuple<T, T, T>
spher2cart(T r, T theta, T phi)
{
  return std::tuple<T, T, T>{ r * std::sin(theta) * std::cos(phi),
                              r * std::sin(theta) * std::sin(phi),
                              r * std::cos(theta) };
}

/**
 * `cart2polar` changes coordinate system from Cartesian to polar.
 *
 * @tparam T Argument type and base for return type (floating-point).
 * @param x \f$x\f$ coordinate in Cartesian coordinate system.
 * @param y \f$y\f$ coordinate in Cartesian coordinate system.
 * @returns Tuple consisting of coordinates of \f$(r, \phi)\f$ point in polar
 * coordinate system.
 *
 * @note \f$\phi \in [0, 2\pi )\f$.
 *
 * Example:
 * @include cart2polar.cc
 *
 * Result:
 * @verbinclude cart2polar.out
 */
template<std::floating_point T>
std::tuple<T, T>
cart2polar(T x, T y)
{
  return std::tuple<T, T>{ std::hypot(x, y), detail::angle(x, y) };
}

/**
 * `polar2cart` changes coordinate system from polar to Cartesian.
 *
 * @tparam T Argument type and base for return type (floating-point).
 * @param r \f$r\f$ coordinate in polar coordinate system.
 * @param phi \f$\phi\f$ coordinate in polar coordinate system.
 * @returns Tuple consisting of coordinates of \f$(x, y)\f$ point in Cartesian
 * coordinate system.
 *
 * @note \f$\phi \in [0, 2\pi )\f$.
 *
 * Example:
 * @include cart2polar.cc
 *
 * Result:
 * @verbinclude cart2polar.out
 */
template<std::floating_point T>
std::tuple<T, T>
polar2cart(T r, T phi)
{
  return std::tuple<T, T>{ r * std::cos(phi), r * std::sin(phi) };
}

////////////
// Domain //
////////////

/**
 * `domain` is a type representing domain in form of N-dimensional
 * parallelepiped.
 *
 * @tparam T Base type.
 * @tparam N Domain dimension.
 *
 * Example:
 * @include domain.cc
 *
 * Result (might be empty):
 * @verbinclude domain.out
 */
template<typename T, std::size_t N>
using domain = std::array<range<T>, N>;

/**
 * If `T` is some specialization of `domain` then `is_domain` provides member
 * constant `value` equal to `true`. Otherwise `value` is `false`.
 *
 * Example:
 * @include domain.cc
 *
 * Result (might be empty):
 * @verbinclude domain.out
 */
template<typename T>
struct is_domain : std::false_type
{
};

/**
 * Please see documentation for `is_domain<T>`.
 */
template<typename T, std::size_t N>
struct is_domain<domain<T, N>> : std::true_type
{
};

/**
 * `is_domain_v` is helper variable template for `is_domain`.
 *
 * Example:
 * @include domain.cc
 *
 * Result (might be empty):
 * @verbinclude domain.out
 */
template<typename T>
inline constexpr bool is_domain_v = is_domain<T>::value;

/**
 * `set_of_departure` specifies that `T` is some specialization of `domain`.
 *
 * Example:
 * @include domain.cc
 *
 * Result (might be empty):
 * @verbinclude domain.out
 */
template<typename T>
concept set_of_departure = is_domain_v<T>;

/**
 * `contains` checks if argument `p` is within domain `d` and returns `true` in
 * that case. Otherwise it returns `false`.
 *
 * @tparam T Domain base type.
 * @tparam N Domain dimensionality.
 * @param d Domain.
 * @param p Point to be checked.
 * @returns Boolean value describing whether point `p` is within domain `d`.
 *
 * Example:
 * @include contains.cc
 *
 * Result (might be empty):
 * @verbinclude contains.out
 */
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

/**
 * `uniform_domain` creates domain, where constraints on each direction are
 * identical, i.e. domain is of form of hypercube.
 *
 * @tparam T Domain base type.
 * @tparam N Domain dimensionality.
 * @param r Interval for hypercube construction.
 * @returns `N`-dimensional hypercube with edge of `r`.
 *
 * Example:
 * @include uniform_domain.cc
 *
 * Result (might be empty):
 * @verbinclude uniform_domain.out
 */
template<typename T, std::size_t N>
constexpr domain<T, N>
uniform_domain(const range<T>& r)
{
  domain<T, N> res{};
  std::generate_n(std::begin(res), N, [&]() { return r; });
  return res;
}

/**
 * `uniform_domain<T, N>(lo, hi)` is equivalent to `uniform_domain<T,
 * N>(range<T>{ lo, hi})` call. Please see documentation for `uniform_domain`
 * for argument of type `range`.
 *
 * Example:
 * @include uniform_domain.cc
 *
 * Result (might be empty):
 * @verbinclude uniform_domain.out
 */
template<typename T, std::size_t N>
constexpr domain<T, N>
uniform_domain(T lo, T hi)
{
  return uniform_domain<T, N>(range<T>{ lo, hi });
}

/**
 * `uniform` checks whether domain is of form of hypercube.
 *
 * @tparam T Domain base type.
 * @tparam N Domain dimensionality.
 * @param d Domain to be checked.
 * @returns Boolean value of check result.
 *
 * Example:
 * @include uniform_domain.cc
 *
 * Result (might be empty):
 * @verbinclude uniform_domain.out
 */
template<typename T, std::size_t N>
constexpr bool
uniform(const domain<T, N>& d)
{
  const auto x0 = N == 0 ? range<T>{} : d[0];
  return std::ranges::all_of(d, [&](const auto& x) { return x0 == x; });
}

/**
 * `self_adaptive_variation_domain` creates domain for self-adaptive mutation
 * based on domain for ordinary variation.
 *
 * @tparam T Domain base type.
 * @tparam N Domain dimensionality.
 * @param d Domain to use.
 * @param lo Minimum value for \f$\sigma_i\f$ on each direction.
 * @returns Domain with dimensionality of `2 * N`, where first part is equal to
 * `d` and the second one consits of ranges of form `range{ lo, 0.5 *
 * std::max(std::fabs(d[i].min()), std::fabs(d[i].max())) }`.
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
template<typename T, std::size_t N>
requires std::floating_point<T>
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

/**
 * `chain` represents genetic chain.
 *
 * Example:
 * @include chain.cc
 *
 * Result (might be empty):
 * @verbinclude chain.out
 */
template<typename T, std::size_t N>
using chain = std::array<T, N>;

/**
 * `chain_min` returns object of type `chain` filled at each `i` position with
 * `d[i].min()` value.
 *
 * @tparam T Chain base type.
 * @tparam N Chain length.
 * @param d Domain.
 * @returns Chain based on `d`.
 *
 * Example:
 * @include chain.cc
 *
 * Result (might be empty):
 * @verbinclude chain.out
 */
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

/**
 * `g_floating_point` specifies that `genotype` has floating-point
 * representation.
 *
 * @tparam T Floating-point type of representation.
 * @tparam N Genotype length.
 * @tparam D Pointer to the genotype domain.
 *
 * Example:
 * @include g_floating_point.cc
 *
 * Result (might be empty):
 * @verbinclude g_floating_point.out
 */
template<typename T, std::size_t N, const domain<T, N>* D>
requires std::floating_point<T>
struct g_floating_point
{
  static_assert(D != nullptr);
  static_assert(N > 0);

  /**
   * `g_floating_point::type` is floating-point type used for representing gene
   * values.
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  using type = T;

  /**
   * `size` returns domain size, i.e. `N`.
   *
   * @returns Genotype length (domain size).
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  static constexpr std::size_t size() { return N; }

  /**
   * `constraints` returns domain, i.e. `*D`.
   *
   * @returns Domain.
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  static constexpr const domain<type, size()>& constraints() { return *D; }

  /**
   * `g_floating_point::chain_t` is genetic chain type used as underlying
   * representation in `genotype`.
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  using chain_t = chain<type, size()>;

  /**
   * `valid` checks whether `c` belongs to the domain and returns `true` in that
   * case. Otherwise returns `false`.
   *
   * @param c Chain to be checked.
   * @returns Boolean value of check result.
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  static bool valid(const chain<type, size()>& c)
  {
    return contains(constraints(), c);
  }

  /**
   * `default_chain` returns chain filled in default way.
   *
   * @returns Default chain.
   *
   * Example:
   * @include g_floating_point.cc
   *
   * Result (might be empty):
   * @verbinclude g_floating_point.out
   */
  static chain_t default_chain() { return chain_min(constraints()); }
};

/**
 * If `T` is some specialization of `g_floating_point` then
 * `is_g_floating_point` provides member constant value equal to true. Otherwise
 * value is false.
 *
 * Example:
 * @include g_floating_point.cc
 *
 * Result (might be empty):
 * @verbinclude g_floating_point.out
 */
template<typename T>
struct is_g_floating_point : std::false_type
{
};

/**
 * Please see documentation for `is_g_floating_point<T>`.
 */
template<typename T, std::size_t N, const domain<T, N>* D>
struct is_g_floating_point<g_floating_point<T, N, D>> : std::true_type
{
};

/**
 * `is_g_floating_point_v` is helper variable template for
 * `is_g_floating_point`.
 *
 * Example:
 * @include g_floating_point.cc
 *
 * Result (might be empty):
 * @verbinclude g_floating_point.out
 */
template<typename T>
inline constexpr bool is_g_floating_point_v = is_g_floating_point<T>::value;

/**
 * `floating_point_representation` specifies that `T` is some specialization of
 * `g_floating_point`.
 *
 * Example:
 * @include g_floating_point.cc
 *
 * Result (might be empty):
 * @verbinclude g_floating_point.out
 */
template<typename T>
concept floating_point_representation = is_g_floating_point_v<T>;

/**
 * `g_integer` specifies that `genotype` has integer representation.
 *
 * @tparam T Integer type of representation (excluding Boolean).
 * @tparam N Genotype length.
 * @tparam D Pointer to the genotype domain.
 *
 * Example:
 * @include g_integer.cc
 *
 * Result (might be empty):
 * @verbinclude g_integer.out
 */
template<typename T, std::size_t N, const domain<T, N>* D>
requires std::integral<T> &&(!std::is_same_v<T, bool>)struct g_integer
{
  static_assert(D != nullptr);
  static_assert(N > 0);

  /**
   * `g_integer::type` is integer non-Boolean type used for representing gene
   * values.
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  using type = T;

  /**
   * `size` returns domain size, i.e. `N`.
   *
   * @returns Genotype length (domain size).
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  static constexpr std::size_t size() { return N; }

  /**
   * `constraints` returns domain, i.e. `*D`.
   *
   * @returns Domain.
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  static constexpr const domain<type, size()>& constraints() { return *D; }

  /**
   * `g_integer::chain_t` is genetic chain type used as underlying
   * representation in `genotype`.
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  using chain_t = chain<type, size()>;

  /**
   * `valid` checks whether `c` belongs to the domain and returns `true` in that
   * case. Otherwise returns `false`.
   *
   * @param c Chain to be checked.
   * @returns Boolean value of check result.
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  static bool valid(const chain<type, size()>& c)
  {
    return contains(constraints(), c);
  }

  /**
   * `default_chain` returns chain filled in default way.
   *
   * @returns Default chain.
   *
   * Example:
   * @include g_integer.cc
   *
   * Result (might be empty):
   * @verbinclude g_integer.out
   */
  static chain_t default_chain() { return chain_min(constraints()); }
};

/**
 * If `T` is some specialization of `g_integer` then `is_g_integer` provides
 * member constant value equal to true. Otherwise value is false.
 *
 * Example:
 * @include g_integer.cc
 *
 * Result (might be empty):
 * @verbinclude g_integer.out
 */
template<typename T>
struct is_g_integer : std::false_type
{
};

/**
 * Please see documentation for `is_g_integer<T>`.
 */
template<typename T, std::size_t N, const domain<T, N>* D>
struct is_g_integer<g_integer<T, N, D>> : std::true_type
{
};

/**
 * `is_g_integer_v` is helper variable template for `is_g_integer`.
 *
 * Example:
 * @include g_integer.cc
 *
 * Result (might be empty):
 * @verbinclude g_integer.out
 */
template<typename T>
inline constexpr bool is_g_integer_v = is_g_integer<T>::value;

/**
 * `integer_representation` specifies that `T` is some specialization of
 * `g_integer`.
 *
 * Example:
 * @include g_integer.cc
 *
 * Result (might be empty):
 * @verbinclude g_integer.out
 */
template<typename T>
concept integer_representation = is_g_integer_v<T>;

/**
 * `g_binary` specifies that `genotype` has binary representation.
 *
 * @tparam N Genotype length.
 *
 * Example:
 * @include g_binary.cc
 *
 * Result (might be empty):
 * @verbinclude g_binary.out
 */
template<std::size_t N>
struct g_binary
{
  static_assert(N > 0);

  /**
   * `g_binary::type` is Boolean type used for representing gene values.
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  using type = bool;

  /**
   * `size` returns domain size, i.e. `N`.
   *
   * @returns Genotype length (domain size).
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  static constexpr std::size_t size() { return N; }

  /**
   * `constraints` returns domain, i.e. \f$\{0, 1\}^N\f$.
   *
   * @returns Domain.
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  static constexpr const domain<type, size()> constraints()
  {
    return domain<type, size()>{};
  }

  /**
   * `g_binary::chain_t` is genetic chain type used as underlying
   * representation in `genotype`.
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  using chain_t = chain<type, size()>;

  /**
   * `valid` checks whether its argument belongs to the domain.
   *
   * @returns Boolean value of check result.
   *
   * @note Result is equal to `true` by definition.
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  static bool valid(const chain<type, size()>&) { return true; }

  /**
   * `default_chain` returns chain filled in default way.
   *
   * @returns Default chain.
   *
   * Example:
   * @include g_binary.cc
   *
   * Result (might be empty):
   * @verbinclude g_binary.out
   */
  static chain_t default_chain() { return chain_min(constraints()); }
};

/**
 * If `T` is some specialization of `g_binary` then `is_g_binary` provides
 * member constant value equal to true. Otherwise value is false.
 *
 * Example:
 * @include g_binary.cc
 *
 * Result (might be empty):
 * @verbinclude g_binary.out
 */
template<typename T>
struct is_g_binary : std::false_type
{
};

/**
 * Please see documentation for `is_g_binary<T>`.
 */
template<std::size_t N>
struct is_g_binary<g_binary<N>> : std::true_type
{
};

/**
 * `is_g_binary_v` is helper variable template for `is_g_binary`.
 *
 * Example:
 * @include g_binary.cc
 *
 * Result (might be empty):
 * @verbinclude g_binary.out
 */
template<typename T>
inline constexpr bool is_g_binary_v = is_g_binary<T>::value;

/**
 * `binary_representation` specifies that `T` is some specialization of
 * `g_binary`.
 *
 * Example:
 * @include g_binary.cc
 *
 * Result (might be empty):
 * @verbinclude g_binary.out
 */
template<typename T>
concept binary_representation = is_g_binary_v<T>;

/**
 * `g_permutation` specifies that `genotype` has permutation representation.
 *
 * @tparam T Integer type of representation (excluding Boolean).
 * @tparam N Genotype length.
 * @tparam M The lowest permuted number.
 *
 * @note Numbers from \f$\{M, \dots , M + N - 1\}\f$ set are permuted.
 *
 * Example:
 * @include g_permutation.cc
 *
 * Result (might be empty):
 * @verbinclude g_permutation.out
 */
template<typename T, std::size_t N, T M>
requires std::integral<T> &&(!std::is_same_v<T, bool>)struct g_permutation
{
  /**
   * `g_permutation::type` is integer non-Boolean type used for representing
   * gene values.
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  using type = T;

  /**
   * `size` returns domain size, i.e. `N`.
   *
   * @returns Genotype length (domain size).
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  static constexpr std::size_t size() { return N; }

  /**
   * `constraints` returns domain, i.e. \f$\{M, \dots , M + N - 1\}^N\f$.
   *
   * @returns Domain.
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  static constexpr const domain<type, size()> constraints()
  {
    return uniform_domain<type, size()>(M, M + N - 1);
  }

  /**
   * `g_permutation::chain_t` is genetic chain type used as underlying
   * representation in `genotype`.
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  using chain_t = chain<type, size()>;

  /**
   * `valid` checks whether `c` belongs to the domain (incl. check of the
   * permutation condition) and returns `true` in that case. Otherwise returns
   * `false`.
   *
   * @param c Chain to be checked.
   * @returns Boolean value of check result.
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  static bool valid(const chain<type, size()>& c)
  {
    const auto i = iota<type, size()>(M);
    return contains(constraints(), c) &&
           std::is_permutation(std::begin(c), std::end(c), std::begin(i));
  }

  /**
   * `default_chain` returns chain filled in default way.
   *
   * @returns Default chain.
   *
   * Example:
   * @include g_permutation.cc
   *
   * Result (might be empty):
   * @verbinclude g_permutation.out
   */
  static chain_t default_chain() { return iota<type, size()>(M); }
};

/**
 * If `T` is some specialization of `g_permutation` then `is_g_permutation`
 * provides member constant value equal to true. Otherwise value is false.
 *
 * Example:
 * @include g_permutation.cc
 *
 * Result (might be empty):
 * @verbinclude g_permutation.out
 */
template<typename T>
struct is_g_permutation : std::false_type
{
};

/**
 * Please see documentation for `is_g_permutation<T>`.
 */
template<typename T, std::size_t N, T M>
struct is_g_permutation<g_permutation<T, N, M>> : std::true_type
{
};

/**
 * `is_g_permutation_v` is helper variable template for `is_g_permutation`.
 *
 * Example:
 * @include g_permutation.cc
 *
 * Result (might be empty):
 * @verbinclude g_permutation.out
 */
template<typename T>
inline constexpr bool is_g_permutation_v = is_g_permutation<T>::value;

/**
 * `permutation_representation` specifies that `T` is some specialization of
 * `g_permutation`.
 *
 * Example:
 * @include g_permutation.cc
 *
 * Result (might be empty):
 * @verbinclude g_permutation.out
 */
template<typename T>
concept permutation_representation = is_g_permutation_v<T>;

/**
 * `chromosome_representation` specifies that `T` is some specialization of
 * one of allowed representations, i.e. `g_floating_point`, `g_integer`,
 * `g_binary` or `g_permutation`.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename T>
concept chromosome_representation =
  floating_point_representation<T> || integer_representation<T> ||
  binary_representation<T> || permutation_representation<T>;

/**
 * `genotype` is central type of the library---it allows genotype creation and
 * manipulation.
 *
 * @tparam R Type satisfying `chromosome_representation` concept.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename R>
requires chromosome_representation<R>
class genotype
{
public:
  /**
   * `genotype::chain_t` is genetic chain type (underlying genotype
   * representation).
   *
   * Example:
   * @include genotype_data.cc
   *
   * Result:
   * @verbinclude genotype_data.out
   */
  using chain_t = chain<typename R::type, R::size()>;

  /**
   * `genotype::const_iterator` is constant iterator to access underlying
   * representation.
   *
   * Example:
   * @include genotype_begin_end.cc
   *
   * Result:
   * @verbinclude genotype_begin_end.out
   */
  using const_iterator = typename chain_t::const_iterator;

  /**
   * `genotype::gene_t` is type of gene (e.g. floating-point).
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  using gene_t = typename R::type;

  /**
   * `genotype::genotype_t` is type of genotype representation dispatch tag.
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  using genotype_t = R;

  /**
   * `genotype::size` returns domain size.
   *
   * @returns Genotype length (domain size).
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  static constexpr std::size_t size() { return R::size(); }

  /**
   * `genotype::constraints` returns domain.
   *
   * @returns Domain.
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  static constexpr const domain<gene_t, size()> constraints()
  {
    return R::constraints();
  }

  /**
   * `genotype::uniform_domain` states whether `genotype` domain is uniform,
   * i.e. its domain is of form \f$X_0^N\f$.
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  static constexpr bool uniform_domain = uniform(constraints());

  /**
   * `genotype::valid` checks whether `c` belongs to the domain and returns
   * `true` in that case. Otherwise returns `false`.
   *
   * @param c Chain to be checked.
   * @returns Boolean value of check result.
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  static bool valid(const chain_t& c) { return R::valid(c); }

public:
  /**
   * `genotype::genotype` constructor creates object initialized with default
   * genetic chain.
   *
   * Example:
   * @include genotype.cc
   *
   * Result (might be empty):
   * @verbinclude genotype.out
   */
  genotype()
    : chain_{ R::default_chain() }
  {
  }

  /**
   * `genotype::genotype` constructor creates object initialized with genetic
   * chain passed as its argument.
   *
   * @param c Genetic chain to be used for initialization.
   *
   * @throws std::invalid_argument Exception is raised if `c` does not belong to
   * the domain or it does not fulfill condition of given representation (cf.
   * permutation condition).
   */
  explicit genotype(const chain_t& c)
    : chain_{ c }
  {
    if (!valid(c)) {
      throw std::invalid_argument{ "invalid chain" };
    }
  }

  /**
   * Default copy constructor `genotype::genotype`.
   */
  genotype(const genotype&) = default;

  /**
   * Default move constructor `genotype::genotype`.
   */
  genotype(genotype&&) = default;

  /**
   * Default assignment operator `genotype::operator=`.
   */
  genotype& operator=(const genotype&) = default;

  /**
   * Default move assignment operator `genotype::operator=`.
   */
  genotype& operator=(genotype&&) = default;

  /**
   * `genotype::value` returns gene value at \em locus i.
   *
   * @param i Gene \em locus.
   * @returns Gene value.
   *
   * Example:
   * @include genotype_value.cc
   *
   * Result:
   * @verbinclude genotype_value.out
   */
  gene_t value(std::size_t i) const { return chain_[i]; }

  /**
   * `genotype::value` changes gene value to `v` at \em locus `i`.
   *
   * @param i Gene \em locus.
   * @param v New gene value.
   * @returns Reference to `*this`.
   *
   * @throws std::invalid_argument Exception is raised if new gene value is
   * outside permitted interval for given \em locus `i`.
   *
   * @note This method is not available for permutation representation.
   *
   * Example:
   * @include genotype_value.cc
   *
   * Result:
   * @verbinclude genotype_value.out
   */
  template<typename S = R,
           typename = std::enable_if_t<!permutation_representation<S>>>
  genotype& value(std::size_t i, gene_t v)
  {
    if (!constraints()[i].contains(v)) {
      throw std::invalid_argument{ "bad value" };
    }
    chain_[i] = v;
    return *this;
  }

  /**
   * `genotype::random_reset` changes each gene value randomly using uniform
   * random distribution with intervals defined by domain.
   *
   * @returns Reference to `*this`.
   *
   * @note This overload is also available for permutation representation and
   * draws new permutation in that case.
   *
   * Example:
   * @include genotype_random_reset.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude genotype_random_reset.out
   */
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

  /**
   * `genotype::random_reset` changes gene at \em locus `i` value randomly using
   * uniform random distribution with interval defined by domain.
   *
   * @param i Gene \em locus.
   * @returns Reference to `*this`.
   *
   * @note This overload is not available for permutation representation.
   *
   * Example:
   * @include genotype_random_reset.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude genotype_random_reset.out
   */
  template<typename S = R,
           typename = std::enable_if_t<!permutation_representation<S>>>
  genotype& random_reset(std::size_t i)
  {
    const auto& c = constraints();
    chain_[i] = random_U<gene_t>(c[i].min(), c[i].max());
    return *this;
  }

  /**
   * `genotype::random` returns random genotype.
   *
   * @return Random genotype.
   *
   * Example:
   * @include genotype_random.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude genotype_random.out
   */
  static genotype random() { return genotype{}.random_reset(); }

  /**
   * `genotype::operator<=>` performs default lexicographical comparison with
   * use of genotypes' genetic chain.
   *
   * @param g Genotype to be compared with `*this`.
   * @returns Ordering (cf. `std::strong_ordering`, `std::weak_ordering`,
   * `std::partial_ordering`).
   *
   * @note Comparisons of genotypes with floating-point representation does not
   * include tolerance.
   *
   * Example:
   * @include genotype_cmp.cc
   *
   * Result (might be empty):
   * @verbinclude genotype_cmp.cc
   */
  auto operator<=>(const genotype& g) const { return chain_ <=> g.chain_; }

  /**
   * `genotype::operator==` performs comparison with use of genotypes' genetic
   * chain.
   *
   * @param g Genotype to be compared with `*this`.
   * @returns Boolean value `true` if chains are equal and false, otherwise.
   *
   * @note Comparisons of genotypes with floating-point representation does not
   * include tolerance.
   *
   * Example:
   * @include genotype_cmp.cc
   *
   * Result (might be empty):
   * @verbinclude genotype_cmp.cc
   */
  bool operator==(const genotype& g) const { return chain_ == g.chain_; }

  /**
   * `genotype::data` returns constant reference to the underlying genetic
   * chain.
   *
   * @returns Constant reference to the genetic chain.
   *
   * Example:
   * @include genotype_data.cc
   *
   * Result:
   * @verbinclude genotype_data.out
   */
  const chain_t& data() const { return chain_; }

  /**
   * `genotype::begin` returns constant iterator to the begin of genetic chain.
   *
   * @returns Constant iterator to the begin of genetic chain.
   *
   * Example:
   * @include genotype_begin_end.cc
   *
   * Result:
   * @verbinclude genotype_begin_end.out
   */
  const_iterator begin() const { return chain_.begin(); }

  /**
   * `genotype::end` returns constant iterator to the end of genetic chain.
   *
   * @returns Constant iterator to the end of genetic chain.
   *
   * @note The word \em end means past-the-last element.
   *
   * Example:
   * @include genotype_begin_end.cc
   *
   * Result:
   * @verbinclude genotype_begin_end.out
   */
  const_iterator end() const { return chain_.end(); }

private:
  chain_t chain_;
};

/**
 * If `T` is some specialization of `genotype` then `is_genotype`
 * provides member constant value equal to true. Otherwise value is false.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename T>
struct is_genotype : std::false_type
{
};

/**
 * Please see documentation for `is_genotype<T>`.
 */
template<typename T>
struct is_genotype<genotype<T>> : std::true_type
{
};

/**
 * `is_genotype_v` is helper variable template for `is_genotype`.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename T>
inline constexpr bool is_genotype_v = is_genotype<T>::value;

/**
 * `chromosome` specifies that `T` is some specialization of `genotype`.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename G>
concept chromosome = is_genotype_v<G>;

/**
 * `floating_point_chromosome` specifies that `T` is some floating-point type
 * specialization of `genotype`.
 *
 * Example:
 * @include genotype.cc
 *
 * Result (might be empty):
 * @verbinclude genotype.out
 */
template<typename G>
concept floating_point_chromosome =
  chromosome<G> && floating_point_representation<typename G::genotype_t>;

/**
 * `integer_chromosome` specifies that `T` is some integer type specialization
 * of `genotype`.
 *
 * Example:
 * @include integer_chromosome.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude integer_chromosome.out
 */
template<typename G>
concept integer_chromosome =
  chromosome<G> && integer_representation<typename G::genotype_t>;

/**
 * `binary_chromosome` specifies that `T` is some binary type specialization
 * of `genotype`.
 *
 * Example:
 * @include binary_chromosome.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude binary_chromosome.out
 */
template<typename G>
concept binary_chromosome =
  chromosome<G> && binary_representation<typename G::genotype_t>;

/**
 * `permutation_chromosome` specifies that `T` is some permutation type
 * specialization of `genotype`.
 *
 * Example:
 * @include permutation_chromosome.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude permutation_chromosome.out
 */
template<typename G>
concept permutation_chromosome =
  chromosome<G> && permutation_representation<typename G::genotype_t>;

/**
 * `uniform_chromosome` specifies that `T` is some specialization of `genotype`
 * satisfying uniformity condition.
 *
 * Example:
 * @include permutation_chromosome.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude permutation_chromosome.out
 */
template<typename G>
concept uniform_chromosome = chromosome<G> && G::uniform_domain;

/**
 * `genotype_constraints` specifies that `F` is some predicate that states which
 * genotypes are proper.
 *
 * @note Genotype \f$g \in X_0 \times \dots \times X_{c-1} = \prod_{i=0}^{c-1}
 * X_i\f$, where \f$X_i\f$ is equal to \f$\mathbb{B} = \{{\rm false}, {\rm
 * true}\}\f$ or is bounded subset of set of real numbers \f$\mathbb{R}\f$ or
 * integer numbers \f$\mathbb{Z}\f$.
 *
 * @note Proper genotype \f$g \in G = \{(x_0, \dots , x_{c-1}) \in
 * \prod_{i=0}^{c-1} X_i \mid Q(x_0, \dots , x_{c-1}) \}\f$, where \f$Q\f$ is
 * predicate stating which genotype is proper.
 *
 * @note Set defined with use of predicate is called \em extension \em of \em
 * predicate.
 *
 * Example:
 * @include genotype_constraints.cc
 *
 * Result (might be empty):
 * @verbinclude genotype_constraints.out
 */
template<typename F, typename G>
concept genotype_constraints = std::predicate<F, G> && chromosome<G>;

/**
 * `constraints_satisfied` is predicate stating that all genotypes are proper.
 *
 * Example:
 * @include genotype_constraints.cc
 *
 * Result (might be empty):
 * @verbinclude genotype_constraints.out
 */
template<typename G>
requires chromosome<G>
const auto constraints_satisfied = [](const G&) { return true; };

/**
 * `operator<<` prints genotype to the stream.
 *
 * @param os Stream to use.
 * @param g Genotype to be printed.
 * @returns Reference to the `os` stream.
 *
 * Example:
 * @include genotype_stream.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude genotype_stream.out
 */
template<typename G>
requires chromosome<G> std::ostream&
operator<<(std::ostream& os, const G& g)
{
  for (std::size_t i = 0; i < G::size(); ++i) {
    os << g.value(i) << (i + 1 < G::size() ? " " : "");
  }
  return os;
}

/**
 * `operator<<` prints genotype to the stream.
 *
 * @param os Stream to use.
 * @param g Genotype to be printed.
 * @returns Reference to the `os` stream.
 *
 * @note This overload is dedicated for floating-point genotypes.
 *
 * Example:
 * @include genotype_stream.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude genotype_stream.out
 */
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

/**
 * Namespace `std` is opened only for the purpose of injecting `hash` into it.
 */
namespace std {

/**
 * `std::hash` for `genotype`.
 *
 * @note This `std::hash` specialization allows interoperability with STL
 * unordered associative containers like `std::unordered_map`.
 *
 * @note `std::hash` is injected into `std` namespace of programs using this
 * library.
 */
template<typename G>
requires quile::chromosome<G>
struct hash<G>
{
  /**
   * `std::hash::operator()` calculates hash function value for genotype `g`.
   *
   * @param g Genotype.
   * @returns Hash function value.
   */
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

} // namespace std

namespace quile {

////////////////
// Population //
////////////////

/**
 * `population` is a sequence of genotypes implemented as a `std::vector`
 * sequence container.
 *
 * Example:
 * @include population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude population.out
 */
template<typename G>
requires chromosome<G>
using population = std::vector<G>;

/**
 * If `T` is some specialization of `population` then `is_population`
 * provides member constant value equal to true. Otherwise value is false.
 *
 * Example:
 * @include population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude population.out
 */
template<typename T>
struct is_population : std::false_type
{
};

/**
 * Please see documentation for `is_population<T>`.
 */
template<typename G>
struct is_population<population<G>> : std::true_type
{
};

/**
 * `is_population_v` is helper variable template for `is_population`.
 *
 * Example:
 * @include population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude population.out
 */
template<typename G>
inline constexpr bool is_population_v = is_population<G>::value;

/**
 * `genetic_pool` specifies that `T` is some specialization of `population`.
 *
 * Example:
 * @include population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude population.out
 */
template<typename G>
concept genetic_pool = is_population_v<G>;

/**
 * `populate_0_fn` can be used for first generation creation.
 *
 * Example:
 * @include random_population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude random_population.out
 */
template<typename G>
requires chromosome<G>
using populate_0_fn = std::function<population<G>(std::size_t)>;

/**
 * `populate_1_fn` can be used for parents selection.
 *
 * Example:
 * @include selection.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude selection.out
 */
template<typename G>
requires chromosome<G>
using populate_1_fn =
  std::function<population<G>(std::size_t, const population<G>&)>;

/**
 * `populate_2_fn` can be used for selection to the next generation.
 *
 * Example:
 * @include selection.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude selection.out
 */
template<typename G>
requires chromosome<G>
using populate_2_fn = std::function<
  population<G>(std::size_t, const population<G>&, const population<G>&)>;

/**
 * `generations` is a sequence of populations.
 *
 * @note For memory optimization purpose the `std::deque` was chosen instead of
 * `std::vector`---this implementation allows for erasing of the oldest
 * generation.
 *
 * Example:
 * @include generations.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude generations.out
 */
template<typename G>
requires chromosome<G>
using generations = std::deque<population<G>>;

//////////////////////////////
// Mutation & recombination //
//////////////////////////////

/**
 * `mutation` specifies that `M` instance applied to `genotype` returns
 * object convertible to `population`.
 *
 * Example:
 * @include mutation.cc
 *
 * Result:
 * @verbinclude mutation.out
 */
template<typename M, typename G>
concept mutation = requires(M m, G g)
{
  {
    m(g)
    } -> std::convertible_to<population<G>>;
}
&&chromosome<G>;

/**
 * `mutation_fn` is a callable object which can be invoked on `genotype` and
 * returns `population`.
 *
 * Example:
 * @include mutation.cc
 *
 * Result:
 * @verbinclude mutation.out
 */
template<typename G>
requires chromosome<G>
using mutation_fn = std::function<population<G>(const G&)>;

/**
 * `recombination` specifies that `M` instance applied to two objects of type
 * `genotype` returns object convertible to `population`.
 *
 * Example:
 * @include recombination.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude recombination.out
 */
template<typename R, typename G>
concept recombination = requires(R r, G g)
{
  {
    r(g, g)
    } -> std::convertible_to<population<G>>;
}
&&chromosome<G>;

/**
 * `recombination_fn` is a callable object which can be invoked on two objects
 * of type `genotype` and returns `population`.
 *
 * Example:
 * @include recombination.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude recombination.out
 */
template<typename G>
requires chromosome<G>
using recombination_fn = std::function<population<G>(const G&, const G&)>;

/**
 * `unary_identity` is an identity mutation.
 *
 * @tparam G Some `genotype` specialization.
 * @param g Genotype.
 * @returns Population consisting of `g`.
 *
 * Example:
 * @include identity.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude identity.out
 */
template<typename G>
requires chromosome<G> population<G>
unary_identity(const G& g)
{
  return population<G>{ g };
}

/**
 * `binary_identity` is an identity recombination.
 *
 * @tparam G Some `genotype` specialization.
 * @param g0 Genotype (first parent).
 * @param g1 Genotype (second parent)
 * @returns Population consisting of `g0` and `g1`.
 *
 * Example:
 * @include identity.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude identity.out
 */
template<typename G>
requires chromosome<G> population<G>
binary_identity(const G& g0, const G& g1)
{
  return population<G>{ g0, g1 };
}

/**
 * `variation` represents variation operator.
 *
 * @tparam G Some `genotype` specialization.
 *
 * @note At the moment library supports only canonical forms of variations
 * (unary and binary).
 *
 * Example:
 * @include variation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude variation.out
 */
template<typename G>
requires chromosome<G>
class variation
{
public:
  /**
   * `variation::variation` constructor creates object representing variation
   * consisting of recombination `r` with mutation `m` applied separately to
   * each child coming from `r`.
   *
   * @param m Mutation.
   * @param r Recombination.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
  variation(const mutation_fn<G>& m, const recombination_fn<G>& r)
    : m_{ m }
    , r_{ r }
  {
  }

  /**
   * `variation::variation` constructor creates identity variation.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
  variation()
    : variation{ unary_identity<G>, binary_identity<G> }
  {
  }

  /**
   * `variation::variation` constructor creates variation equal to mutation `m`.
   *
   * @param m Mutation.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
  explicit variation(const mutation_fn<G>& m)
    : variation{ m, binary_identity<G> }
  {
  }

  /**
   * `variation::variation` constructor creates variation equal to recombination
   * `r`.
   *
   * @param r Recombination.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
  explicit variation(const recombination_fn<G>& r)
    : variation{ unary_identity<G>, r }
  {
  }

  /**
   * `variation::operator()` applies variation to genotypes `g0` and `g1`.
   *
   * @param g0 Genotype.
   * @param g1 Genotype.
   * @returns Population resulting from application of variation to genotypes.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
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

  /**
   * `variation::operator()` applies variation to consecutive pairs of genotypes
   * in population `p`.
   *
   * @param p Population consisting of pairs of parents.
   * @returns Populative consisting of cumulative offspring.
   *
   * @throws std::invalid_argument Exception is raised if population size is
   * odd.
   *
   * Example:
   * @include variation.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude variation.out
   */
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

/**
 * `stochastic_mutation` creates stochastic mutation consisting of `m` applied
 * with probability `p`.
 *
 * @tparam G Some `genotype` specialization.
 * @param m Mutation.
 * @param p Probability.
 * @returns Stochastic mutation.
 *
 * Example:
 * @include stochastic_variation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude stochastic_variation.out
 */
template<typename G>
requires chromosome<G>
auto
stochastic_mutation(const mutation_fn<G>& m, probability p)
{
  return [=](const G& g) { return success(p) ? m(g) : population<G>{ g }; };
}

/**
 * `stochastic_recombination` creates stochastic recombination consisting of `r`
 * applied with probability `p`.
 *
 * @tparam G Some `genotype` specialization.
 * @param r Recombination.
 * @param p Probability.
 * @returns Stochastic recombination.
 *
 * Example:
 * @include stochastic_variation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude stochastic_variation.out
 */
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

/**
 * `termination_condition` specifies that `F` is some predicate that states when
 * evolution should be finished.
 */
template<typename F, typename G>
concept termination_condition = std::predicate<F, std::size_t, generations<G>>;

/**
 * `termination_condition_fn` is a callable object which states when evolution
 * should be finished.
 */
template<typename G>
requires chromosome<G>
using termination_condition_fn =
  std::function<bool(std::size_t, const generations<G>&)>;

/**
 * `evolution` executes evolutionary process.
 *
 * @tparam G Some `genotype` specialization.
 * @param v Variation.
 * @param first_generation First generation.
 * @param p1 Parents selection mechanism.
 * @param p2 Selection to the next generation mechanism.
 * @param tc Termination condition.
 * @param parents_sz Size of the parents multiset (should be even).
 * @param max_history Number of generations kept in memory and returned to the
 * caller. Default zero value is special and means keeping and returning all
 * generations.
 * @returns Generations produced during evolution (cf. `max_history` argument).
 */
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

/**
 * `evolution` executes evolutionary process.
 *
 * @tparam G Some `genotype` specialization.
 * @param v Variation.
 * @param p0 Mechanism for first generation creation.
 * @param p1 Parents selection mechanism.
 * @param p2 Selection to the next generation mechanism.
 * @param tc Termination condition.
 * @param generation_sz Generation size.
 * @param parents_sz Size of the parents multiset (should be even).
 * @param max_history Number of generations kept in memory and returned to the
 * caller. Default zero value is special and means keeping and returning all
 * generations.
 * @returns Generations produced during evolution (cf. `max_history` argument).
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
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

/**
 * `fitness` is a fitness function value type.
 *
 * Example:
 * @include selection.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude selection.out
 */
using fitness = double;

/**
 * `fitnesses` is a sequential container of fitness values.
 *
 * Example:
 * @include fitness_db.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude fitness_db.out
 */
using fitnesses = std::vector<fitness>;

/**
 * `fitness_function` describes how fit genotype is.
 *
 * @note From implementation point of view the `fitness_function` should be
 * thread-safe.
 *
 * Example:
 * @include fitness_db.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude fitness_db.out
 */
template<typename G>
requires chromosome<G>
using fitness_function = std::function<fitness(const G&)>;

/**
 * `incalculable` is a special value which can be used when given gentotype is
 * not proper (cf. `genotype_constraints`) or to signal some problem in
 * `fitness_function` (e.g. non-convergence).
 */
const fitness incalculable = -std::numeric_limits<fitness>::infinity();

/**
 * `fitness_db` is an intermediary object to fitness function values database.
 *
 * @tparam G Some `genotype` specialization.
 *
 * @note Intermediary objects own database through the `std::shared_ptr`.
 *
 * Example:
 * @include fitness_db.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude fitness_db.out
 */
template<typename G>
requires chromosome<G>
class fitness_db
{
private:
  /**
   * `fitness_db::database` is an unordered map with genotypes as its keys and
   * fitness function values as its values.
   */
  using database = std::unordered_map<G, fitness>;

public:
  /**
   * `fitness_db::const_iterator` is a constant iterator to the underlying
   * `database` container.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
  using const_iterator = typename database::const_iterator;

public:
  /**
   * `fitness_db::fitness_db` constructor creates intermediary object to fitness
   * function values database.
   *
   * @param f Fitness function.
   * @param gc Predicate defining proper genotypes.
   * @param thread_sz Number of threads for concurrent fitness function values
   * calculations. Default value is equal to
   * `std::thread::hardware_concurrency()`.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
  explicit fitness_db(
    const fitness_function<G>& f,
    const genotype_constraints<G> auto& gc,
    unsigned int thread_sz = std::thread::hardware_concurrency())
    : function_{ [=](const G& g) { return gc(g) ? f(g) : incalculable; } }
    , thread_sz_{ thread_sz }
  {
  }

  /**
   * Default copy constructor `fitness_db::fitness_db`.
   */
  fitness_db(const fitness_db&) = default;

  /**
   * Default assignment operator `fitness_db::operator=`.
   */
  fitness_db& operator=(const fitness_db&) = default;

  /**
   * `fitness_db::operator()` returns fitness function value for genotype `g`
   * from database. If the value is not yet available, it is calculated,
   * inserted to the database and then returned to the caller.
   *
   * @param g Genotype for which fitness function value is needed.
   * @returns Fitness function value for genotype `g`.
   *
   * @note This method is non-concurrent.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
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

  /**
   * `fitness_db::operator()` returns fitness function values for genotypes from
   * population from database. If some values are not yet available, they are
   * calculated, inserted to the database and then all values are returned to
   * the caller.
   *
   * @param p Population for which fitness function values are needed.
   * @returns Fitness function values for genotypes from population `p` in order
   * corresponding to the order of genotypes in population itself.
   *
   * @note This method is potentially concurrent.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
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

  /**
   * `fitness_db::size` returns number of keys (genotypes) in database.
   *
   * @returns Number of database keys.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
  std::size_t size() const { return fitness_values_->size(); }

  /**
   * `fitness_db::begin` returns constant iterator to the begin of database.
   *
   * @returns Constant iterator to the begin of database.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
  const_iterator begin() const { return fitness_values_->begin(); }

  /**
   * `fitness_db::end` returns constant iterator to the end of database.
   *
   * @returns Constant iterator to the end of database.
   *
   * @note The word \em end means past-the-last element.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
  const_iterator end() const { return fitness_values_->end(); }

  /**
   * `fitness_db::rank_order` returns all genotypes from database in descending
   * order of fitness function value.
   *
   * @returns Population consisting of all genotypes from database.
   *
   * @note `rank_order()[0]` gives the best genotype for non-empty database.
   *
   * Example:
   * @include fitness_db.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude fitness_db.out
   */
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

/**
 * `print` prints to the stream `os` information about each genotype from each
 * generation accompanied with optional information about fitness function
 * value.
 *
 * @tparam G Some `genotype` specialization.
 * @param os Stream to print on.
 * @param gs Generations.
 * @param fd Pointer to the fitness function values database. Default value is
 * equal to `nullptr`.
 *
 * @note If `fd` is equal to `nullptr` then information about fitness function
 * values is not printed.
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
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

/**
 * `selection_probabilities` is a sequential container intended for keeping
 * selection probabilities values.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
using selection_probabilities = std::vector<probability>;

/**
 * `selection_probabilities_fn` is a callable object which can be invoked on
 * population and returns corresponding selection probabilities for each
 * genotype from population.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
template<typename G>
requires chromosome<G>
using selection_probabilities_fn =
  std::function<selection_probabilities(const population<G>&)>;

/**
 * `cumulative_probabilities` serves for calculation of cumulative selection
 * probabilities, which can be used later in roulette wheel selection or
 * stochastic universal sampling.
 *
 * @tparam G Some `genotype` specialization.
 * @param spf Selection probabilities function.
 * @param p Population.
 * @returns Cumulative selection probabilities.
 *
 * @note \em Cumulative selection probability for index `i` means sum of
 * selection probabilities from index `0` up to `i`. Cumulative selection
 * probability for last genotype in population is by definition equal to `1`.
 *
 * @note It is guaranteed that last cumulative probability in returned
 * container is equal to `1` without any numerical error.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
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

/**
 * `select_different_than` returns copy of container `c` with deleted elements
 * equal to `t`.
 *
 * @tparam C Container type.
 * @tparam T Type convertible to C container element type.
 * @param c Container.
 * @param t Element to be dropped in returned container.
 * @param require_nonempty_result Flag for possible exception.
 * @returns Container with values different than `t`.
 *
 * @throws std::runtime_error Exception is raised if `require_nonempty_result`
 * is `true` and returned container is empty.
 */
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

/**
 * `select_calculable` returns container with fitness function values equal to
 * `fs`, but with `incalculable` entries skipped.
 *
 * @param fs Fitness function values container.
 * @param require_nonempty_result Flag for possible exception.
 * @returns Container with values different than `incalculable`.
 *
 * @throws std::runtime_error Exception is raised if `require_nonempty_result`
 * is `true` and returned container is empty.
 */
inline fitnesses
select_calculable(const fitnesses& fs, bool require_nonempty_result = false)
{
  return select_different_than(fs, incalculable, require_nonempty_result);
}

/**
 * `fitness_proportional_selection` is fitness proportional selection (a.k.a.
 * fitness \em proportionate selection) with windowing procedure (FPS).
 *
 * @tparam G Some `genotype` specialization.
 *
 * @note This implementation has workarounds for population of equally fit
 * genotypes and populations containing genotypes which fitnesses cannot be
 * calculated. Please note that there should be at least one genotype, which
 * fitness can be calculated.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
template<typename G>
class fitness_proportional_selection
{
public:
  /**
   * `fitness_proportional_selection::fitness_proportional_selection`
   * constructor creates FPS mechanism for database represented by intermediary
   * object `ff`.
   *
   * @param ff Fitness database intermediary object.
   *
   * Example:
   * @include spf.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude spf.out
   */
  explicit fitness_proportional_selection(const fitness_db<G>& ff)
    : ff_{ ff }
  {
  }

  /**
   * `fitness_proportional_selection::operator()` returns selection
   * probabilities for population `p`.
   *
   * @param p Population.
   * @returns FPS selection probabilities for population `p`.
   *
   * @throws std::runtime_error Exception is raised if fitness function
   * evaluates to `incalculable` for all genotypes from `p`.
   *
   * Example:
   * @include spf.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude spf.out
   */
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

/**
 * `detail::advance_cpy` wraps `std::advance` and returns iterator.
 *
 * @tparam It Iterator type.
 * @param it Iterator.
 * @param n Distance to advance.
 * @returns Advanced iterator.
 */
template<typename It>
It
advance_cpy(It it, std::size_t n)
{
  std::advance(it, n);
  return it;
}

/**
 * `detail::rank` returns ranking position of element in range after stable
 * sort.
 *
 * @tparam It Iterator type.
 * @tparam Compare Type of comparison mechanism.
 * @param first Range begin.
 * @param last Range end.
 * @param comp Comparison function object.
 * @returns Sequence container with ranking positions.
 */
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

/**
 * `detail::id` performs identity operation between types.
 *
 * @tparam T Destination type.
 * @tparam U Source type.
 * @param u Argument.
 * @returns Argument after transformation to destination type.
 */
template<typename T, typename U>
T
id(U u)
{
  return u;
}

} // namespace detail

/**
 * `linear_ranking_selection` creates linear pressure mechanism for ranking
 * selection.
 *
 * @param s So-called \em s parameter.
 * @returns Linear pressure mechanism.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
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

/**
 * `exponential_ranking_selection` is exponential pressure mechanism for ranking
 * selection.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
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

/**
 * `ranking_selection` is ranking selection (RS) mechanism.
 *
 * @note This implementation has workarounds for populations containing
 * genotypes which fitnesses cannot be calculated. Please note that there should
 * be at least one genotype, which fitness can be calculated.
 *
 * Example:
 * @include spf.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude spf.out
 */
template<typename G>
class ranking_selection
{
private:
  using probability_fn = std::function<probability(std::size_t, std::size_t)>;

public:
  /**
   * `ranking_selection::ranking_selection` constructor creates RS.
   *
   * @param ff Fitness function database intermediary object.
   * @param pf Selection pressure mechanism (linear or exponential).
   *
   * Example:
   * @include spf.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude spf.out
   */
  ranking_selection(const fitness_db<G>& ff, const probability_fn& pf)
    : ff_{ ff }
    , pf_{ pf }
  {
  }

  /**
   * `ranking_selection::operator()` returns selection probabilities for
   * population `p`.
   *
   * @param p Population.
   * @returns RS selection probabilities for population `p`.
   *
   * @throws std::runtime_error Exception is raised if fitness function
   * evaluates to `incalculable` for all genotypes from `p`.
   *
   * Example:
   * @include spf.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude spf.out
   */
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

/**
 * `detail::generate` creates population of size `lambda` filled with result of
 * function `f`.
 *
 * @tparam G Some `genotype` specialization.
 * @param lambda Result population size.
 * @param f Function returning genotype.
 * @returns Population of size `lambda` where each member is result of `f`.
 */
template<typename G>
requires chromosome<G> population<G>
generate(std::size_t lambda, const std::function<G()>& f)
{
  population<G> res{};
  std::generate_n(std::back_inserter(res), lambda, f);
  return res;
}

} // namespace detail

/**
 * `adapter` implements \em flatten function between `populate_2_fn` mechanism
 * and `populate_1_fn` mechanism.
 *
 * @tparam G Some `genotype` specialization.
 * @param fn Mechanism of `populate_1_fn` type.
 * @returns Mechanism of `populate_2_fn` type, which applies `fn` to two
 * flattened populations.
 *
 * @note `adapter` can be useful when used with roulette wheel selection or
 * stochastic universal sampling as a method of selection to the next
 * generation.
 *
 * Example:
 * @include selection.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude selection.out
 */
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

/**
 * `random_population` returns random population of size `lambda`, where each
 * member genotype satisfies predicate `C`.
 *
 * @tparam C Proper genotype predicate.
 * @tparam G Some `genotype` specialization.
 * @param lambda Size of returned population.
 * @returns Random population.
 *
 * Example:
 * @include random_population.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude random_population.out
 */
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

/**
 * `roulette_wheel_selection` is roulette wheel selection (a.k.a. roulette wheel
 * \em algorithm, RWA).
 *
 * @tparam G Some `genotype` specialization.
 */
template<typename G>
requires chromosome<G>
class roulette_wheel_selection
{
public:
  /**
   * `roulette_wheel_selection::roulette_wheel_selection` constructor creates
   * RWA with selection probability function `spf`.
   *
   * @param spf Selection probability function.
   */
  explicit roulette_wheel_selection(const selection_probabilities_fn<G>& spf)
    : spf_{ spf }
  {
  }

  /**
   * `roulette_wheel_selection::operator()` draws `lambda` genotypes from
   * population `p` according to the RWA.
   *
   * @param lambda Size of the returned population.
   * @param p Source population.
   * @returns Population consisting of genotypes drawn from `p`.
   */
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

/**
 * `stochastic_universal_sampling` is stochastic universal sampling (SUS).
 *
 * @tparam G Some `genotype` specialization.
 *
 * Example:
 * @include selection.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude selection.out
 */
template<typename G>
requires chromosome<G>
class stochastic_universal_sampling
{
public:
  /**
   * `stochastic_universal_sampling::stochastic_universal_sampling` constructor
   * creates SUS with selection probability function `spf`.
   *
   * @param spf Selection probability function.
   *
   * Example:
   * @include selection.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude selection.out
   */
  explicit stochastic_universal_sampling(
    const selection_probabilities_fn<G>& spf)
    : spf_{ spf }
  {
  }

  /**
   * `stochastic_universal_sampling::operator()` draws `lambda` genotypes from
   * population `p` according to the SUS.
   *
   * @param lambda Size of the returned population.
   * @param p Source population.
   * @returns Population consisting of genotypes drawn from `p`.
   *
   * Example:
   * @include selection.cc
   *
   * Result (might be different due to randomness):
   * @verbinclude selection.out
   */
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

/**
 * `generational_survivor_selection` is generational survivor selection
 * mechanism.
 *
 * @tparam G Some `genotype` specialization.
 * @param sz Generation / offspring size.
 * @param generation Current generation.
 * @param offspring Offspring.
 * @returns Offspring.
 *
 * @throws std::invalid_argument Exception is raised if `generation` size or
 * `offspring` size is different from `sz`.
 */
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

/**
 * `max` returns maximum fitness function value from `fs`.
 *
 * @param fs Fitness function values container.
 * @returns Maximum value.
 *
 * @throws std::runtime_error Exception is raised if all `fs` elements are equal
 * to `incalculable`.
 */
inline fitness
max(const fitnesses& fs)
{
  const fitnesses calc{ select_calculable(fs, true) };
  return *std::ranges::max_element(calc);
}

/**
 * `max` returns maximum fitness function value for population `p` and database
 * intermediary object `ff`.
 *
 * @tparam G Some `genotype` specialization.
 * @param p Population.
 * @param ff Database intermediary object.
 * @returns Maximum value.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from `p`.
 */
template<typename G>
requires chromosome<G> fitness
max(const population<G>& p, const fitness_db<G>& ff)
{
  return max(ff(p));
}

/**
 * `max` returns maximum fitness function values for each generation from `gs`
 * and database intermediary object `ff`.
 *
 * @tparam G Some `genotype` specialization.
 * @param gs Sequence of generations.
 * @param ff Database intermediary object.
 * @returns Maximum values corresponding to each generation.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from at least one generation.
 */
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

/**
 * `min` returns minimum fitness function value from `fs`.
 *
 * @param fs Fitness function values container.
 * @returns Minimum value.
 *
 * @throws std::runtime_error Exception is raised if all `fs` elements are equal
 * to `incalculable`.
 */
inline fitness
min(const fitnesses& fs)
{
  const fitnesses calc{ select_calculable(fs, true) };
  return *std::ranges::min_element(calc);
}

/**
 * `min` returns minimum fitness function value for population `p` and database
 * intermediary object `ff`.
 *
 * @tparam G Some `genotype` specialization.
 * @param p Population.
 * @param ff Database intermediary object.
 * @returns Minimum value.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from `p`.
 */
template<typename G>
requires chromosome<G> fitness
min(const population<G>& p, const fitness_db<G>& ff)
{
  return min(ff(p));
}

/**
 * `min` returns minimum fitness function values for each generation from `gs`
 * and database intermediary object `ff`.
 *
 * @tparam G Some `genotype` specialization.
 * @param gs Sequence of generations.
 * @param ff Database intermediary object.
 * @returns Minimum values corresponding to each generation.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from at least one generation.
 */
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

/**
 * `max_iterations_termination` returns condition, which terminates algorithm
 * after performing `max` loop iterations.
 *
 * @tparam G Some `genotype` specialization.
 * @param max Number of genetic algorithm loop iteration to perform (number of
 * generations).
 * @returns Predicate terminating genetic algorithm after `max` iterations.
 */
template<typename G>
termination_condition_fn<G>
max_iterations_termination(std::size_t max)
{
  return [=](std::size_t i, const generations<G>&) { return i == max; };
}

/**
 * `max_fitness_improvement_termination` returns condition, which terminates
 * algorithm after reaching fitness function \em plateau. The algorithm is
 * terminated if after `n` last generations fitness function maximum has not
 * improved \em relatively more than `frac` with respect to the whole
 * evolutionary process.
 *
 * @tparam G Some `genotype` specialization.
 * @param ff Database intermediary object.
 * @param n Number of \em last generations.
 * @param frac \em Plateau \em flatness.
 * @returns Predicate terminating genetic algorithm after reaching fitness
 * function \em plateau.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from at least one generation.
 * @note This condition is not intended for use with `evolution` argument
 * `max_history` different than `0`.
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
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

/**
 * `max_fitness_improvement_termination_2` returns condition, which terminates
 * algorithm after reaching fitness function \em plateau. The algorithm is
 * terminated if after `n` last generations fitness function maximum has not
 * improved \em absolutely more than `delta` with respect to the whole
 * evolutionary process.
 *
 * @tparam G Some `genotype` specialization.
 * @param ff Database intermediary object.
 * @param n Number of \em last generations.
 * @param delta \em Plateau \em flatness.
 * @returns Predicate terminating genetic algorithm after reaching fitness
 * function \em plateau.
 *
 * @throws std::runtime_error Exception is raised if fitness function evaluates
 * to `incalculable` for all genotypes from at least one generation.
 * @note This condition is not intended for use with `evolution` argument
 * `max_history` different than `0`.
 */
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

/**
 * `threshold_termination` returns condition, which terminates algorithm if at
 * least one genotype fulfills predicate `thr`.
 *
 * @tparam G Some `genotype` specialization.
 * @tparam F Predicate type.
 * @param thr Predicate identifying searched genotype.
 * @returns Predicate terminating genetic algorithm after genotype satisfying
 * `thr` predicate is found.
 */
template<typename G, typename F>
requires chromosome<G> && std::predicate<F, G> termination_condition_fn<G>
threshold_termination(const F& thr)
{
  return [=](std::size_t, const generations<G>& gs) {
    return gs.empty() ? false : std::ranges::any_of(gs.back(), thr);
  };
}

/**
 * `fitness_threshold_termination` returns condition, which terminates algorithm
 * if at least one genotype reaches `thr` fitness function value with absolute
 * precision `eps`.
 *
 * @tparam G Some `genotype` specialization.
 * @param fd Database intermediary object.
 * @param thr Fitness function value to achieve.
 * @param eps Fitness function value absolute precision.
 * @returns Predicate terminating genetic algorithm after genotype reaching
 * fitness function value `thr` with absolute precision `eps` is found.
 */
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

/**
 * `Gaussian_mutation` returns Gaussian mutation operator with standard
 * deviation `sigma` and gene mutation probability `p`.
 *
 * @tparam G Some `genotype` specialization.
 * @param sigma Standard deviation.
 * @param p Gene mutation probability.
 * @returns Gaussian mutation operator.
 *
 * Example:
 * @include Gaussian_mutation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude Gaussian_mutation.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
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

/**
 * `self_adaptive_mutation` returns self adaptive mutation operator with
 * parameters `a0` and `a1`.
 *
 * @tparam G Some `genotype` specialization.
 * @param a0 Self adaptive mutation parameter.
 * @param a1 Self adaptive mutation parameter.
 * @returns Self adaptive mutation operator.
 *
 * @note Due to documentation processing problem the above template declaration
 * is incorrect. Corrected declaration:
 * @code
 * template<typename G>
 * requires floating_point_chromosome<G> && (G::size() % 2 == 0)
 * auto self_adaptive_mutation(typename G::gene_t a0, typename G::gene_t a1)
 * @endcode
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
template<typename G>
requires floating_point_chromosome<G>
  /**
   * \cond
   */
  &&(G::size() % 2 == 0)
  // This unfortunately cannot be processed properly by documentation system.
  /**
   * \endcond
   */
  auto self_adaptive_mutation(typename G::gene_t a0, typename G::gene_t a1)
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

/**
 * `swap_mutation` is swap mutation.
 *
 * @tparam G Some `genotype` specialization.
 * @param g Genotype.
 * @returns Population containing mutated genotype.
 *
 * Example:
 * @include variation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude variation.out
 */
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

/**
 * `random_reset` returns random reset mutation with parameter `p`.
 *
 * @tparam G Some `genotype` specialization.
 * @param p Gene mutation probability.
 * @returns Random reset mutation operator.
 */
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

/**
 * `bit_flipping` returns bit-flipping mutation with parameter `p`.
 *
 * @tparam G Some `genotype` specialization.
 * @param p Gene mutation probability.
 * @returns Bit-flipping mutation operator.
 */
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

/**
 * `arithmetic_recombination` is arithmetic recombination.
 *
 * @tparam Some `G` specialization.
 * @param g0 First parent.
 * @param g1 Secong parent.
 * @returns Population containing one offspring genotype.
 *
 * Example:
 * @include self_adaptive.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude self_adaptive.out
 *
 * @note Evolution result is saved in separate file (not included).
 */
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

/**
 * `single_arithmetic_recombination` is single arithmetic recombination.
 *
 * @tparam Some `G` specialization.
 * @param g0 First parent.
 * @param g1 Secong parent.
 * @returns Population containing two offspring genotypes.
 *
 * Example:
 * @include recombination.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude recombination.out
 */
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

/**
 * `one_point_xover` is one-point crossover recombination.
 *
 * @tparam Some `G` specialization.
 * @param g0 First parent.
 * @param g1 Secong parent.
 * @returns Population containing two offspring genotypes.
 *
 * Example:
 * @include recombination.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude recombination.out
 */
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

/**
 * `cut_n_crossfill` is cut-and-crossfill recombination.
 *
 * @tparam Some `G` specialization.
 * @param g0 First parent.
 * @param g1 Secong parent.
 * @returns Population containing two offspring genotypes.
 *
 * Example:
 * @include variation.cc
 *
 * Result (might be different due to randomness):
 * @verbinclude variation.out
 */
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

/**
 * `test_functions` contains mechanisms for floating-point test functions.
 */
namespace test_functions {

/**
 * `test_functions::point` is point in N-dimensional space.
 *
 * @tparam T Floating-point type.
 * @tparam N Space dimension.
 */
template<std::floating_point T, std::size_t N>
using point = std::array<T, N>;

/**
 * `test_functions::distance` returns distance between two points.
 *
 * @tparam T Floating-point type.
 * @tparam N Space dimension.
 * @param p0 First point.
 * @param p1 Second point.
 * @returns Distance between `p0` and `p1`.
 */
template<std::floating_point T, std::size_t N>
T
distance(const point<T, N>& p0, const point<T, N>& p1)
{
  T res = .0;
  for (std::size_t i = 0; i < N; ++i) {
    res += square(p0[i] - p1[i]);
  }
  return std::sqrt(res);
}

/**
 * `test_functions::coordinates` converts 2D-point to `std::tuple` containing
 * point coordinates.
 *
 * @tparam T Floating-point type.
 * @param p 2D-point.
 * @returns Corresponding tuple.
 */
template<std::floating_point T>
std::tuple<T, T>
coordinates(const point<T, 2>& p)
{
  return std::tuple<T, T>{ p[0], p[1] };
}

/**
 * `test_functions::coordinates` converts 3D-point to `std::tuple` containing
 * point coordinates.
 *
 * @tparam T Floating-point type.
 * @param p 3D-point.
 * @returns Corresponding tuple.
 */
template<std::floating_point T>
std::tuple<T, T, T>
coordinates(const point<T, 3>& p)
{
  return std::tuple<T, T, T>{ p[0], p[1], p[2] };
}

/**
 * `test_functions::uniform_point` returns point in N-dimensional space, where
 * each coordinate has the same value `v`.
 *
 * @tparam T Floating-point type.
 * @tparam Space dimension.
 * @param v Coordinate value.
 * @returns Point in N-dimensional space.
 */
template<std::floating_point T, std::size_t N>
point<T, N>
uniform_point(T v)
{
  point<T, N> res{};
  std::ranges::generate(res, [=]() -> T { return v; });
  return res;
}

/**
 * `test_functions::test_function` is floating-point test function.
 *
 * @tparam T Floating-point type.
 * @tparam N Space dimension.
 */
template<std::floating_point T, std::size_t N>
class test_function
{
public:
  /**
   * `test_functions::test_function::function` is underlying test function
   * type.
   */
  using function = std::function<T(const point<T, N>&)>;

  /**
   * `test_functions::test_function::domain_fn` is test function domain.
   */
  using domain_fn = std::function<domain<T, N>()>;

  /**
   * `test_functions::test_function::point_fn` is solution generating function
   * type.
   */
  using point_fn = std::function<point<T, N>()>;

public:
  /**
   * `test_functions::test_function::test_function` creates test function.
   *
   * @param name Test function name.
   * @param fn The test function representation.
   * @param d The test function domain.
   * @param p_min Function generating solution minimizing the test function.
   */
  test_function(const std::string& name,
                const function& fn,
                const domain_fn& d,
                const point_fn& p_min)
    : name_{ name }
    , fn_{ fn }
    , d_{ d }
    , p_min_{ p_min }
  {
  }

  /**
   * `test_functions::test_function::name` returns test function name.
   *
   * @returns Name of the test function.
   */
  std::string name() const { return name_; }

  /**
   * `test_functions::test_function::operator()` returns test function value
   * at point `p`.
   *
   * @param p Point.
   * @returns Test function value at point `p`.
   */
  T operator()(const point<T, N>& p) const { return fn_(p); }

  /**
   * `test_functions::test_function::function_domain` returns test function
   * domain.
   *
   * @returns Test function domain.
   */
  domain<T, N> function_domain() const { return d_(); }

  /**
   * `test_functions::test_function::p_min` returns point minimizing test
   * function over its domain.
   *
   * @returns Point minimizing test function over its domain.
   */
  point<T, N> p_min() const { return p_min_(); }

private:
  std::string name_;
  function fn_;
  domain_fn d_;
  point_fn p_min_;
};

/**
 * `test_functions::Ackley` is Ackley test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * -20 \exp\left( \frac{-0{.}02}{\sqrt{n}} \sqrt{\sum_{i = 0}^{n - 1} x_i^2}
 * \right) - \exp\left( \frac{1}{n} \sum_{i = 0}^{n - 1} \cos \left( 2 \pi x_i
 * \right) \right) + 20 + e
 * \f]
 */
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

/**
 * `test_functions::Alpine` is Alpine test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \sum_{i = 0}^{n - 1} \left| x_i \sin x_i  + 0{.}1 x_i \right|
 * \f]
 */
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

/**
 * `test_functions::Aluffi_Pentini` is Aluffi-Pentini test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \frac{1}{4} x_0^4 - \frac{1}{2} x_0^2 + \frac{1}{10} x_0 + \frac{1}{2}
 * x_1^2
 * \f]
 */
template<std::floating_point T>
const test_function<T, 2> Aluffi_Pentini{
  "Aluffi-Pentini",
  [](const point<T, 2>& p) {
    const auto [x, y] = coordinates(p);
    return ((.25 * x * x - .5) * x + .1) * x + 0.5 * y * y;
  },
  []() { return uniform_domain<T, 2>(-10., 10.); },
  []() {
    return point<T, 2>{ [q = 0.1](int k) -> T {
                         return 2. * std::sqrt(3.) *
                                std::cos(
                                  std::acos(-3. * std::sqrt(3.) * q / 2.) / 3. -
                                  2. * std::numbers::pi_v<T> * k / 3.) /
                                3.;
                       }(2),
                        0. };
  }
};

/**
 * `test_functions::Booth` is Booth test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * (x_0 + 2x_1 - 7)^2 + (2x_0 + x_1 - 5)^2
 * \f]
 */
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

/**
 * `test_functions::Colville` is Colville test function.
 *
 * \f{eqnarray*}{
 * f^*\left(\vec{x}\right) & = &
 * 100 \left( x_0 - x_1^2 \right)^2 + \left( 1 - x_0 \right)^2 + 90 \left( x_3
 * - x_2^2 \right)^2 + \left( 1 - x_2 \right)^2 \\
 * & & +\ 10{.}1 \left( x_1 - 1
 * \right)^2 + \left( x_3 - 1 \right)^2 + 19{.}8 \left( x_1 - 1 \right) \left(
 * x_3 - 1 \right)
 * \f}
 */
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

/**
 * `test_functions::Easom` is Easom test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * -\cos x_0 \cdot \cos x_1 \cdot \exp \left( -(x_0 - \pi )^2 - (x_1 - \pi )^2
 * \right)
 * \f]
 */
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

/**
 * `test_functions::exponential` is exponential test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * -\exp\left( -\frac{1}{2} \sum_{i = 0}^{n - 1} x_i^2 \right)
 * \f]
 */
template<std::floating_point T, std::size_t N>
const test_function<T, N> exponential{
  "exponential",
  [](const point<T, N>& p) {
    return -std::exp(
      -.5 * std::transform_reduce(
              std::begin(p), std::end(p), T{ 0. }, std::plus<T>{}, square<T>));
  },
  []() { return uniform_domain<T, N>(-1., 1.); },
  []() { return uniform_point<T, N>(0.); }
};

/**
 * `test_functions::Goldstein_Price` is Goldstein-Price test function.
 *
 * \f{eqnarray*}{
 * f^*\left(\vec{x}\right) & = &
 * \left( 1 + \left( x_0 + x_1 + 1 \right)^2 \left( 19 - 14 x_0 + 3 x_0^2
 * - 14
 * x_1 + 6 x_0 x_1 + 3 x_1^2 \right) \right) \\
 * & & \cdot \left( 30 + \left( 2 x_0 -
 * 3 x_1 \right)^2 \left( 18 - 32 x_0 + 12 x_0^2 + 48 x_1 - 36 x_0 x_1 + 27
 * x_1^2 \right) \right)
 * \f}
 */
template<std::floating_point T>
const test_function<T, 2> Goldstein_Price{
  "Goldstein-Price",
  [](const point<T, 2>& p) {
    const auto [x, y] = coordinates(p);
    const auto [x2, y2] = std::tuple<T, T>{ x * x, y * y };
    const auto xy = x * y;
    return (1. + square(x + y + 1.) *
                   (19. - 14. * x + 3. * x2 - 14. * y + 6. * xy + 3. * y2)) *
           (30. + square(2. * x - 3. * y) *
                    (18. - 32. * x + 12. * x2 + 48. * y - 36. * xy + 27. * y2));
  },
  []() { return uniform_domain<T, 2>(-2., 2.); },
  []() {
    return point<T, 2>{ 0., -1. };
  }
};

/**
 * `test_functions::Hosaki` is Hosaki test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \left( 1 - 8 x_0 + 7 x_0^2 - \frac{7}{3} x_0^3 + \frac{1}{4} x_0^4 \right)
 * x_1^2 \exp (-x_1)
 * \f]
 */
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

/**
 * `test_functions::Leon` is Leon test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * 100 \left( x_1 - x_0^2 \right)^2 + \left( 1 - x_0 \right)^2
 * \f]
 */
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

/**
 * `test_functions::Matyas` is Matyas test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * 0{.}26 \left( x_0^2 + x_1^2 \right) - 0{.}48 x_0 x_1
 * \f]
 */
template<std::floating_point T>
const test_function<T, 2> Matyas{ "Matyas",
                                  [](const point<T, 2>& p) {
                                    const auto [x, y] = coordinates(p);
                                    return .26 * (x * x + y * y) - .48 * x * y;
                                  },
                                  []() {
                                    return uniform_domain<T, 2>(-10., 10.);
                                  },
                                  []() { return uniform_point<T, 2>(0.); } };

/**
 * `test_functions::Mexican_hat` is Mexican hat test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * -20 \frac{\sin g(x_0, x_1)}{g(x_0, x_1)}, \, g(x_0, x_1) = 0{.}1 +
 * \sqrt{(x_0 - 4)^2 + (x_1 - 4)^2}
 * \f]
 */
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

/**
 * `test_functions::Miele_Cantrell` is Miele-Cantrell test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \left( \exp (-x_0) - x_1 \right)^4 + 100 \left( x_1 - x_2 \right)^6 +
 * \tan^4 (x_2 - x_3) + x_0^8
 * \f]
 */
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

/**
 * `test_functions::Rosenbrock` is Rosenbrock test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \sum_{i = 0}^{n - 2} \left( 100 \left( x_{i + 1} - x_i^2 \right)^2 + \left(
 * x_i - 1 \right)^2 \right)
 * \f]
 */
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

/**
 * `test_functions::Schwefel` is Schwefel test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \sum_{i = 0}^{n - 1} \left( \sum_{j = 0}^i x_i \right)^2
 * \f]
 */
template<std::floating_point T, std::size_t N>
const test_function<T, N> Schwefel{ "Schwefel",
                                    [](const point<T, N>& p) {
                                      T res = 0.;
                                      for (T sum = 0.; auto x : p) {
                                        res += square(sum += x);
                                      }
                                      return res;
                                    },
                                    []() {
                                      return uniform_domain<T, N>(-100., 100.);
                                    },
                                    []() { return uniform_point<T, N>(0.); } };

/**
 * `test_functions::sphere` is sphere test function.
 *
 * \f[
 * f^*\left(\vec{x}\right) =
 * \sum_{i = 0}^{n - 1} x_i^2
 * \f]
 */
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
