#ifndef MITHRIL_SRC_SHEET_H
#define MITHRIL_SRC_SHEET_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <quile/quile.h>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

namespace mithril {

namespace detail {

template<std::integral T>
auto
division_condition(T n)
{
  return [=](T i) { return (n / i) * i == n; };
}

template<std::integral T>
auto
divisors(T n)
{
  return std::views::iota(T{ 1 }, n + 1) |
         std::views::filter(division_condition(n));
}

std::ranges::range auto
range_product(std::ranges::range auto r0, std::ranges::range auto r1)
{
  return r0 | std::views::transform([=](auto p) mutable {
           return r1 | std::views::transform(
                         [=](auto q) mutable { return std::make_pair(p, q); });
         }) |
         std::views::join;
}

auto
sorted_pairs_of_divisors(int m, int n)
{
  std::vector<std::pair<int, int>> v{};
  std::ranges::transform(range_product(divisors(m), divisors(n)),
                         std::back_inserter(v),
                         std::identity{});
  std::ranges::sort(
    v, [](auto a, auto b) { return a.first * a.second < b.first * b.second; });
  return v;
}

template<typename T>
std::pair<T, T>
invert_pair(const std::pair<T, T>& p)
{
  return std::pair<T, T>{ p.second, p.first };
}

} // namespace detail

template<binary_chromosome G>
class intermediate_type
{
private:
  using subcontainer = std::vector<bool>;
  using container = std::vector<subcontainer>;

public:
  using reference = typename subcontainer::reference;
  using const_reference = typename subcontainer::const_reference;
  using size_type = typename subcontainer::size_type;
  using value_type = typename subcontainer::value_type;
  struct m_dir;
  struct n_dir;

public:
  intermediate_type(size_type m, size_type n)
    : c_(n, subcontainer(m, false)) // std::vector {}-ctor is problematic.
  {}

  intermediate_type(const intermediate_type&) = default;
  intermediate_type(intermediate_type&&) = default;
  intermediate_type& operator=(const intermediate_type&) = default;
  intermediate_type& operator=(intermediate_type&&) = default;

  std::pair<size_type, size_type> size() const
  {
    return std::pair{ c_.size() ? c_[0].size() : 0, c_.size() };
  }

  reference at(size_type im, size_type in) { return c_.at(in).at(im); }

  const_reference at(size_type im, size_type in) const
  {
    return c_.at(in).at(im);
  }

  void resize(size_type m, size_type n) { return resize(m, n, false); }

  void resize(size_type m, size_type n, const value_type& v)
  {
    std::ranges::transform(c_, c_.begin(), [=](auto& cc) { cc.resize(m, v); });
    c_.resize(n, subcontainer(m, v));
  }

  static intermediate_type create(const G& g, size_type m, size_type n)
  {
    assert(m * n == G::size());
    intermediate_type<G> res{ m, n };
    for (std::size_t i = 0; auto x : g) {
      res.at(i % m, i / m) = x;
      ++i;
    }
    return res;
  }

  template<typename T>
  requires std::is_same_v<T, m_dir> || std::is_same_v<T, n_dir>
  bool periodic(size_type ds)
  {
    const bool b = std::is_same_v<T, m_dir>;
    const auto [k, l] = b ? this->size() : detail::invert_pair(this->size());
    assert(ds <= k && detail::division_condition(k)(ds));
    bool res = true;
    for (size_type i = 0; i < ds && res; ++i) {
      for (size_type s = i + ds; s < k && res; s += ds) {
        for (size_type j = 0; j < l && res; ++j) {
          res &= b ? (this->at(i, j) == this->at(s, j))
                   : (this->at(j, i) == this->at(j, s));
        }
      }
    }
    return res;
  }

private:
  size_type m_size() const { return n_size() ? c_[0].size() : 0; }
  size_type n_size() const { return c_.size(); }

private:
  container c_;
};

} // namespace mithril

#endif // MITHRIL_SRC_SHEET_H
