#ifndef EVENSTAR_SRC_NANOWIRE_H
#define EVENSTAR_SRC_NANOWIRE_H

#include "../../common/pwx.h"
#include <algorithm>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <quile/quile.h>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

namespace evenstar {

template<typename G>
requires quile::floating_point_chromosome<G> constexpr std::size_t
number_of_atoms(bool flat)
{
  return flat ? (1 + G::size()) / 2 : 1 + G::size() / 3;
}

template<std::size_t N_A>
constexpr std::size_t
genotype_size(bool flat)
{
  return flat ? 2 * N_A - 1 : N_A == 1 ? 1 : 3 * (N_A - 1);
}

namespace detail {

template<std::floating_point T, std::size_t N_A>
quile::domain<T, genotype_size<N_A>(true)>
domain_flat(const quile::range<T>& bond)
{
  const T y_max = (N_A - 1) * bond.max();
  const quile::range<T> y{ -y_max, y_max };
  const quile::range<T> dz{ 0., bond.max() };
  quile::domain<T, genotype_size<N_A>(true)> res{};
  res[0] = dz;
  for (std::size_t i = 1; i < genotype_size<N_A>(true); i += 2) {
    res[i] = y;
    res[i + 1] = dz;
  }
  return res;
}

template<std::floating_point T, std::size_t N_A>
quile::domain<T, genotype_size<N_A>(false)>
domain_buckled(const quile::range<T>& bond)
{
  const quile::range<T> rho{ 0., (N_A - 1) * bond.max() };
  const quile::range<T> phi{ 0.,
                             std::nextafter(2 * std::numbers::pi_v<T>, 0.) };
  const quile::range<T> dz{ 0., bond.max() };
  quile::domain<T, genotype_size<N_A>(false)> res{};
  res[0] = dz;
  if (N_A > 1) {
    res[1] = rho;
    res[2] = dz;
  }
  for (std::size_t i = 3; i < genotype_size<N_A>(false); i += 3) {
    res[i] = rho;
    res[i + 1] = phi;
    res[i + 2] = dz;
  }
  return res;
}

} // namespace detail

template<std::floating_point T, std::size_t N_A, bool Flat>
requires(N_A > 0) quile::domain<T, genotype_size<N_A>(Flat)> construct_domain(
  const quile::range<T>& bond)
{
  if constexpr (Flat) {
    return detail::domain_flat<T, N_A>(bond);
  } else {
    return detail::domain_buckled<T, N_A>(bond);
  }
}

namespace detail {

// Mystic Rose is a complete graph with vertices placed on the points of
// a regular polygon. This function returns edges needed to construct this
// graph, e.g. for (0, 1, 2) it returns ((0, 1), (0, 2), (1, 2)).
template<typename C>
std::vector<std::tuple<typename C::value_type, typename C::value_type>>
mystic_rose_edges(const C& c)
{
  using T = typename C::value_type;
  std::vector<std::tuple<T, T>> res{};
  for (std::size_t i = 0; const auto& x : c) {
    for (const auto& y : c | std::views::drop(++i)) {
      res.push_back(std::make_tuple(x, y));
    }
  }
  return res;
}

} // namespace detail

// This function template checks if all atoms are separated from each other.
template<std::floating_point T>
bool
atoms_not_too_close_pbc(const pwx_positions& ps, T h, T min_distance)
{
  return std::ranges::all_of(
           detail::mystic_rose_edges(ps),
           [=](const auto& t) { return pwx_distance(t) >= min_distance; }) &&
         std::ranges::all_of(detail::mystic_rose_edges(ps),
                             [=](const auto& t) {
                               return pwx_distance_pbc(t, 0., 0., h) >=
                                      min_distance;
                             }) &&
         (h >= min_distance);
}

// This function template checks if all atoms are connected, i.e. form a wire.
template<std::floating_point T>
bool
all_atoms_connected_pbc(pwx_positions ps, T h, T max_distance)
{
  const std::size_t n{ ps.size() };
  for (std::size_t i = 0; i < n; ++i) {
    ps.push_back(pwx_position{ ps[i].symbol, ps[i].x, ps[i].y, ps[i].z + h });
  }
  assert(ps.size() == 2 * n);
  boost::adjacency_matrix<boost::undirectedS> g{ ps.size() };
  std::vector<std::size_t> v(ps.size()); // Do not use {...} here!
  std::iota(v.begin(), v.end(), 0);
  for (auto [i, j] : detail::mystic_rose_edges(v)) {
    if (ps[i].distance(ps[j]) <= max_distance) {
      boost::add_edge(i, j, g);
    }
  }
  std::vector<std::size_t> c(boost::num_vertices(g)); // And here.
  return 1 == boost::connected_components(g, c.data());
}

pwx_positions
adjust_positions(const pwx_positions& ps);

namespace detail {

template<typename G>
requires quile::floating_point_chromosome<G> std::tuple<pwx_positions, double>
geometry_flat(const G& g, const std::string& atom_symbol)
{
  // n = number of atoms in unit cell:
  // a) n > 0: dz_n
  // b) n > 1: dz_n, (y_i, dz_i) for i = 1, ..., n - 1
  // Note: G::size() == 2 * n - 1
  using T = typename G::gene_t;
  std::size_t i = 0;
  const T dz_n = g.value(i++);
  T z = 0.;
  pwx_positions res{ pwx_position{ atom_symbol, 0., 0., z } }; // 0
  while (i < G::size()) {                                      // 1, ..., n - 1
    const auto y = g.value(i++);
    z += g.value(i++);
    res.push_back(pwx_position{ atom_symbol, 0., y, z });
  }
  assert(i == G::size() && res.size() == number_of_atoms<G>(true));
  return std::tuple<pwx_positions, double>{ adjust_positions(res), z + dz_n };
}

template<typename G>
requires quile::floating_point_chromosome<G> std::tuple<pwx_positions, double>
geometry_buckled(const G& g, const std::string& atom_symbol)
{
  // n = number of atoms in unit cell:
  // a) n > 0: dz_n
  // b) n > 1: dz_n, rho_1, dz_1
  // c) n > 2: dz_n, rho_1, dz_1, (rho_i, phi_i, dz_i) for i = 2, ..., n - 1
  // Note: G::size() == 1 && n == 1 || g.size() == 3 * (n - 1) && n > 1
  using T = typename G::gene_t;
  std::size_t i = 0;
  const T dz_n = g.value(i++);
  T z = 0.;
  pwx_positions res{ pwx_position{ atom_symbol, 0., 0., z } }; // 0
  if (number_of_atoms<G>(false) > 1) {                         // 1
    const auto [x, y] = quile::polar2cart(g.value(i++), 0.);
    z += g.value(i++);
    res.push_back(pwx_position{ atom_symbol, x, y, z });
  }
  while (i < G::size()) { // 2, ..., n - 1
    const auto rho = g.value(i++);
    const auto [x, y] = quile::polar2cart(rho, g.value(i++));
    z += g.value(i++);
    res.push_back(pwx_position{ atom_symbol, x, y, z });
  }
  assert(i == g.size() && res.size() == number_of_atoms<G>(false));
  return std::tuple<pwx_positions, double>{ adjust_positions(res), z + dz_n };
}

} // namespace detail

template<typename G>
requires quile::floating_point_chromosome<G> std::tuple<pwx_positions, double>
geometry(const G& g, const std::string& atom_symbol, bool flat)
{
  return flat ? detail::geometry_flat<G>(g, atom_symbol)
              : detail::geometry_buckled<G>(g, atom_symbol);
}

} // namespace evenstar

#endif // EVENSTAR_SRC_NANOWIRE_H
