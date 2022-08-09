#ifndef MITHRIL_SRC_NANOTUBE_H
#define MITHRIL_SRC_NANOTUBE_H

#include <algorithm>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>
#include <cstddef>
#include <iterator>
#include <quile/quile.h>
#include <ranges>
#include <sstream>
#include <tuple>
#include <vector>

namespace mithril {

// Finding neighbors in hexagonal lattice with periodic boundary conditions.
class hex_lattice_pbc
{
public:
  hex_lattice_pbc(std::size_t n_phi, std::size_t n_z);
  std::size_t right(std::size_t i) const;
  std::size_t left(std::size_t i) const;
  std::size_t up_right(std::size_t i) const;
  std::size_t up_left(std::size_t i) const { return left(up_right(i)); }
  std::size_t down_right(std::size_t i) const;
  std::size_t down_left(std::size_t i) const { return left(down_right(i)); }
  std::vector<std::size_t> neighbors(std::size_t i) const;

private:
  std::size_t n_phi_;
  std::size_t n_z_;
};

// Finding neighbors in hexagonal lattice without periodic boundary conditions.
class hex_lattice_ord
{
public:
  hex_lattice_ord(std::size_t n_phi, std::size_t n_z);
  std::size_t right(std::size_t i) const;
  std::size_t left(std::size_t i) const;
  std::size_t up_right(std::size_t i) const;
  std::size_t up_left(std::size_t i) const;
  std::size_t down_right(std::size_t i) const;
  std::size_t down_left(std::size_t i) const;
  std::vector<std::size_t> neighbors(std::size_t i) const;

private:
  std::size_t n_phi_;
  std::size_t n_z_;
};

// Atoms' indices.
template<quile::binary_chromosome G>
std::vector<std::size_t>
atoms(const G& g)
{
  std::vector<std::size_t> res{};
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (g.value(i)) {
      res.push_back(i);
    }
  }
  return res;
}

// Number of atoms encoded by genotype.
template<quile::binary_chromosome G>
std::size_t
number_of_atoms(const G& g)
{
  return std::ranges::count(g, true);
}

// Neighbor atoms of index i (periodic boundary condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::vector<std::size_t>
neighbor_atoms_pbc(const G& g, std::size_t i)
{
  assert(i < G::size());
  static const hex_lattice_pbc hl{ n_phi, n_z };
  std::vector<std::size_t> res{};
  std::ranges::copy_if(hl.neighbors(i),
                       std::back_inserter(res),
                       [&](std::size_t j) { return g.value(j); });
  return res;
}

// Neighbor atoms of index i (without periodic boundary condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::vector<std::size_t>
neighbor_atoms_ord(const G& g, std::size_t i)
{
  assert(i < G::size());
  static const hex_lattice_ord hl{ n_phi, n_z };
  std::vector<std::size_t> res{};
  std::ranges::copy_if(hl.neighbors(i),
                       std::back_inserter(res),
                       [&](std::size_t j) { return g.value(j); });
  return res;
}

// Number of neighbor atoms of index i (periodic boundary condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::size_t
number_of_neighbor_atoms_pbc(const G& g, std::size_t i)
{
  assert(i < G::size());
  static const hex_lattice_pbc hl{ n_phi, n_z };
  return std::ranges::count(
    hl.neighbors(i), true, [&](std::size_t j) { return g.value(j); });
}

// Number of neighbor atoms of index i (without periodic boundary condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::size_t
number_of_neighbor_atoms_ord(const G& g, std::size_t i)
{
  assert(i < G::size());
  static const hex_lattice_ord hl{ n_phi, n_z };
  return std::ranges::count(
    hl.neighbors(i), true, [&](std::size_t j) { return g.value(j); });
}

// Motif decomposition (n_0, n_1, ..., n_6).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::vector<std::size_t>
decomposition(const G& g)
{
  std::vector<std::size_t> res{ 0, 0, 0, 0, 0, 0, 0 };
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (g.value(i)) {
      ++res[number_of_neighbor_atoms_pbc<G, n_phi, n_z>(g, i)];
    }
  }
  return res;
}

// Energy from decomposition model.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
double
energy_from_model(const G& g, const double* decomposition_values)
{
  double res{ 0. };
  for (std::size_t i = 0; auto n : decomposition<G, n_phi, n_z>(g)) {
    res += n * decomposition_values[i++];
  }
  return res / number_of_atoms(g);
}

// Predicate testing whether atoms are connected within unit cell.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
bool
atoms_connected_in_unit_cell(const G& g)
{
  static const hex_lattice_ord hl{ n_phi, n_z };
  boost::adjacency_matrix<boost::undirectedS> am{ G::size() };
  for (auto i : atoms(g)) {
    for (auto j : neighbor_atoms_ord<G, n_phi, n_z>(g, i)) {
      boost::add_edge(i, j, am);
    }
  }
  std::vector<std::size_t> c(boost::num_vertices(am));
  return number_of_atoms(g) + boost::connected_components(am, c.data()) ==
         1 + G::size();
}

// Predicate testing whether at least one atom at unit cell boundary along
// nanotube has at least one neighbor (periodic boundary condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
bool
adjacency_at_unit_cell_boundary_along_nanotube(const G& g)
{
  assert(G::size() == 2 * n_phi * n_z);
  static const hex_lattice_pbc hl{ n_phi, n_z };
  for (std::size_t i = 2 * n_z - 1; i < G::size(); i += 2 * n_z) {
    if (g.value(i) && (g.value(hl.up_right(i)) || g.value(hl.down_right(i)))) {
      return true;
    }
  }
  for (std::size_t i = n_z - 1; i < G::size(); i += 2 * n_z) {
    if (g.value(i) && g.value(hl.right(i))) {
      return true;
    }
  }
  return false;
}

// Predicate testing whether at least one atom at unit cell boundary at
// circumference of nanotube has at least one neighbor (periodic boundary
// condition).
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
bool
adjacency_at_unit_cell_boundary_at_circumference(const G& g)
{
  assert(G::size() == 2 * n_phi * n_z);
  static const hex_lattice_pbc hl{ n_phi, n_z };
  for (std::size_t i = G::size() - n_z; i < G::size(); ++i) {
    if (g.value(i) && (g.value(hl.up_left(i)) || g.value(hl.up_right(i)))) {
      return true;
    }
  }
  return false;
}

// operator>> reads genotype from the stream.
// TODO: Move this to the Quilë library.
template<quile::chromosome G>
std::istream&
operator>>(std::istream& is, G& g)
{
  for (std::size_t i = 0; i < G::size(); ++i) {
    typename G::gene_t v{};
    is >> v;
    g.value(i, v);
  }
  return is;
}

// Rotates nanotube representation by positive value delta.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
G
rotate(const G& g, std::size_t delta)
{
  typename G::chain_t res{};
  for (std::size_t i = 0; i < G::size(); ++i) {
    res[(i + 2 * delta * n_z) % (2 * n_phi * n_z)] = g.value(i);
  }
  return G{ res };
}

// Translates nanotube representation by positive value delta.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
G
translate(const G& g, std::size_t delta)
{
  typename G::chain_t res{};
  for (std::size_t i = 0; i < G::size(); ++i) {
    res[((i + delta) % n_z) + (i / n_z) * n_z] = g.value(i);
  }
  return G{ res };
}

// Three-way comparison with unsigned numbers interpretation of function
// arguments (\sum_i 2^i x_i).
template<quile::binary_chromosome G>
auto
compare_as_unsigned_numbers(const G& g0, const G& g1)
{
  auto c0{ g0.data() };
  auto c1{ g1.data() };
  std::reverse(std::begin(c0), std::end(c0));
  std::reverse(std::begin(c1), std::end(c1));
  return c0 <=> c1;
}

// Finds minimum of the nanotube abstract class representation.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
G
find_min_element_of_abstract_class(const G& g)
{
  G min{ g };
  for (std::size_t d_n_phi = 0; d_n_phi < n_phi; ++d_n_phi) {
    for (std::size_t d_n_z = 0; d_n_z < n_z; ++d_n_z) {
      G h{ rotate<G, n_phi, n_z>(translate<G, n_phi, n_z>(g, d_n_z), d_n_phi) };
      if (compare_as_unsigned_numbers(h, min) < 0) {
        min = h;
      }
    }
  }
  return min;
}

namespace detail {

// Computes divisors of a number.
std::vector<std::size_t>
divisors(std::size_t n);

} // namespace detail

// Computes d_n_phi.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::size_t
d_n_phi(const G& g)
{
  static const auto divs{ detail::divisors(n_phi) };
  std::size_t res{ 1 };
  for (auto d : divs) {
    if (g == rotate<G, n_phi, n_z>(g, n_phi / d) && d > res) {
      res = d;
    }
  }
  return res;
}

// Computes d_n_z.
template<quile::binary_chromosome G, std::size_t n_phi, std::size_t n_z>
std::size_t
d_n_z(const G& g)
{
  static const auto divs{ detail::divisors(n_z) };
  std::size_t res{ 1 };
  for (auto d : divs) {
    if (g == translate<G, n_phi, n_z>(g, n_z / d) && d > res) {
      res = d;
    }
  }
  return res;
}

} // namespace mithril

#endif // MITHRIL_SRC_NANOTUBE_H
