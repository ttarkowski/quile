#include "nanotube.h"
#include <cstddef>
#include <vector>

mithril::hex_lattice_pbc::hex_lattice_pbc(std::size_t n_phi, std::size_t n_z)
  : n_phi_{ n_phi }
  , n_z_{ n_z }
{}

std::size_t
mithril::hex_lattice_pbc::right(std::size_t i) const
{
  return i + 1 - (i % n_z_ != n_z_ - 1 ? 0 : n_z_);
}

std::size_t
mithril::hex_lattice_pbc::left(std::size_t i) const
{
  return i + (i % n_z_ ? 0 : n_z_) - 1;
}

std::size_t
mithril::hex_lattice_pbc::up_right(std::size_t i) const
{
  return ((i / n_z_) % 2 == 0 ? i + n_z_
                              : i + 1 + (i % n_z_ != n_z_ - 1 ? n_z_ : 0)) %
         (2 * n_phi_ * n_z_);
}

std::size_t
mithril::hex_lattice_pbc::down_right(std::size_t i) const
{
  i += 2 * n_phi_ * n_z_;
  return ((i / n_z_) % 2 == 0 ? i - n_z_
                              : i + 1 - n_z_ * (i % n_z_ != n_z_ - 1 ? 1 : 2)) %
         (2 * n_phi_ * n_z_);
}

std::vector<std::size_t>
mithril::hex_lattice_pbc::neighbors(std::size_t i) const
{
  return std::vector{ up_left(i), up_right(i),  left(i),
                      right(i),   down_left(i), down_right(i) };
}

mithril::hex_lattice_ord::hex_lattice_ord(std::size_t n_phi, std::size_t n_z)
  : n_phi_{ n_phi }
  , n_z_{ n_z }
{}

std::size_t
mithril::hex_lattice_ord::right(std::size_t i) const
{
  return i + (i % n_z_ != n_z_ - 1 ? 1 : 0);
}

std::size_t
mithril::hex_lattice_ord::left(std::size_t i) const
{
  return i - (i % n_z_ ? 1 : 0);
}

std::size_t
mithril::hex_lattice_ord::up_right(std::size_t i) const
{
  const std::size_t aux{ i / n_z_ };
  return i + (aux == 2 * n_phi_ - 1                      ? 0
              : (aux % 2 == 1) && (i % n_z_ == n_z_ - 1) ? 0
              : (aux % 2 == 1) && (i % n_z_ != n_z_ - 1) ? n_z_ + 1
                                                         : n_z_);
}

std::size_t
mithril::hex_lattice_ord::up_left(std::size_t i) const
{
  const std::size_t aux{ i / n_z_ };
  return i + (aux == 2 * n_phi_ - 1               ? 0
              : (aux % 2 == 0) && (i % n_z_ == 0) ? 0
              : (aux % 2 == 0) && (i % n_z_ != 0) ? n_z_ - 1
                                                  : n_z_);
}

std::size_t
mithril::hex_lattice_ord::down_right(std::size_t i) const
{
  const std::size_t aux{ i / n_z_ };
  return i - (aux == 0                                   ? 0
              : (aux % 2 == 1) && (i % n_z_ == n_z_ - 1) ? 0
              : (aux % 2 == 1) && (i % n_z_ != n_z_ - 1) ? n_z_ - 1
                                                         : n_z_);
}

std::size_t
mithril::hex_lattice_ord::down_left(std::size_t i) const
{
  const std::size_t aux{ i / n_z_ };
  return i - (aux == 0                            ? 0
              : (aux % 2 == 0) && (i % n_z_ == 0) ? 0
              : (aux % 2 == 0) && (i % n_z_ != 0) ? n_z_ + 1
                                                  : n_z_);
}

std::vector<std::size_t>
mithril::hex_lattice_ord::neighbors(std::size_t i) const
{
  std::vector<std::size_t> res{};
  for (const auto f : { &hex_lattice_ord::up_left,
                        &hex_lattice_ord::up_right,
                        &hex_lattice_ord::left,
                        &hex_lattice_ord::right,
                        &hex_lattice_ord::down_left,
                        &hex_lattice_ord::down_right }) {
    if (const auto aux = (this->*f)(i); aux != i) {
      res.push_back(aux);
    }
  }
  return res;
}

std::vector<std::size_t>
mithril::detail::divisors(std::size_t n)
{
  std::vector<std::size_t> res{};
  for (std::size_t i = 1; i <= n; ++i) {
    if (n % i == 0) {
      res.push_back(i);
    }
  }
  return res;
}
