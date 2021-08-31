#ifndef COMMON_PWX_H
#define COMMON_PWX_H

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

struct pwx_atom
{
  std::string symbol;
  double mass;
  std::string pp;
};

struct pwx_position
{
  std::string symbol;
  double x;
  double y;
  double z;

  bool operator==(const pwx_position&) const = default;

  auto different() const
  {
    return [q = *this](const pwx_position& p) { return p != q; };
  }

  double distance(const pwx_position& p) const
  {
    return std::hypot(x - p.x, y - p.y, z - p.z);
  }
};

inline double
pwx_distance(const std::tuple<pwx_position, pwx_position>& t)
{
  return std::get<0>(t).distance(std::get<1>(t));
}

double
pwx_distance_pbc(const std::tuple<pwx_position, pwx_position>& t,
		 double x, double y, double z);

using pwx_atoms = std::vector<pwx_atom>;

using pwx_positions = std::vector<pwx_position>;

std::string
pwx_fixed(double d);

std::string
pwx_scientific(double d);

std::string
pwx_control(const std::string& calculation, const std::string& outdir_postfix);

std::string
pwx_system(int ibrav,
           int nat,
           int ntyp,
           double tot_charge,
           double degauss,
           double ecutwfc);

std::string
pwx_electrons(int electron_maxstep, double mixing_beta);

std::string
pwx_cell_parameters_triangle_60deg(double a,
                                   std::size_t m,
                                   std::size_t n,
                                   double z);

std::string
pwx_cell_parameters_diag(double x, double y, double z);

std::string
pwx_atomic_species(const pwx_atoms& as);

std::string
pwx_atomic_positions(const pwx_positions& ps);

std::string
pwx_k_points(int nk1, int nk2, int nk3, int sk1, int sk2, int sk3);

std::string
pwx_unique_filename();

#endif // COMMON_PWX_H
