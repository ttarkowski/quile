#include "pwx.h"
#include <atomic>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

double
pwx_distance_pbc(const std::tuple<pwx_position, pwx_position>& t,
                 double x, double y, double z)
{
  pwx_position p{ std::get<1>(t) };
  p.x += x;
  p.y += y;
  p.z += z;
  return std::get<0>(t).distance(p);
}


std::string
pwx_fixed(double d)
{
  std::ostringstream oss{};
  oss << std::fixed << std::setprecision(9) << d;
  return oss.str();
}

std::string
pwx_scientific(double d)
{
  std::ostringstream oss{};
  oss << std::scientific << std::setprecision(9) << d;
  std::string s{ oss.str() };
  return s.replace(s.find('e'), 1, "D");
}

std::string
pwx_control(const std::string& calculation, const std::string& outdir_postfix)
{
  std::ostringstream oss{};
  oss << "&CONTROL\n"
      << "calculation = '" << calculation << "'\n"
      << "prefix = 'dft'\n"
      << "pseudo_dir = './'\n"
      << "outdir = 'results_" << outdir_postfix << "'\n"
      << "/\n";
  return oss.str();
}

std::string
pwx_system(int ibrav,
           int nat,
           int ntyp,
           double tot_charge,
           double degauss,
           double ecutwfc)
{
  std::ostringstream oss{};
  oss << "&SYSTEM\n"
      << "ibrav = " << ibrav << "\n"
      << "nat = " << nat << "\n"
      << "ntyp = " << ntyp << "\n"
      << "tot_charge = " << pwx_scientific(tot_charge) << "\n"
      << "occupations = 'smearing'\n"
      << "smearing = 'methfessel-paxton'\n"
      << "degauss = " << pwx_scientific(degauss) << "\n"
      << "ecutwfc = " << pwx_scientific(ecutwfc) << "\n"
      << "/\n";
  return oss.str();
}

std::string
pwx_electrons(int electron_maxstep, double mixing_beta)
{
  std::ostringstream oss{};
  oss << "&ELECTRONS\n"
      << "electron_maxstep = " << electron_maxstep << "\n"
      << "mixing_beta = " << pwx_scientific(mixing_beta) << "\n"
      << "/\n";
  return oss.str();
}

std::string
pwx_cell_parameters_triangle_60deg(double a,
                                   std::size_t m,
                                   std::size_t n,
                                   double z)
{
  std::ostringstream oss{};
  oss << "CELL_PARAMETERS angstrom\n"
      << pwx_fixed(a * m) << " " << pwx_fixed(0) << " " << pwx_fixed(0) << "\n"
      << pwx_fixed(a * n / 2.) << " " << pwx_fixed(a * n * std::sqrt(3) / 2.)
      << " " << pwx_fixed(0) << "\n"
      << pwx_fixed(0) << " " << pwx_fixed(0) << " " << pwx_fixed(z) << "\n"
      << "\n";
  return oss.str();
}

std::string
pwx_cell_parameters_diag(double x, double y, double z)
{
  std::ostringstream oss{};
  oss << "CELL_PARAMETERS angstrom\n"
      << pwx_fixed(x) << " " << pwx_fixed(0) << " " << pwx_fixed(0) << "\n"
      << pwx_fixed(0) << " " << pwx_fixed(y) << " " << pwx_fixed(0) << "\n"
      << pwx_fixed(0) << " " << pwx_fixed(0) << " " << pwx_fixed(z) << "\n"
      << "\n";
  return oss.str();
}

std::string
pwx_atomic_species(const pwx_atoms& as)
{
  std::ostringstream oss{};
  oss << "ATOMIC_SPECIES\n";
  for (const auto& a : as) {
    oss << a.symbol << " " << pwx_fixed(a.mass) << " " << a.pp << "\n";
  }
  oss << "\n";
  return oss.str();
}

std::string
pwx_atomic_positions(const pwx_positions& ps)
{
  std::ostringstream oss{};
  oss << "ATOMIC_POSITIONS angstrom\n";
  for (auto p : ps) {
    oss << p.symbol << " " << pwx_fixed(p.x) << " " << pwx_fixed(p.y) << " "
        << pwx_fixed(p.z) << "\n";
  }
  oss << "\n";
  return oss.str();
}

std::string
pwx_k_points(int nk1, int nk2, int nk3, int sk1, int sk2, int sk3)
{
  std::ostringstream oss{};
  oss << "K_POINTS automatic\n"
      << nk1 << ' ' << nk2 << ' ' << nk3 << ' ' << sk1 << ' ' << sk2 << ' '
      << sk3 << '\n';
  return oss.str();
}

std::string
pwx_unique_filename()
{
  static std::atomic_size_t i{ 0 };
  std::ostringstream oss{};
  oss << "calc_" << i++ << ".in";
  return oss.str();
}
