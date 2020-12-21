#include "pwx.h"
#include <atomic>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

std::string
evenstar::pwx_fixed(double d)
{
  std::ostringstream oss{};
  oss << std::fixed << std::setprecision(9) << d;
  return oss.str();
}

std::string
evenstar::pwx_scientific(double d)
{
  std::ostringstream oss{};
  oss << std::scientific << std::setprecision(9) << d;
  std::string s{ oss.str() };
  return s.replace(s.find('e'), 1, "D");
}

std::string
evenstar::pwx_control(const std::string& outdir_postfix)
{
  std::ostringstream oss{};
  oss << "&CONTROL\n"
      << "calculation = 'scf'\n"
      << "prefix = 'dft'\n"
      << "pseudo_dir = './'\n"
      << "outdir = 'results-" << outdir_postfix << "'\n"
      << "/\n";
  return oss.str();
}

std::string
evenstar::pwx_system(int nat,
                     int ntyp,
                     double tot_charge,
                     double degauss,
                     double ecutwfc)
{
  std::ostringstream oss{};
  oss << "&SYSTEM\n"
      << "ibrav = 0\n"
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
evenstar::pwx_electrons(int electron_maxstep, double mixing_beta)
{
  std::ostringstream oss{};
  oss << "&ELECTRONS\n"
      << "electron_maxstep = " << electron_maxstep << "\n"
      << "mixing_beta = " << pwx_scientific(mixing_beta) << "\n"
      << "/\n";
  return oss.str();
}

std::string
evenstar::pwx_cell_parameters_diag(double x, double y, double z)
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
evenstar::pwx_atomic_species(const pwx_atoms& as)
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
evenstar::pwx_atomic_positions(const pwx_positions& ps)
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
evenstar::pwx_k_points(int k)
{
  std::ostringstream oss{};
  oss << "K_POINTS automatic\n"
      << "1 1 " << k << " 0 0 1\n";
  return oss.str();
}

std::string
evenstar::pwx_unique_filename()
{
  static std::atomic_size_t i{ 0 };
  std::ostringstream oss{};
  oss << "stripe-" << i++ << ".in";
  return oss.str();
}
