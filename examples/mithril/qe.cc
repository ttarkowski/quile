#include "src/nanotube.h"
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <quile/quile.h>

using namespace quile;
using namespace mithril;

namespace {

const std::size_t n_phi = N_PHI;
const std::size_t n_z = N_Z;
const std::size_t c = 2 * n_phi * n_z;
using G = genotype<g_binary<2 * n_phi * n_z>>;

template<binary_chromosome G>
G
read_genotype(std::istream& is)
{
  G g{};
  is >> g;
  return g;
}

const double a = 1.675;              // lattice constant (A)
const double buffer_thickness = 15.; // buffer thickness (A)

double rho_coord(std::size_t)
{
  const double p = a * std::sqrt(6.);
  const double q = 4. * std::sqrt(1. - std::cos(pi<double> / n_phi));
  return p / q;
}

double
phi_coord(std::size_t i)
{
  return pi<double> * (i / n_z) / n_phi;
}

double
z_coord(std::size_t i)
{
  return (i % n_z) * a + (((i / n_z) % 2 == 0) ? 0. : a / 2.);
}

double
rho()
{
  return rho_coord(0);
}

template<binary_chromosome G>
void
qe_input(std::ostream& os, const G& g)
{
  os << "&CONTROL\n"
     << "calculation = 'vc-relax'\n"
     << "nstep = 200\n"
     << "prefix = 'dft'\n"
     << "pseudo_dir = './'\n"
     << "outdir = 'results'\n"
     << "/\n"
     << "&SYSTEM\n"
     << "ibrav = 0\n"
     << "nat = " << number_of_atoms(g) << '\n'
     << "ntyp = 1\n"
     << "tot_charge = 0.000000000D+00\n"
     << "occupations = 'smearing'\n"
     << "smearing = 'methfessel-paxton'\n"
     << "degauss = 1.000000000D-02\n"
     << "ecutwfc = 6.000000000D+01\n"
     << "/\n"
     << "&ELECTRONS\n"
     << "electron_maxstep = 100\n"
     << "mixing_beta = 7.000000000D-01\n"
     << "/\n"
     << "&IONS\n"
     << "/\n"
     << "&CELL\n"
     << "cell_dofree = 'z'\n"
     << "/\n"
     << "CELL_PARAMETERS angstrom\n";
  os << std::fixed << std::setprecision(9) << 2 * rho() + buffer_thickness;
  os << " 0.000000000 0.000000000\n"
     << "0.000000000 ";
  os << std::fixed << std::setprecision(9) << 2 * rho() + buffer_thickness;
  os << " 0.000000000\n"
     << "0.000000000 0.000000000 ";
  os << std::fixed << std::setprecision(9) << n_z * a << '\n';
  os << '\n'
     << "ATOMIC_SPECIES\n"
     << "B 10.811000000 B.pbe-n-kjpaw_psl.1.0.0.UPF\n"
     << '\n'
     << "ATOMIC_POSITIONS angstrom\n";
  for (std::size_t i = 0; i < G::size(); ++i) {
    if (g.value(i)) {
      const auto [x, y] = polar2cart(rho_coord(i), phi_coord(i));
      const auto z = z_coord(i);
      os << "B" << ' ';
      os << std::fixed << std::setprecision(9)
         << x + rho() + buffer_thickness / 2. << ' ';
      os << std::fixed << std::setprecision(9)
         << y + rho() + buffer_thickness / 2. << ' ';
      os << std::fixed << std::setprecision(9) << z << '\n';
    }
  }
  os << '\n'
     << "K_POINTS automatic\n"
     << "1 1 16 0 0 1\n";
}

}

int
main()
{
  qe_input(std::cout, read_genotype<G>(std::cin));
}
