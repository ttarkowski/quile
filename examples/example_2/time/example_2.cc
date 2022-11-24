#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <quile/quile.h>

using namespace quile;

using type = int;
const std::size_t n = NUMBER;

fitness f(const auto& chessboard)
{
  assert(n == chessboard.size());
  int counts = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const int d = j - i;
      if (chessboard[i] == chessboard[j] ||
          chessboard[i] + d == chessboard[j] ||
          chessboard[i] == chessboard[j] + d) {
        ++counts;
      }
    }
  }
  return -counts;
}

using G = genotype<g_permutation<type, n, 0>>;

std::string Forsyth_Edwards_Notation(const G& g)
{
  std::string res{};
  for (std::size_t i = 0; auto x : g) {
    if (x != 0) {
      res += std::to_string(x);
    }
    res += 'Q';
    if (x != (n - 1)) {
      res += std::to_string(n - 1 - x);
    }
    if (++i != n) {
      res += '/';
    }
  }
  res += " w - - 0 0";
  return res;
}

int main()
{
  const fitness_function<G> ff = [](const G& g) { return f(g.data()); };
  const fitness_db<G> fd{ ff, constraints_satisfied<G> };
  const fitness_proportional_selection<G> fps{ fd };

  const auto p0 = random_population<constraints_satisfied<G>, G>;
  const auto p1 = stochastic_universal_sampling<G>{ fps };
  const auto p2 = adapter<G>(stochastic_universal_sampling<G>{ fps });

  const std::size_t generation_sz{ 20 };
  const std::size_t parents_sz{ 2 };
  const auto tc = fitness_threshold_termination<G>(fd, 0., 0.01);

  const variation<G> v{ swap_mutation<G>, cut_n_crossfill<G> };

  const auto e = evolution<G>(v, p0, p1, p2, tc, generation_sz, parents_sz, 1);
  std::ofstream file{ "solution.dat" };
  file << Forsyth_Edwards_Notation(fd.rank_order()[0]) << '\n';
  std::cout << fd.size() << " genotypes created in " << e.size()
            << " generations.\n";
}
