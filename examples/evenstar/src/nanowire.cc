#include "nanowire.h"
#include "../../common/pwx.h"
#include <algorithm>
#include <iterator>

pwx_positions
evenstar::adjust_positions(const pwx_positions& ps)
{
  pwx_positions res{};
  const auto min_x = std::ranges::min_element(ps, {}, &pwx_position::x)->x;
  const auto min_y = std::ranges::min_element(ps, {}, &pwx_position::y)->y;
  std::ranges::transform(
    ps, std::back_inserter(res), [min_x, min_y](const pwx_position& p) {
      return pwx_position{ p.symbol, p.x - min_x, p.y - min_y, p.z };
    });
  return res;
}
