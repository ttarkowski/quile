#include <quile/quile.h>

constexpr auto d0 = quile::uniform_domain<int, 42>(quile::range{ 0, 100 });
constexpr auto d1 = quile::uniform_domain<int, 42>(0, 100);

static_assert(quile::uniform(d0));
static_assert(quile::uniform(d1));

int
main()
{}
