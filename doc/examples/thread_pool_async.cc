#include <chrono>
#include <iostream>
#include <quile/quile.h>
#include <thread>
#include <vector>

int
main()
{
  using namespace std::chrono_literals;
  quile::thread_pool tp{ 4 };
  std::vector<std::future<void>> v0{};
  for (int i = 0; i < 10; ++i) {
    v0.push_back(tp.async<void>(std::launch::async, [i]() {
      std::cout << i << '\n';
      std::this_thread::sleep_for((i % 2 + 1) * 1s);
    }));
  }
  for (auto& x : v0) {
    x.get();
  }
}
