#ifndef COMMON_MAP_H
#define COMMON_MAP_H

#include <mutex>
#include <unordered_map>

template<typename Key, typename T>
class thread_safe_unordered_map
{
private:
  using map_t = std::unordered_map<Key, T>;
  static std::mutex mtx;

public:
  using iterator = typename map_t::iterator;
  using value_type = typename map_t::value_type;

public:
  thread_safe_unordered_map() = default;

  void insert_or_modify(const Key& key, const T& t)
  {
    const std::lock_guard<std::mutex> lg{ mtx };
    m_[key] = t;
  }

  const T& at(const Key& key) const { return m_.at(key); }

private:
  map_t m_;
};

template<typename Key, typename T>
std::mutex thread_safe_unordered_map<Key, T>::mtx{};

#endif // COMMON_MAP_H
