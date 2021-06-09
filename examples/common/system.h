#ifndef COMMON_SYSTEM_H
#define COMMON_SYSTEM_H

#include <string>
#include <tuple>

std::tuple<std::string, std::string>
execute(const std::string& command);

#endif // COMMON_SYSTEM_H
