#include "system.h"
#include <boost/asio/io_context.hpp>
#include <boost/process/async.hpp>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>
#include <future>

std::tuple<std::string, std::string>
execute(const std::string& command)
{
  boost::asio::io_context ioc{};
  std::future<std::string> out{};
  std::future<std::string> err{};
  boost::process::child c{ command.c_str(),
                           boost::process::std_in.close(),
                           boost::process::std_out > out,
                           boost::process::std_err > err,
                           ioc };
  ioc.run();
  return std::tuple<std::string, std::string>{ out.get(), err.get() };
}
