#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <cstdint>


union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};

#endif
