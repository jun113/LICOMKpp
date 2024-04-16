#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_

#include "def-undef.h"
#include <cmath>

namespace CppConstantMod {

//constexpr double PI = 4.0 * atan(1.0);
#define PI (double)(4.0 * atan(1.0))
constexpr double G      = 9.806;

constexpr double C0     = 0.0;
constexpr double C1     = 1.0;
constexpr double C2     = 2.0;
constexpr double C3     = 3.0;
constexpr double C4     = 4.0;
constexpr double C5     = 5.0;
constexpr double C8     = 8.0;
constexpr double C10    = 10.0;
constexpr double C16    = 16.0;
constexpr double C1000  = 1000.0;
constexpr double C10000 = 10000.0;
constexpr double C1P5   = 1.5;
constexpr double P33    = C1 / C3;
constexpr double P5     = 0.5;
constexpr double P25    = 0.25;
constexpr double P125   = 0.125;
constexpr double P001   = 0.001;

} // namespace CppConstantMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_
