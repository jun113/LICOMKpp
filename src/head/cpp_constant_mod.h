#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_

#include "def-undef.h"
#include <cmath>

namespace CppConstantMod {

// constexpr double PI = 4.0 * std::atan(1.0);
// #define PI (double)(4.0 * atan(1.0))
const double PI = 4.0 * std::atan(1.0);
const double DEGTORAD   = PI / 180.0;
constexpr double OMEGA      = 0.7292e-4;
constexpr double G          = 9.806;
                           
constexpr double C0         = 0.0;
constexpr double C1         = 1.0;
constexpr double C2         = 2.0;
constexpr double C3         = 3.0;
constexpr double C4         = 4.0;
constexpr double C5         = 5.0;
constexpr double C8         = 8.0;
constexpr double C10        = 10.0;
constexpr double C16        = 16.0;
constexpr double C1000      = 1000.0;
constexpr double C10000     = 10000.0;
constexpr double C1P5       = 1.5;
constexpr double P33        = C1 / C3;
constexpr double P5         = 0.5;
constexpr double P25        = 0.25;
constexpr double P125       = 0.125;
constexpr double P001       = 0.001;
constexpr double VERY_SMALL = 1.0e-15;
constexpr double KARMAN     = 0.4;

#ifdef CANUTO2010
constexpr double KO    = 1.66;
constexpr double SIG_T = 0.72;
const double PI0_1 = 1.0 / (std::sqrt(27.0 / 5.0 * 
    static_cast<double>(std::pow(static_cast<long double>(KO), 3))) 
        * (1.0 + 1.0 / SIG_T));
constexpr double PI0_2 = 1.0 / 3.0;
#endif

} // namespace CppConstantMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_
