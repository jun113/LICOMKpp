#ifndef KOKKOS_ATHREAD_MATH_HPP_
#define KOKKOS_ATHREAD_MATH_HPP_

#include "Kokkos_Athread_MathematicalFunctionsSlave.h"

namespace Kokkos {

template<typename T> T fmax(T a, T b) {
	return a > b ? a : b;
}

template<typename T> T fmin(T a, T b) {
	return a < b ? a : b;
}

template<typename T> T fabs (T x) {
  return x < (T)0.0 ? -x : x;
}

template<typename T> T sqrt (T x) {
	return athread_slave_sqrt(x);
}

template<typename T> T log (T x) {
	return athread_slave_log(x);
}

template<typename T> T atan (T x) {
	return athread_slave_atan_double(x);
}
template<> double atan (double x) {
	return athread_slave_atan_double(x);
}
template<typename T> T sin (T x) {
	return athread_slave_sin_double(x);
}

template<> double sin (double x) {
	return athread_slave_sin_double(x);
}

template<typename T> T cos (T x) {
	return athread_slave_cos_double(x);
}
template<> double cos (double x) {
	return athread_slave_cos_double(x);
}

} // namespace Kokkos

#endif // KOKKOS_ATHREAD_MATH_HPP_