#ifndef KOKKOS_ATHREAD_MATHEMATICALFUNCTIONS_H_
#define KOKKOS_ATHREAD_MATHEMATICALFUNCTIONS_H_

extern "C" double athread_slave_sqrt(double);
extern "C" double athread_slave_log(double);

extern "C" double athread_slave_atan_double(double);
extern "C" double athread_slave_sin_double(double);
extern "C" double athread_slave_cos_double(double);

// extern double fmax(double, double);
// extern double fmin(double, double);

#endif // KOKKOS_ATHREAD_MATHEMATICALFUNCTIONS_H_