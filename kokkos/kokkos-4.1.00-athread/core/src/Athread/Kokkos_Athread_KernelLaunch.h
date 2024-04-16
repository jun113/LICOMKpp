#ifndef KOKKOS_ATHREAD_KERNELLAUNCH_H_
#define KOKKOS_ATHREAD_KERNELLAUNCH_H_

#include "Kokkos_Athread_ParamWrap.h"

extern "C" void launch_register_kernel();

extern "C" void athread_parallel_for_launch_1D (struct AthreadParamWrap*);
extern "C" void athread_parallel_for_launch_2D (struct AthreadParamWrap*);
extern "C" void athread_parallel_for_launch_3D (struct AthreadParamWrap*);
extern "C" void athread_parallel_for_launch_4D (struct AthreadParamWrap*);

extern "C" void athread_parallel_reduce_launch_1D (struct AthreadParamWrap*);
extern "C" void athread_parallel_reduce_launch_2D (struct AthreadParamWrap*);
extern "C" void athread_parallel_reduce_launch_3D (struct AthreadParamWrap*);

#endif // KOKKOS_ATHREAD_KERNELLAUNCH_H_