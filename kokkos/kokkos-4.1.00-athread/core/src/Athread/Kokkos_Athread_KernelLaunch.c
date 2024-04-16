#include "Kokkos_Athread_ParamWrap.h"
#include "simd.h"
#include "athread.h"

// DMA: 128B, sizeof(uint) = 4B
// 128 / 4 = 32
// extern uint* g_athread_functor_key __attribute__ ((aligned(64)));
// uintv16: 16 * uint
// intv16* g_athread_functor_key_simd __attribute__ ((aligned(64))) = 
//     (intv16 *)libc_aligned_malloc(2 * sizeof(intv16));

extern void SLAVE_FUN(register_kernel)();
void launch_register_kernel () {
	athread_spawn(register_kernel, 0);
	athread_join();
	return ;
}

// parallel_for
extern void SLAVE_FUN(parallel_for_1D)(struct AthreadParamWrap*);
void athread_parallel_for_launch_1D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_for_1D, param);
	athread_join();
	return ;
}

extern void SLAVE_FUN(parallel_for_2D)(struct AthreadParamWrap*);
void athread_parallel_for_launch_2D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_for_2D, param);
	athread_join();
	return ;
}

extern void SLAVE_FUN(parallel_for_3D)(struct AthreadParamWrap*);
void athread_parallel_for_launch_3D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_for_3D, param);
	athread_join();
	return ;
}
extern void SLAVE_FUN(parallel_for_4D)(struct AthreadParamWrap*);
void athread_parallel_for_launch_4D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_for_4D, param);
	athread_join();
	return ;
}

// parallel_reduce
extern void SLAVE_FUN(parallel_reduce_1D)(struct AthreadParamWrap*);
void athread_parallel_reduce_launch_1D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_reduce_1D, param);
	athread_join();
	return ;
}
extern void SLAVE_FUN(parallel_reduce_2D)(struct AthreadParamWrap*);
void athread_parallel_reduce_launch_2D (struct AthreadParamWrap* param) {
	athread_spawn(parallel_reduce_2D, param);
	athread_join();
	return ;
}