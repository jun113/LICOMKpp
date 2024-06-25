#include "head/def-undef.h"
#include "head/cpp_param_mod.h"
#include "head/cpp_pconst_mod.h"
#include "head/cpp_extern_functions.h"
#include "head/fortran_extern_functions.h"
#include "head/licom_test_time.hpp"

#include "head/kokkos_config.hpp"

#ifndef LICOM_ENABLE_FORTRAN
#include "Kokkos_Core.hpp"
#endif // LICOM_ENABLE_FORTRAN

#define LICOM_ENABLE_NCU
#undef  LICOM_ENABLE_NCU

#define LICOM_ENABLE_RCU
#undef  LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_NCU
#include "nvToolsExt.h"
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
#include "roctx.h"
#endif // LICOM_ENABLE_RCU

#include <cstdio>

void stepon(const int &mm) {
  using CppPconstMod::ii;
  using CppPconstMod::jj;
  using CppPconstMod::ncc;
  using CppPconstMod::iday;
  using CppPconstMod::mon0;
  using CppPconstMod::iyfm;
  using CppPconstMod::dts_accum;
  using CppPconstMod::curr_ymd_licom;
  using CppPconstMod::num_step_per_day;

#ifdef LICOM_ENABLE_TEST_TIME
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_TIME 

  // if (ii == 1 && mm == 1 && iday == 1) {
  //  energy_();
  // }

//  if (ii == 5) {
//   exit(0);
//  }

  num_step_per_day = num_step_per_day + 1;

  curr_ymd_licom = iyfm * 10000 + mon0 * 100 + iday;   
//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("stepon\n");
// }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.readyt_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("readyt");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("readyt");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
  readyt_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
  cpp_readyt();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
  cuda_readyt();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
  hip_readyt();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
  // readyt_();
  kokkos_readyt();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.readyt_stop();
#endif // LICOM_ENABLE_TEST_TIME 
//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl readyt ok\n");
// }

  for (jj = 1; jj <= ncc; ++jj) {

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.readyc_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("readyc");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("readyc");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
    readyc_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
    cpp_readyc();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
    cuda_readyc();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
    hip_readyc();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
    // readyc_();
    kokkos_readyc();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.readyc_stop();
#endif // LICOM_ENABLE_TEST_TIME 

//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl readyc ok\n");
// }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.barotr_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("barotr");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("barotr");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
    barotr_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
    cpp_barotr();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
    cuda_barotr();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
    hip_barotr();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
    // barotr_();
    kokkos_barotr();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.barotr_stop();
#endif // LICOM_ENABLE_TEST_TIME 

//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl barotr ok\n");
// }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.bclinc_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("bclinc");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("bclinc");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
    bclinc_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
    cpp_bclinc();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
    cuda_bclinc();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
    hip_bclinc();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
    // bclinc_();
    kokkos_bclinc();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.bclinc_stop();
#endif // LICOM_ENABLE_TEST_TIME 

//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl bclinc ok\n");
// }
  // daily_update_d2h();
  // energy_();
  // -----------------------
} // End for (jj = 1; jj <= ncc; ++jj)

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.tracer_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("tracer");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("tracer");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
  tracer_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
  cpp_tracer();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
  cuda_tracer();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
  hip_tracer();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
  // tracer_();
  kokkos_tracer();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.tracer_stop();
#endif // LICOM_ENABLE_TEST_TIME 

//   fortran_mpi_barrier_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl tracer ok\n");
// }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.icesnow_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("icesnow");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("icesnow");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_FORTRAN
  icesnow_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
  cpp_icesnow();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
  cuda_icesnow();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
  hip_icesnow();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
  // icesnow_();
  kokkos_icesnow();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.icesnow_stop();
#endif // LICOM_ENABLE_TEST_TIME 

  daily_update_d2h();
  energy_();
//**********************

  if (ii == 1) {
#ifdef LICOM_ENABLE_TEST_TIME
  my_time.convadj_start();
#endif // LICOM_ENABLE_TEST_TIME 

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePushA ("convadj");
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePushA ("convadj");
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_FORTRAN
    convadj_();
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
    cpp_convadj();
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
    cuda_convadj();
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
    hip_convadj();
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
    // convadj_();
    kokkos_convadj();
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_NCU
  cudaDeviceSynchronize();
  nvtxRangePop();
#endif // LICOM_ENABLE_NCU

#ifdef LICOM_ENABLE_RCU
  hipDeviceSynchronize();
  roctxRangePop;
#endif // LICOM_ENABLE_RCU

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.convadj_stop();
#endif // LICOM_ENABLE_TEST_TIME 
  }
  if (dts_accum) {
    accumm_();
  }
  return ;
}
