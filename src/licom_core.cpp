#include "head/def-undef.h"
#include "head/cpp_blocks.h"
#include "head/cpp_domain.h"
#include "head/cpp_msg_mod.h"
#include "head/cpp_param_mod.h"
#include "head/cpp_pconst_mod.h"
#include "head/cpp_extern_functions.h"
#include "head/fortran_extern_functions.h"

#include "head/licom_test_time.hpp"

#ifdef LICOM_ENABLE_KOKKOS
#include "Kokkos_Core.hpp"
#include "head/kokkos_config.hpp"
#endif // LICOM_ENABLE_KOKKOS

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <cstring>

namespace MyTest {
TestTime::MyTime my_time;
} // namespace MyTest


extern "C" void licom_core_() {

  using CppMsgMod::nproc;

  using CppParamMod::mytid;

  using CppPconstMod::ii;

  using CppPconstMod::nss;
  using CppPconstMod::imd;
  using CppPconstMod::iday;
  using CppPconstMod::mon0;
  using CppPconstMod::month;
  using CppPconstMod::nmonth;
  using CppPconstMod::dd_licom;
  using CppPconstMod::dts_accum;
  using CppPconstMod::number_day;
  using CppPconstMod::number_month;
  using CppPconstMod::num_step_per_day;

  using MyTest::my_time;
  //---------------------------------------------

#ifndef LICOM_ENABLE_FORTRAN
  get_block_info_();
  using CppPOPHaloMod::pop_halo_create_from_fortran;

  pop_halo_create_from_fortran(CppDomain::POP_haloClinic_C);
  if (CppParamMod::mytid == 0) {
    printf ("In LICOM core\n");
    printf ("KM = %d, JMT = %d, IMT = %d\n",
        CppParamMod::KM, CppParamMod::JMT, CppParamMod::IMT);
    printf ("JST = %d, JET = %d\n",
        CppParamMod::JST, CppParamMod::JET);
    printf ("NY_BLOCK = %d, NX_BLOCK = %d\n",
        CppParamMod::NY_BLOCK, CppParamMod::NX_BLOCK);
    printf ("block ib = %d, ie = %d, jb = %d, je = %d\n",
        CppBlocks::ib, CppBlocks::ie,
        CppBlocks::jb, CppBlocks::je);
  }

  if (CppParamMod::mytid == 0) {
    printf("Create the data struct of POP Halo in C/C++ code done.\n");
  }
  // using CppPOPDistributionMod::pop_distribution_create;
  // pop_distribution_create(CppDomain::POP_distrbClinic);
  // using CppPOPHaloMod::pop_halo_create;
  // pop_halo_create(CppDomain::POP_distrbClinic, 
  //     CppDomain::ns_boundary_type, CppDomain::ew_boundary_type, 
  //         CppParamMod::IMT_GLOBAL, CppDomain::POP_haloClinic);

  // ========
#else
  // allocate_readyt_();
  // allocate_readyc_();
  // allocate_tracer_();
#endif // LICOM_ENABLE_FORTRAN

  int totalday = 0;

  my_time.testTime_initialize();

  if (mytid == 0) printf("********************************\n");
#ifdef LICOM_ENABLE_FORTRAN
  if (mytid == 0) {
    printf("stepon begin, version: %s (Fortran) %d CORES\n", LICOM_RES, nproc);
  }
#endif // LICOM_ENABLE_FORTRAN
#ifdef LICOM_ENABLE_CPP
  if (mytid == 0) {
    printf("stepon begin, version: %s (CPP)\n", LICOM_RES);
  }
#if defined(LOWRES)
  cpp_init_jra_daily_low();
#elif defined(HIGHRES) || defined(SUPHIGH)
  cpp_init_jra_daily_high();
#endif // defined(HIGHRES) || defined(SUPHIGH)
  if (CppParamMod::mytid == 0) {
    printf("Initialize jra daily in C/C++ code done.\n");
  }
#endif // LICOM_ENABLE_CPP
#ifdef LICOM_ENABLE_CUDA
  if (mytid == 0) {
    printf("stepon begin, version: %s (CUDA)\n", LICOM_RES);
  }
#endif // LICOM_ENABLE_CUDA
#ifdef LICOM_ENABLE_HIP
  if (mytid == 0) {
    printf("stepon begin, version: %s (HIP)\n", LICOM_RES);
  }
#endif // LICOM_ENABLE_HIP
#ifdef LICOM_ENABLE_KOKKOS
  if (mytid==0) {
    printf("stepon begin, version: %s (Kokkos)\n", LICOM_RES);
  }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_start("Kokkos and jra init");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.start_once();

  auto kokkos_init_args = Kokkos::InitializationSettings()
      .set_print_configuration(true)
      .set_map_device_id_by("mpi_rank")
      .set_tune_internals(true);
  
  // auto kokkos_init_args = Kokkos::InitializationSettings()
  //     .set_print_configuration(true)
  //     .set_device_id(mytid % 4)
  //     .set_tune_internals(true);

      //.set_device_id(mytid % 4)
      //.set_map_device_id_by("mpi_rank")

  Kokkos::initialize(kokkos_init_args);

#ifdef LICOM_ENABLE_KOKKOS
  atexit(Kokkos::finalize);
#endif // LICOM_ENABLE_KOKKOS

#ifdef LICOM_ENABLE_TEST_TIME
  // atexit(my_time.testTime_finalize);
#endif // LICOM_ENABLE_TEST_TIME

  kokkos_init_view();

  if (CppParamMod::mytid == 0) {
    printf("Initialize Kokkos Views done.\n");
  }

#if defined(LOWRES)
  // kokkos_init_jra_daily_high();
  kokkos_init_jra_daily_low();
#elif defined(HIGHRES) || defined(SUPHIGH)
  kokkos_init_jra_daily_high();
#endif // defined(HIGHRES) || defined(SUPHIGH)
  if (CppParamMod::mytid == 0) {
    printf("Initialize jra daily in Kokkos code done.\n");
  }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_stop("Kokkos and jra init");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.end_once();

  if (mytid == 0) {
    auto HOST = Kokkos::DefaultHostExecutionSpace();
    auto DEVICE = Kokkos::DefaultExecutionSpace();
    HOST.print_configuration(std::cout, true);
    DEVICE.print_configuration(std::cout, true);
    printf("Kokkos and jra initialize time: %.3f s\n", my_time.t_once);
  }
#endif // LICOM_ENABLE_KOKKOS

// =============================================
// Loop
#ifdef COUP 

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.fence();
#endif // LICOM_ENABLE_TEST_TIME
  my_time.start_stepon();
  for (ii = 1; ii <= nss/num_cpl; ++ii) {
#else // COUP

  for (int mm = 1; mm <= CppPconstMod::number; ++mm) {

    mon0 = (number_month - 1) % 12 + 1;
    imd  = nmonth[mon0 - 1];
    if (mytid == 0) {
      printf("========================\n");
      printf("test-mm: imd = %d, mon0 = %d, month = %d\n", imd, mon0, month);
    }
    for (iday = number_day; iday <= imd; ++iday) {
      if (mytid == 0) {
        printf("-------------------------------------\n");
      }
#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_start("daily");
      my_time.testTime_start("jra daily");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.start_daily();
      my_time.start_once();
#if (defined LICOM_ENABLE_FORTRAN)
      jra_daily_();
#endif
#if (defined LICOM_ENABLE_CPP)
#if defined(LOWRES)
      cpp_jra_daily_low(iday);
#elif defined(HIGHRES) || defined(SUPHIGH)
      cpp_jra_daily_high(iday);
#endif // defined(HIGHRES) || defined(SUPHIGH)
#endif // (defined LICOM_ENABLE_CPP)
#if (defined LICOM_ENABLE_KOKKOS)
#if defined(LOWRES)
      // kokkos_jra_daily_high(iday);
      kokkos_jra_daily_low(iday);
#elif defined(HIGHRES) || defined(SUPHIGH)
      kokkos_jra_daily_high(iday);
#endif // defined(HIGHRES) || defined(SUPHIGH)
#endif // (defined LICOM_ENABLE_KOKKOS)
#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_stop("jra daily");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.end_once();
      if (mytid == 0) {
        printf("jra_daily time: %.3f s\n", my_time.t_once);
      }

#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_start("stepon");
      //my_time.testTime_start("daily_h2d");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.start_stepon();

      //my_time.start_daily_h2d();

      //daily_update_h2d();

#ifdef LICOM_ENABLE_TEST_TIME
      //my_time.testTime_stop("daily_h2d");
#endif // LICOM_ENABLE_TEST_TIME
      //my_time.end_daily_h2d();

      for (ii = 1; ii <= nss; ++ii) {
#endif // CPUP
        // fortran_mpi_barrier_();
        // if (CppParamMod::mytid == 0) {
        //   printf("ii = %d\n", ii);
        // }
        stepon(mm);
      } // loop ii

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_start("daily_d2h");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.start_daily_d2h();
      // Copy back to do addps, nextstep
      daily_update_d2h();
#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_stop("daily_d2h");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.end_daily_d2h();
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

      num_step_per_day = nss;

      energy_();

#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_stop("stepon");
#endif // LICOM_ENABLE_TEST_TIME
      my_time.end_stepon();
      if (mytid == 0) {
        printf("stepon time: %.3f s (SYPD: %.2f), loops: %d\n", 
            my_time.t_stepon, 86400.0 / (my_time.t_stepon * 365.0), nss);
      }

      if (num_step_per_day == nss) {
        num_step_per_day = 0;
        addps_();
      }
      if ((!dts_accum) && (num_step_per_day == 0)) {
        accumm_();
      }
      if (iday == imd && num_step_per_day == 0) {
        ++month;
        dd_licom = 1;
      }
#ifdef COUP
      call shr_sys_flush(6);
#endif // COUP
      //LOGMSG
      if (num_step_per_day == 0) {

#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_start("write file");
#endif // LICOM_ENABLE_TEST_TIME

        my_time.start_once();

        ssaveins_();

#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_stop("write file");
#endif // LICOM_ENABLE_TEST_TIME
        my_time.end_once();

        if (mytid == 0) {
          printf("write_file time: %.3f s\n", my_time.t_once);
        }
#ifdef LICOM_ENABLE_KOKKOS
        kokkos_nextstep();
#else
        nextstep_();
#endif
      }

#ifdef LICOM_ENABLE_TEST_TIME
      my_time.testTime_stop("daily");
#endif // LICOM_ENABLE_TEST_TIME

      my_time.end_daily();

      if (mytid == 0) {
        printf("day: %.3f s, SYPD: %.2f , MM: %d, iday: %d\n", 
            my_time.t_daily, 86400.0 / (my_time.t_daily * 365.0), mm, iday);
      }
      //if (iday == 7) {
      totalday++;
      if (totalday == 2) {
        my_time.fence();
#ifdef LICOM_ENABLE_KOKKOS
        //Kokkos::finalize();
#endif // LICOM_ENABLE_KOKKOS
#ifdef LICOM_ENABLE_TEST_TIME
        my_time.testTime_finalize();
#endif // LICOM_ENABLE_TEST_TIME
        my_time.fence();
        exit(0);
      }
//--------------------------
    } // loop: iday
  }   // loop: mm
//-----------------------------
#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_finalize();
#endif // LICOM_ENABLE_TEST_TIME
#ifdef LICOM_ENABLE_KOKKOS
  Kokkos::finalize();
#endif // LICOM_ENABLE_KOKKOS

  return ;
}
