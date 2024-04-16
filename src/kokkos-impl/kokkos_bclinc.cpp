#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_bclinc.hpp"
//---------------------------------------------
//    BCLINC
void kokkos_bclinc() {

using Kokkos::parallel_for;
using Kokkos::MDRangePolicy;

#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BCLINC
#define LICOM_ENABLE_TEST_BCLINC
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_BCLINC
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_BCLINC

  double aa = 0.0;
  if (isc != 0) { 
    aa = 0.5;
  }

  parallel_for ("bclinc_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorBclinc1());

  parallel_for ("bclinc_2", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorBclinc2());

  parallel_for ("bclinc_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorBclinc3());

  parallel_for ("bclinc_4", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc4(aa));

  parallel_for ("bclinc_5", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorBclinc5());

  parallel_for ("bclinc_6", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc6());

  parallel_for ("bclinc_7", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorBclinc7());

  parallel_for ("bclinc_8", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorBclinc8());

  if (isc < 1) {

    parallel_for ("bclinc_9", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorBclinc9());

    parallel_for ("bclinc_10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclinc10());

    parallel_for ("bclinc_11", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclinc11());

#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate uv");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_2(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_u, CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_u,
        0, 2, KM, JMT, IMT);

    gpu_get_halo_transpose_bclinc (*p_v_v, CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_v,
        0, 2, KM, JMT, IMT);
#elif (defined KOKKOS_ENABLE_ATHREAD)
    athread_get_halo_transpose_double_host ((*p_v_u).data(), CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
          CppDomain::POP_haloClinic_C, 
          CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
          CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_u).data(), 
        0, 2, KM, JMT, IMT);
 
    athread_get_halo_transpose_double_host ((*p_v_v).data(), CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
          CppDomain::POP_haloClinic_C, 
          CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
          CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_v).data(), 
        0, 2, KM, JMT, IMT);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_u).data(), KM, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    CppPOPHaloMod::pop_halo_update((*p_v_v).data(), KM, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#endif
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc haloupdate uv");
#endif // LICOM_ENABLE_TEST_BCLINC

    parallel_for ("bclinc_12", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc12());

    parallel_for ("bclinc_13", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc13());

  } else {

#define BCLINC_MERGED_HALO
#undef  BCLINC_MERGED_HALO

#ifndef BCLINC_MERGED_HALO
  // Original
  {
    parallel_for("bclinc_14", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclinc14());
    // -----------------------------
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_3(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_wka, CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_wka,
        0, 2, KM, JMT, IMT);
#elif (defined KOKKOS_ENABLE_ATHREAD)
    // Athread LDM
    athread_get_halo_transpose_double_host ((*p_v_wka).data(), CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);
 
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
 
    athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_wka).data(), 
        0, 2, KM, JMT, IMT);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_wka).data(), KM, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#endif
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
    // -----------------------------
    parallel_for ("bclinc_15", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc15());

    parallel_for ("bclinc_16", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc16());

    parallel_for ("bclinc_17", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc17());
    // -----------------------------
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_3(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_wka, CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_wka,
        0, 2, KM, JMT, IMT);
#elif (defined KOKKOS_ENABLE_ATHREAD)
    // Athread LDM
    athread_get_halo_transpose_double_host ((*p_v_wka).data(), CppPOPHaloMod::arrCommPriorK,
        2, 2, KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_wka).data(), 
        0, 2, KM, JMT, IMT);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_wka).data(), KM, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#endif
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
    // -----------------------------

    parallel_for ("bclinc_18", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc18());

    parallel_for ("bclinc_19", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc19());
  }
#else // BCLINC_MERGED_HALO
  // Merged
  {
    parallel_for("bclinc_14", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{2, 2, 0}, koArr3D{IMT-2, JMT-2, KM}, tile3D), FunctorBclinc141());

    parallel_for("bclinc_15", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{IMT-2, JMT-2}, tile2D), FunctorBclinc151());
    parallel_for("bclinc_16", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{IMT-2, JMT-2}, tile2D), FunctorBclinc161());
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate");
#endif // LICOM_ENABLE_TEST_BCLINC
    pop_haloupdate_bclinc_33(KM, 2);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc haloupdate");
#endif // LICOM_ENABLE_TEST_BCLINC

    parallel_for ("bclinc_17", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{IMT, JMT}, tile2D), FunctorBclinc171());

    parallel_for ("bclinc_18", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{IMT, JMT, KM}, tile3D), FunctorBclinc181());
    parallel_for ("bclinc_19", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{IMT, JMT, KM}, tile3D), FunctorBclinc191());
  }
#endif // BCLINC_MERGED_HALO

  }
  ++isc;

//   fortran_mpi_barrier_();

  // if (CppParamMod::mytid == 0) {
  //   std::cout<<"FunctorBclinc1:  "<<sizeof(FunctorBclinc1)<<std::endl;
  //   std::cout<<"FunctorBclinc2:  "<<sizeof(FunctorBclinc2)<<std::endl;
  //   std::cout<<"FunctorBclinc3:  "<<sizeof(FunctorBclinc3)<<std::endl;
  //   std::cout<<"FunctorBclinc4:  "<<sizeof(FunctorBclinc4)<<std::endl;
  //   std::cout<<"FunctorBclinc5:  "<<sizeof(FunctorBclinc5)<<std::endl;
  //   std::cout<<"FunctorBclinc6:  "<<sizeof(FunctorBclinc6)<<std::endl;
  //   std::cout<<"FunctorBclinc7:  "<<sizeof(FunctorBclinc7)<<std::endl;
  //   std::cout<<"FunctorBclinc8:  "<<sizeof(FunctorBclinc8)<<std::endl;
  //   std::cout<<"FunctorBclinc9:  "<<sizeof(FunctorBclinc9)<<std::endl;
  //   std::cout<<"FunctorBclinc10: "<<sizeof(FunctorBclinc10)<<std::endl;
  //   std::cout<<"FunctorBclinc11: "<<sizeof(FunctorBclinc11)<<std::endl;
  //   std::cout<<"FunctorBclinc12: "<<sizeof(FunctorBclinc12)<<std::endl;
  //   std::cout<<"FunctorBclinc13: "<<sizeof(FunctorBclinc13)<<std::endl;
  //   std::cout<<"FunctorBclinc14: "<<sizeof(FunctorBclinc14)<<std::endl;
  //   std::cout<<"FunctorBclinc15: "<<sizeof(FunctorBclinc15)<<std::endl;
  //   std::cout<<"FunctorBclinc16: "<<sizeof(FunctorBclinc16)<<std::endl;
  //   std::cout<<"FunctorBclinc17: "<<sizeof(FunctorBclinc17)<<std::endl;
  //   std::cout<<"FunctorBclinc18: "<<sizeof(FunctorBclinc18)<<std::endl;
  //   std::cout<<"FunctorBclinc19: "<<sizeof(FunctorBclinc19)<<std::endl;

  //   std::cout<<"FunctorBclinc141: "<<sizeof(FunctorBclinc141)<<std::endl;
  //   std::cout<<"FunctorBclinc151: "<<sizeof(FunctorBclinc151)<<std::endl;
  //   std::cout<<"FunctorBclinc161: "<<sizeof(FunctorBclinc161)<<std::endl;
  //   std::cout<<"FunctorBclinc171: "<<sizeof(FunctorBclinc171)<<std::endl;
  //   std::cout<<"FunctorBclinc181: "<<sizeof(FunctorBclinc181)<<std::endl;
  //   std::cout<<"FunctorBclinc191: "<<sizeof(FunctorBclinc191)<<std::endl;
  // }

  return ;
}
//--------------------------------------
// End bclinc
//--------------------------------------
#endif // LICOM_ENABLE_KOKKOS
