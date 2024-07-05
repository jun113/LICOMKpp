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

#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate uv");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_2(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_u, CppPOPHaloMod::arrCommPriorK,
        0, 4, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_u,
        0, 4, KM, JMT, IMT);

    gpu_get_halo_transpose_bclinc (*p_v_v, CppPOPHaloMod::arrCommPriorK,
        0, 4, KM, JMT, IMT);
    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_v,
        0, 4, KM, JMT, IMT);
// #elif (defined KOKKOS_ENABLE_ATHREAD)
//     athread_get_halo_transpose_double_host ((*p_v_u).data(), CppPOPHaloMod::arrCommPriorK,
//         2, 2, KM, JMT, IMT);
//     pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
//         KM, JMT, IMT, 
//           CppDomain::POP_haloClinic_C, 
//           CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
//           CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
//     athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_u).data(), 
//         0, 2, KM, JMT, IMT);
 
//     athread_get_halo_transpose_double_host ((*p_v_v).data(), CppPOPHaloMod::arrCommPriorK,
//         2, 2, KM, JMT, IMT);
//     pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
//         KM, JMT, IMT, 
//           CppDomain::POP_haloClinic_C, 
//           CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
//           CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
//     athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_v).data(), 
//         0, 2, KM, JMT, IMT);
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
    parallel_for("bclinc_14", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorBclinc14());
    parallel_for("bclinc_15", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclinc15());
    // -----------------------------
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_3(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_wka, CppPOPHaloMod::arrCommPriorK,
        0, 4, KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_wka,
        0, 4, KM, JMT, IMT);
// #elif (defined KOKKOS_ENABLE_ATHREAD)
//     // Athread LDM
//     athread_get_halo_transpose_double_host ((*p_v_wka).data(), CppPOPHaloMod::arrCommPriorK,
//         2, 2, KM, JMT, IMT);
 
//     pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
//         KM, JMT, IMT, 
//         CppDomain::POP_haloClinic_C, 
//         CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
//         CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
 
//     athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_wka).data(), 
//         0, 2, KM, JMT, IMT);
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
    parallel_for ("bclinc_16", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc16());

    parallel_for ("bclinc_17", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc17());

    parallel_for ("bclinc_18", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclinc18());
    // -----------------------------
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate wka");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_3(KM, 2);
    gpu_get_halo_transpose_bclinc (*p_v_wka, CppPOPHaloMod::arrCommPriorK,
        0, 4, KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_wka,
        0, 4, KM, JMT, IMT);
// #elif (defined KOKKOS_ENABLE_ATHREAD)
//     // Athread LDM
//     athread_get_halo_transpose_double_host ((*p_v_wka).data(), CppPOPHaloMod::arrCommPriorK,
//         2, 2, KM, JMT, IMT);

//     pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
//         KM, JMT, IMT, 
//         CppDomain::POP_haloClinic_C, 
//         CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
//         CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

//     athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_wka).data(), 
//         0, 2, KM, JMT, IMT);
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

    parallel_for ("bclinc_19", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclinc19());

    parallel_for ("bclinc_20", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclinc20());
  }
#else // BCLINC_MERGED_HALO
  // Merged
  {
#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_start("bclinc_merged_14");
#endif // LICOM_ENABLE_TEST_BCLINC
    parallel_for("bclinc_merged_14", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorBclincMer14());

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_stop ("bclinc_merged_14");
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_start("bclinc_merged_invtriu");
#endif // LICOM_ENABLE_TEST_BCLINC

#if (!defined KOKKOS_ENABLE_ATHREAD)
    parallel_for("bclinc_merged_15", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclincMer15());

    parallel_for("bclinc_merged_16", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBclincMer16());
#else
  athread_invtriu_host (KM, JMT, IMT, 
      CppPconstMod::dtc2, aidif,
      (*p_v_kmu).data(), 
      (*p_v_work_merged).data(),
      (*p_v_sbcx).data(), 
      (*p_v_bbcx).data(), 
      (*p_v_akmu).data(), 
      (*p_v_odzt).data(), 
      (*p_v_odzp).data(), 
      (*p_v_viv).data());

  athread_invtriu_host (KM, JMT, IMT, 
      CppPconstMod::dtc2, aidif,
      (*p_v_kmu).data(), 
      &((*p_v_work_merged).data()[KM * JMT * IMT]),
      (*p_v_sbcy).data(), 
      (*p_v_bbcy).data(), 
      (*p_v_akmu).data(), 
      (*p_v_odzt).data(), 
      (*p_v_odzp).data(), 
      (*p_v_viv).data());
#endif

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_stop("bclinc_merged_invtriu");
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc haloupdate merged");
#endif // LICOM_ENABLE_TEST_BCLINC
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_bclinc_3(KM, 2);
    // gpu_get_halo_transpose_bclinc (*p_v_wka, CppPOPHaloMod::arrCommPriorK,
    //     2, 2, KM, JMT, IMT);

    // pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
    //     KM, JMT, IMT, 
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    // gpu_put_halo_transpose_bclinc (CppPOPHaloMod::arrCommPriorK, *p_v_wka,
    //     0, 2, KM, JMT, IMT);
#elif (defined KOKKOS_ENABLE_ATHREAD)
    // Athread LDM
    athread_get_halo_transpose_double_host ((*p_v_wka).data(), CppPOPHaloMod::arrCommPriorK,
        0, 4, 2 * KM, JMT, IMT);

    pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
        2 * KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);

    athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_wka).data(), 
        0, 2, 2 * KM, JMT, IMT);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_work_merged).data(), KM << 1, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#endif
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc haloupdate merged");
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_start("bclinc_merged_vinteg");
#endif // LICOM_ENABLE_TEST_BCLINC

#if (!defined KOKKOS_ENABLE_ATHREAD)
    parallel_for ("bclinc_merged_vinteg_17", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBclincMer17());
#else
  athread_vinteg_host (KM, JMT, IMT,
      (*p_v_work_merged).data(),
      (*p_v_work1).data(),
      (*p_v_dzp).data(),
      (*p_v_viv).data(),
      (*p_v_ohbu).data());
  athread_vinteg_host (KM, JMT, IMT,
      &((*p_v_work_merged).data()[KM * JMT * IMT]),
      (*p_v_work2).data(),
      (*p_v_dzp).data(),
      (*p_v_viv).data(),
      (*p_v_ohbu).data());
#endif

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_stop ("bclinc_merged_vinteg");
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_start("bclinc_merged_18");
#endif // LICOM_ENABLE_TEST_BCLINC

    parallel_for ("bclinc_merged_18", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclincMer18());

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_stop ("bclinc_merged_18");
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_start("bclinc_merged_19");
#endif // LICOM_ENABLE_TEST_BCLINC

    parallel_for ("bclinc_merged_19", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBclincMer19());

#ifdef LICOM_ENABLE_TEST_BCLINC
          my_time.testTime_stop ("bclinc_merged_19");
#endif // LICOM_ENABLE_TEST_BCLINC
  }
#endif // BCLINC_MERGED_HALO

  } // End if (isc >= 1)
  isc += 1;

//   fortran_mpi_barrier_();

  return ;
}
//--------------------------------------
// End bclinc
//--------------------------------------
#endif // LICOM_ENABLE_KOKKOS
