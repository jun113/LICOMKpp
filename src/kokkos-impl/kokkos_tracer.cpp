#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_tracer.hpp"

#include <mpi.h>

// TRACER
void kokkos_tracer() {

//   using CppDomain  ::nblocks_clinic;
  using CppParamMod::mytid;
#ifndef ISO
//   using CppDomain::blocks_clinic;
#endif // ISO
  using CppGrid     ::area_t;
//   using CppGrid     ::horiz_grid_opt;
  using CppPconstMod::dts;
  using CppPconstMod::ist;
//   using CppPconstMod::adv_tracer;
  using CppPconstMod::boundary_restore;
  using CppTracerMod::fw_norm2;

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_TRACER
#define LICOM_ENABLE_TEST_TRACER
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_TRACER
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_TRACER

#ifdef ISO
#ifdef LDD97
  ViewDouble4D v_f1("view_f1", IMT, JMT, KM, MAX_BLOCKS_CLINIC);
  ViewDouble4D v_f2("view_f2", IMT, JMT, KM, MAX_BLOCKS_CLINIC);
#endif // LDD97
#endif // ISO

  double aa, c2dtts;

  int flag_adv_tracer = 0;

/*
  // TODO: 419 new bug, wjl 20230328
  const std::string str_adv_tracer(adv_tracer);

  const char* adv_tracer_centered = "centered";
  const char* adv_tracer_flux     = "flux";
  const char* adv_tracer_tspas    = "tspas";
  if (str_trim_cmp(adv_tracer, adv_tracer_centered) == 0) {
    flag_adv_tracer = 1;
    if (ist >= 1) {
      aa = 0.5;
      c2dtts = dts * 2.0;
    } else {
      aa = 0.0;
      c2dtts = dts;
    }
  } else if (str_trim_cmp(adv_tracer, adv_tracer_flux) == 0) {
    flag_adv_tracer = 2;
  } else if (str_trim_cmp(adv_tracer, adv_tracer_tspas) == 0) {
    flag_adv_tracer = 3;
    aa = 0.5;
    c2dtts = dts;
  } else {
    if (mytid == 0) {
      printf ("%s, %d\n", __FILE__, __LINE__);
      printf ("The false advection option for tracer\n");
    }
    exit(0);
  }
*/
  // =======================================================
  // Needed for 419
  flag_adv_tracer = 3;

  if (flag_adv_tracer == 1) {
    if (ist >= 1) {
      aa = 0.5;
      c2dtts = dts * 2.0;
    } else {
      aa = 0.0;
      c2dtts = dts;
    }
  } else if (flag_adv_tracer == 2) {
  } else if (flag_adv_tracer == 3) {
    aa = 0.5;
    c2dtts = dts;
  } else {
    if (mytid == 0) {
      printf ("%s, %d\n", __FILE__, __LINE__);
      printf ("The false advection option for tracer\n");
    }
    exit(0);
  }
  // =======================================================

  parallel_for ("tracer_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer1(aa));

  parallel_for ("tracer_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorTracer2(aa));

  parallel_for ("tracer_3_upwell_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer3());
  parallel_for ("tracer_4_upwell_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 1}, koArr3D{KM, JMT-1, IMT}, tile3D), FunctorTracer4());
  parallel_for ("tracer_5_upwell_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorTracer5());
  parallel_for ("tracer_6_upwell_4", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer6());
#ifdef NODIAG

#ifdef ISO
  // isopyc
  parallel_for("tracer_isopyc_1", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), functor_tracer_isopyc_1());

  parallel_for("tracer_isopyc_2", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), functor_tracer_isopyc_2());

  // k2_3
  parallel_for("tracer_isopyc_k2_3_1", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 1, 1}, koArr3D{JMT-1, KM-1, IMT-1}, tile3D), functor_tracer_isopyc_k2_3_1());

  parallel_for("tracer_isopyc_k2_3_2", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_k2_3_2());

  parallel_for("tracer_isopyc_k2_3_3", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), functor_tracer_isopyc_k2_3_3());

  parallel_for("tracer_isopyc_k2_3_4", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_k2_3_4());

  parallel_for("tracer_isopyc_k2_3_5", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k2_3_5());

  parallel_for("tracer_isopyc_k2_3_6", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k2_3_6());
#ifdef LDD97
  parallel_for("tracer_isopyc_k2_3_7", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), 
          functor_tracer_isopyc_k2_3_7(v_f1, v_f2));

  parallel_for("tracer_isopyc_k2_3_8", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), 
          functor_tracer_isopyc_k2_3_8(v_f1, v_f2));
#else  // LDD97
  parallel_for("tracer_isopyc_k2_3_9", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k2_3_9());
#endif // LDD97
  // End k2_3

  // k1_3
  parallel_for("tracer_isopyc_k1_3_1", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 1, 1}, koArr3D{JMT-1, KM-1, IMT-1}, tile3D), functor_tracer_isopyc_k1_3_1());
    
  parallel_for("tracer_isopyc_k1_3_2", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 0}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_k1_3_2());

  parallel_for("tracer_isopyc_k1_3_3", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), functor_tracer_isopyc_k1_3_3());

  parallel_for("tracer_isopyc_k1_3_4", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{0, 0}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_k1_3_4());

  parallel_for("tracer_isopyc_k1_3_5", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k1_3_5());

  parallel_for("tracer_isopyc_k1_3_6", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 0}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k1_3_6());
#ifdef LDD97
  parallel_for("tracer_isopyc_k1_3_7", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), 
          functor_tracer_isopyc_k1_3_7(v_f1, v_f2));

  parallel_for("tracer_isopyc_k1_3_8", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 0}, koArr3D{JMT-1, KM, IMT-1}, tile3D), 
          functor_tracer_isopyc_k1_3_8(v_f1, v_f2));
#else  // LDD97
  parallel_for("tracer_isopyc_k1_3_9", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 0}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k1_3_9());
#endif // LDD97
  // End k1_3

  // k3_123
  parallel_for("tracer_isopyc_k3_123_1", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 1, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k3_123_1());

  parallel_for("tracer_isopyc_k3_123_2", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 0}, koArr2D{JMT-1, IMT}, tile2D), functor_tracer_isopyc_k3_123_2());
#ifdef LDD97
  parallel_for("tracer_isopyc_k3_123_3", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{IMT, JMT, KM}, tile3D), 
          functor_tracer_isopyc_k3_123_3(v_f1, v_f2));

  parallel_for("tracer_isopyc_k3_123_4", MDRangePolicy<Kokkos::Rank<3>>
      (koArr2D{1, 0, 0}, koArr2D{JMT-1, KM, IMT-1}, tile3D), 
          functor_tracer_isopyc_k3_123_4(v_f1, v_f2));

#else  // LDD97
  parallel_for("tracer_isopyc_k3_123_5", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 0, 1}, koArr3D{JMT-1, KM, IMT-1}, tile3D), functor_tracer_isopyc_k3_123_5());
#endif // LDD97
  // End k3_123 

  // isoadv
  parallel_for("tracer_isopyc_isoadv_1", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{JMT, KM, IMT}, tile3D), functor_tracer_isopyc_isoadv_1());

  parallel_for("tracer_isopyc_isoadv_2", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{0, 0, 0}, koArr3D{JMT, KM+1, IMT}, tile3D), functor_tracer_isopyc_isoadv_2());

  parallel_for("tracer_isopyc_isoadv_3", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 1, 1}, koArr3D{JMT-1, KM-1, IMT-1}, tile3D), functor_tracer_isopyc_isoadv_3());

  parallel_for("tracer_isopyc_isoadv_4", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_isoadv_4());

  parallel_for("tracer_isopyc_isoadv_5", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_isoadv_5());

  parallel_for("tracer_isopyc_isoadv_6", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{1, 1, 1}, koArr3D{JMT-1, KM-1, IMT-1}, tile3D), functor_tracer_isopyc_isoadv_6());

  parallel_for("tracer_isopyc_isoadv_7", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_isoadv_7());

  parallel_for("tracer_isopyc_isoadv_8", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), functor_tracer_isopyc_isoadv_8());

  parallel_for("tracer_isopyc_isoadv_9", MDRangePolicy<Kokkos::Rank<3>>
      (koArr3D{2, 0, 2}, koArr3D{JMT-2, KM-1, IMT-2}, tile3D), functor_tracer_isopyc_isoadv_9());

  parallel_for("tracer_isopyc_isoadv_10", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), functor_tracer_isopyc_isoadv_10());

  parallel_for("tracer_isopyc_isoadv_11", MDRangePolicy<Kokkos::Rank<2>>
      (koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), functor_tracer_isopyc_isoadv_11());
  // End isoadv
  // End isopyc
#endif // ISO

  parallel_for ("tracer_7", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorTracer7());

  //======================
  // !---------------------------------------------
  // !     SOLVE FOR ONE TRACER AT A TIME
  // !---------------------------------------------
  // !     NTRA = 1 => TEMPERATURE
  // !     NTRA = 2 => SALINITY
  for (int n = 0; n < NTRA; ++n) {

    parallel_for ("tracer_8", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorTracer8());

// #ifdef LICOM_ENABLE_TEST_TRACER
//       my_time.testTime_start("advection tracer");
// #endif // LICOM_ENABLE_TEST_TRACER
    // for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // advection_tracer(wkd[iblock], wkb[iblock], ws[iblock], at[iblock][n], adv_tt, iblock, n)
      if (flag_adv_tracer == 1 || flag_adv_tracer == 3) {
        parallel_for ("tracer_9_advcetion_tracer_centered_tspas_1", 
            MDRangePolicy<Kokkos::Rank<3>> (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), 
                FuncAdvTraCenTsp1());
      } else if (flag_adv_tracer == 2) {
        // parallel_for ("tracer_9_advcetion_tracer_flux_1", MDRangePolicy<Kokkos::Rank<3>>
        //     ({0, 0, 1}, {KM, JMT, IMT-1}, tile3D), 
        //         functor_tracer_advection_tracer_3(0));
        // parallel_for("tracer_9_advcetion_tracer_flux_2", MDRangePolicy<Kokkos::Rank<3>>
        //     ({0, 1, 0}, {KM, JMT-1, IMT}, tile3D), 
        //         functor_tracer_advection_tracer_4(0));
      }

      if (flag_adv_tracer == 1) {
        // parallel_for ("tracer_10_advcetion_tracer_centered_2", MDRangePolicy<Kokkos::Rank<2>>
        //     ({2, 2}, {JMT-2, IMT-2}, tile2D), 
        //         functor_tracer_advection_tracer_5(n, 0));
        // parallel_for ("tracer_11_advcetion_tracer_centered_3", MDRangePolicy<Kokkos::Rank<3>>
        //     ({1, 2, 2}, {KM-1, JMT-2, IMT-2}, tile3D), 
        //         functor_tracer_advection_tracer_6(n, 0));
        // parallel_for ("tracer_12_advcetion_tracer_centered_4", MDRangePolicy<Kokkos::Rank<2>>
        //     ({2, 2}, {JMT-2, IMT-2}, tile2D), 
        //         functor_tracer_advection_tracer_7(n, 0));
      } else if (flag_adv_tracer == 2) {
        // parallel_for ("tracer_10_advcetion_tracer_flux_3", MDRangePolicy<Kokkos::Rank<3>>
        //     ({0, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), 
        //         functor_tracer_advection_tracer_8(n, 0));
      } else if (flag_adv_tracer == 3) {
        parallel_for ("tracer_10_advcetion_tracer_tspas_2", 
            MDRangePolicy<Kokkos::Rank<3>>(koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), 
                FuncAdvTraTsp2(n));

        parallel_for ("tracer_11_advcetion_tracer_tspas_3", 
            MDRangePolicy<Kokkos::Rank<3>>(koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), 
                FuncAdvTraTsp3(n));

        // k = 0
        parallel_for ("tracer_12_advcetion_tracer_tspas_4", 
            MDRangePolicy<Kokkos::Rank<2>> (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), 
                FuncAdvTraTsp4(n));
        // 0 < k < KM
        parallel_for ("tracer_13_advction_tracer_5", 
            MDRangePolicy<Kokkos::Rank<3>> (koArr3D{1, 1, 1}, koArr3D{KM-1, JMT-1, IMT-1}, tile3D), 
                FuncAdvTraTsp5(n));
        // k = KM-1
        parallel_for ("tracer_14_advction_tracer_6", 
            MDRangePolicy<Kokkos::Rank<2>> (koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), 
                FuncAdvTraTsp6(n));
      }
      // End advection_tracer
// #ifdef LICOM_ENABLE_TEST_TRACER
//       my_time.testTime_stop("advection tracer");
// #endif // LICOM_ENABLE_TEST_TRACER
    
      parallel_for ("tracer_15", MDRangePolicy<Kokkos::Rank<3>> (
          koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorTracer15());
    // }

#ifdef CANUTO
    parallel_for ("tracer_16", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorTracer16(n));
#endif // CANUTO

#ifdef ISO
    // isoflux(n)
    parallel_for("tracer_isoflux_1", MDRangePolicy<Kokkos::Rank<3>>
        ({1, 1, 1}, {KM-1, JMT-2, IMT-2}, tile3D), functor_tracer_isoflux_1(n));

    parallel_for("tracer_isoflux_2", MDRangePolicy<Kokkos::Rank<2>>
        ({1, 1}, {JMT-2, IMT-2}, tile2D), functor_tracer_isoflux_2(n));

    parallel_for("tracer_isoflux_3", MDRangePolicy<Kokkos::Rank<2>>
        ({0, 0}, {JMT, IMT}, tile2D), functor_tracer_isoflux_3());

    parallel_for("tracer_isoflux_4", MDRangePolicy<Kokkos::Rank<2>>
        ({1, 1}, {JMT-1, IMT-1}, tile2D), functor_tracer_isoflux_4(n));

    parallel_for("tracer_isoflux_5", MDRangePolicy<Kokkos::Rank<2>>
        ({0, 0}, {KM, IMT}, tile2D), functor_tracer_isoflux_5());

    parallel_for("tracer_isoflux_6", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), functor_tracer_isoflux_6(n));

    parallel_for("tracer_isoflux_7", MDRangePolicy<Kokkos::Rank<3>>
        ({1, 1, 1}, {KM-1, JMT-1, IMT-1}, tile3D), functor_tracer_isoflux_7(n));

    parallel_for("tracer_isoflux_8", MDRangePolicy<Kokkos::Rank<2>>
        ({1, 1}, {JMT-1, IMT-1}, tile2D), functor_tracer_isoflux_8(n));

    parallel_for("tracer_isoflux_9", MDRangePolicy<Kokkos::Rank<2>>
        ({0, 0}, {JMT, IMT}, tile2D), functor_tracer_isoflux_9());

    parallel_for("tracer_isoflux_10", MDRangePolicy<Kokkos::Rank<2>>
        ({1, 1}, {JMT-1, IMT-1}, tile2D), functor_tracer_isoflux_10(n));

    parallel_for("tracer_isoflux_11", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile2D), functor_tracer_isoflux_11(n));

    parallel_for("tracer_isoflux_12", MDRangePolicy<Kokkos::Rank<3>>
        ({1, 1, 1}, {KM, JMT-1, IMT-1}, tile2D), functor_tracer_isoflux_12(n));

    parallel_for("tracer_isoflux_13", MDRangePolicy<Kokkos::Rank<2>>
        ({0, 0}, {JMT, IMT}, tile2D), functor_tracer_isoflux_13());

    parallel_for("tracer_isoflux_14", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), functor_tracer_isoflux_14(n));

    parallel_for("tracer_isoflux_15", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), functor_tracer_isoflux_15(n));

    parallel_for("tracer_isoflux_16", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), functor_tracer_isoflux_16(n));

    parallel_for("tracer_isoflux_17", MDRangePolicy<Kokkos::Rank<2>>
        ({0, 0}, {JMT, IMT}, tile2D), functor_tracer_isoflux_17());

    parallel_for("tracer_isoflux_18", MDRangePolicy<Kokkos::Rank<3>>
        ({1, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), functor_tracer_isoflux_18(n));

    parallel_for("tracer_isoflux_19", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), functor_tracer_isoflux_19(n));
#else  // ISO
#ifdef SMAG
    call SMAG3
#else  // SMAG
    // const int iiblock = 0;
    // int block_id = blocks_clinic[iiblock];
    // int local_id = iiblock + 1;
    // const struct block this_block = CppBlocks::get_block(&block_id, &local_id);
    // const int ib = this_block.ib;
    // const int ie = this_block.ie;
    // const int jb = this_block.jb;
    // const int je = this_block.je;
#ifdef BIHAR
    parallel_for ("tracer_17_hdifft_del4_1", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, NY_BLOCK, NX_BLOCK}, tile3D), FunctorTracer17(n));

    parallel_for ("tracer_18_hdifft_del4_2", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorTracer18());
#else  // BIHAR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        // hdifft_del2(k, hdtk, atb[iblock][n][k+1], this_block);
        parallel_for("tracer_hdifft_del2_1", MDRangePolicy<Kokkos::Rank<2>>
            ({0, 0}, {NY_BLOCK, NX_BLOCK}, tile2D), 
                functor_tracer_hdifft_del2_1(k));

        parallel_for("tracer_hdifft_del2_2", MDRangePolicy<Kokkos::Rank<2>>
            ({0, 0}, {NY_BLOCK, NX_BLOCK}, tile2D), 
                functor_tracer_hdifft_del2_2());

        parallel_for("tracer_hdifft_del2_3", MDRangePolicy<Kokkos::Rank<2>>
            ({jb-1, ib-1}, {je, ie}, tile2D), 
                functor_tracer_hdifft_del2_3(n, k, iblock));

        // End hdifft_del2
        parallel_for("tracer_10", MDRangePolicy<Kokkos::Rank<2>>
            ({2, 2}, {JMT-2, IMT-2}, tile2D), 
                functor_tracer_10(n, k, iblock));
      }
    }
#endif // BIHAR
#endif //SMAG
#endif // iSO

    if (n == 0) {
      parallel_for ("tracer_19", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer19());
      parallel_for ("tracer_20", MDRangePolicy<Kokkos::Rank<3>> (
          koArr3D{1, 2, 2}, koArr3D{KM-1, JMT-2, IMT-2}, tile3D), FunctorTracer20());
      parallel_for ("tracer_21", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer21());
    }

    parallel_for ("tracer_22", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer22(n));
    parallel_for ("tracer_23", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{1, 2, 2}, koArr3D{KM-1, JMT-2, IMT-2}, tile3D), FunctorTracer23(n));
    parallel_for ("tracer_24", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer24(n));

    if (n == 1) {
      parallel_for ("tracer_25", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer25());
// #undef SSSNORM
#ifdef SSSNORM
      double err_norm1;
      double err_norm2 = 0.0;

#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
      Kokkos::parallel_reduce ("tracer_26", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer26(), err_norm2);
#else
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          err_norm2 += (*p_v_tarea)(0, j, i) 
              * (*p_v_net)(0, 1, j, i) * (*p_v_vit)(0, 0, j, i);
        }
      }
#endif

#ifdef SPMD

      MPI_Reduce (&err_norm2, &err_norm1, 1, MPI_DOUBLE,
          MPI_SUM, CppParamMod::master_task, 
              CppDomain::POP_haloClinic_C.communicator);
      MPI_Bcast (&err_norm1, 1, MPI_DOUBLE, CppParamMod::master_task,
          CppDomain::POP_haloClinic_C.communicator);

      err_norm2 = - err_norm1 / area_t;
#else  // SPMD
#endif // SPMD
      fw_norm2 = err_norm2;

      parallel_for ("tracer_27", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{JST-1, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer27(err_norm2));
#endif // SSSNORM
    } else {
      // for temperature
      parallel_for ("tracer_28", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer28());
      parallel_for ("tracer_29", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer29());
    }

    //--------------------
    // TODO: GPU and SW needs debug
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer haloupdate net");
#endif // LICOM_ENABLE_TEST_TRACER
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
  pop_haloupdate_tracer_net (NTRA, 3);
#elif (defined KOKKOS_ENABLE_ATHREAD)
  athread_get_halo_transpose_double_host ((*p_v_net).data(), CppPOPHaloMod::arrCommPriorK,
      2, 3, NTRA, JMT, IMT);

  pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
      NTRA, JMT, IMT, 
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);

  athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_net).data(), 
      0, 3, NTRA, JMT, IMT);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_net).data(), NTRA, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer haloupdate net");
#endif // LICOM_ENABLE_TEST_TRACER
    //--------------------
      
    //if (simple_assm) {
    //  if (mytid == 0) {
    //    printf("into restoring, n = %d\n", n+1);
    //  }
    //}

    // boundary_restore = 2
    // printf ("boundary_restore = %d\n", boundary_restore);
    if (boundary_restore == 1) {
      parallel_for("tracer_21", MDRangePolicy<Kokkos::Rank<3>>
          ({1, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), functor_tracer_21_22(n));
      parallel_for("tracer_22", MDRangePolicy<Kokkos::Rank<3>>
          ({1, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), functor_tracer_21_22(n));
    }

    if (flag_adv_tracer == 3) {

      parallel_for ("tracer_30", MDRangePolicy<Kokkos::Rank<3>> (
          koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorTracer30(n));

    } else if (flag_adv_tracer == 1) {

      parallel_for("tracer_24", MDRangePolicy<Kokkos::Rank<3>>
          ({0, 2, 2}, {KM, JMT-2, IMT-2}, tile3D), functor_tracer_24(n, c2dtts));

    } else {
      if (mytid == 0) {
        printf ("%s, %d\n", __FILE__, __LINE__);
        printf("The false advection option for tracer\n");
      }
      exit(0);
    }

    parallel_for("tracer_31", MDRangePolicy<Kokkos::Rank<3>>
        ({0, 0, 0}, {KM, JMT, IMT}, tile3D), FunctorTracer31());

    // invtrit(vtl, stf, wkc, aidif, c2dtts)
    parallel_for ("tracer_32_invtrit", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorTracer32 (c2dtts));
    // End invtrit
  
    //--------------------
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer haloupdate vtl");
#endif // LICOM_ENABLE_TEST_TRACER
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    // pop_haloupdate_tracer_2(KM);
  gpu_get_halo_transpose_tracer (*p_v_vtl, CppPOPHaloMod::arrCommPriorK,
      2, 3, KM, JMT, IMT);
  pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
      KM, JMT, IMT, 
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  gpu_put_halo_transpose_tracer (CppPOPHaloMod::arrCommPriorK, *p_v_vtl,
      0, 3, KM, JMT, IMT);
#elif (defined KOKKOS_ENABLE_ATHREAD)

  athread_get_halo_transpose_double_host ((*p_v_vtl).data(), CppPOPHaloMod::arrCommPriorK,
      2, 3, KM, JMT, IMT);

  pop_halo_update_priority_k (CppPOPHaloMod::arrCommPriorK,
      KM, JMT, IMT, 
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);

  athread_put_halo_transpose_double_host (CppPOPHaloMod::arrCommPriorK, (*p_v_vtl).data(), 
      0, 3, KM, JMT, IMT);

#else
    CppPOPHaloMod::pop_halo_update((*p_v_vtl).data(), KM, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer haloupdate vtl");
#endif // LICOM_ENABLE_TEST_TRACER
    //--------------------

    parallel_for ("tracer_33", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorTracer33 (n, c2dtts));
    parallel_for ("tracer_34", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{1, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorTracer34 (n, c2dtts));

    //const char* horiz_grid_opt_lat_lon  = "lat_lon";
    //if (str_trim_cmp(horiz_grid_opt, horiz_grid_opt_lat_lon) == 0) {
    // wjl 20230328
    // Needed for 419
    if (false) {
      printf ("tracer lat_lon\n");
      if (ist % 180 == 1) {
        // smts(vtl, vit, fil_lat2)
        parallel_for("tracer_smts_1_1", JMT, functor_tracer_smts_1());

        const double fil_lat2 = 63.0;
        parallel_for("tracer_smts_1_2", 
            Kokkos::RangePolicy<>(JMT-2, 2), functor_tracer_smts_2(fil_lat2));

        parallel_for("tracer_smts_1_3", MDRangePolicy<Kokkos::Rank<2>>
            ({2, 0}, {JMT-2, KM}), functor_tracer_smts_3());
        // End smts(vtl, vit, fil_lat2)

        //============================================
        // pop haloupdate
        //--------------------
        //pop_haloupdate_smts_(&errorCode);
        //--------------------
        // End pop haloupdate
        //============================================

        // End smts
      } else {
        if (flag_adv_tracer == 3) {
          parallel_for("tracer_27", MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {JMT, IMT}), functor_tracer_27(n));

          parallel_for("tracer_28", MDRangePolicy<Kokkos::Rank<3>>
              ({1, 0, 0}, {KM, JMT, IMT}), functor_tracer_28(n));

        } else if (flag_adv_tracer == 1) {

          parallel_for("tracer_29", MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {JMT, IMT}), functor_tracer_29(n, c2dtts));

          parallel_for("tracer_30", MDRangePolicy<Kokkos::Rank<3>>
              ({1, 0, 0}, {KM, JMT, IMT}), functor_tracer_30(n));

        } else {
          if (mytid == 0) {
            printf ("%s, %d\n", __FILE__, __LINE__);
            printf("The false advection option for tracer\n");
          }
          exit(0);
        }
        // smts(vtl, vit, fil_lat1)
        parallel_for("tracer_smts_2_1", JMT, functor_tracer_smts_1());

        const double fil_lat1 = 63.0;
        parallel_for("tracer_smts_2_2", 
            Kokkos::RangePolicy<>(2, JMT-2), functor_tracer_smts_2(fil_lat1));

        parallel_for("tracer_smts_2_3", MDRangePolicy<Kokkos::Rank<2>>
            ({2, 0}, {JMT-2, KM}), functor_tracer_smts_3());

        //============================================
        // pop haloupdate
        //--------------------
        //pop_haloupdate_smts_(&errorCode);
        //--------------------
        // End pop haloupdate
        //============================================

        // End smts
        if (flag_adv_tracer == 3) {

          parallel_for("tracer_31", MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {JMT, IMT}), functor_tracer_31(n));

          parallel_for("tracer_32", MDRangePolicy<Kokkos::Rank<3>>
              ({1, 0, 0}, {KM, JMT, IMT}), functor_tracer_32(n));

        } else if (flag_adv_tracer == 1) {

          parallel_for("tracer_33", MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {JMT, IMT}), functor_tracer_33(n, c2dtts));

          parallel_for("tracer_34", MDRangePolicy<Kokkos::Rank<3>>
              ({1, 0, 0}, {KM, JMT, IMT}), functor_tracer_34(n));

        } else {
          if (mytid == 0) {
            printf ("%s, %d\n", __FILE__, __LINE__);
            printf("The false advection option for tracer\n");
          }
          exit(0);
        }
      }
    } // false
    // =============

    if (flag_adv_tracer == 3) {
      parallel_for ("tracer_35", MDRangePolicy<Kokkos::Rank<3>> (
          koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorTracer35(n, c2dtts));

    } else if (flag_adv_tracer == 1) {

      parallel_for("tracer_35", MDRangePolicy<Kokkos::Rank<3>>
          (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), functor_tracer_35(n, c2dtts));

      if (ist >= 1) {
        parallel_for("tracer_37", MDRangePolicy<Kokkos::Rank<3>>
            (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), functor_tracer_37(n));
      }

      parallel_for("tracer_38", MDRangePolicy<Kokkos::Rank<3>>
          (koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), functor_tracer_38(n));

    } else {
      if (mytid == 0) {
        printf ("%s, %d\n", __FILE__, __LINE__);
        printf("The false advection option for tracer\n");
      }
      exit(0);
    }
//-------------------
  } // End Loop N

#else  // NODIAG
  parallel_for("tracer_39", MDRangePolicy<Kokkos::Rank<3>>
      ({0, 0, 0}, {KM+1, JMT, IMT}), functor_tracer_39());

  parallel_for("tracer_40", MDRangePolicy<Kokkos::Rank<3>>
      ({0, 0, 0}, {KM, JMT, IMT}), functor_tracer_40());

#endif // NODIAG

  ist += 1;

  return ;
}
#endif // LICOM_ENABLE_KOKKOS
