#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_readyc.hpp"

void kokkos_readyc () {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  using CppPconstMod::adv_momentum;

#ifdef BCKMEX
  ViewDouble3D v_diff_back("view_diff_back", IMT, JMT, MAX_BLOCKS_CLINIC);
  ViewDouble3D v_diff_back_sh("view_diff_back_sh", IMT, JMT, MAX_BLOCKS_CLINIC);
  ViewDouble3D v_diff_back_nn("view_diff_back_nn", IMT, JMT, MAX_BLOCKS_CLINIC);
#endif // BCKMEX

  parallel_for ("readyc_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc1());

  // parallel_for ("readyc_2", MDRangePolicy<Kokkos::Rank<3>> (
  //     koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc2());

  parallel_for ("readyc_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc3());

  parallel_for ("readyc_4", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KMM1, JMT, IMT}, tile3D), FunctorReadyc4());

#ifdef BCKMEX
  parallel_for("readyc_7", MDRangePolicy<Kokkos::Rank<2>>
      ({0, 0}, {IMT, JMT}), functor_readyc_7(v_diff_back, v_diff_back_sh, v_diff_nh));
  parallel_for("readyc_8", MDRangePolicy<Kokkos::Rank<2>>
      ({0, 0}, {IMT, JMT}), functor_readyc_8(v_diff_back, v_diff_back_sh, v_diff_nh));
#endif // BCKMEX

#ifdef CANUTO
  // parallel_for ("readyc_5", MDRangePolicy<Kokkos::Rank<2>> (
  //     koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorReadyc5());
  parallel_for ("readyc_51", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorReadyc51());
  // turb_2
  parallel_for ("readyc_52", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorReadyc52());
  // End turb_2
  parallel_for ("readyc_53", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM, JMT-1, IMT-1}, tile3D), FunctorReadyc53());
  parallel_for ("readyc_54", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorReadyc54());
  parallel_for ("readyc_55", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 1, 1}, koArr3D{KM-1, JMT-1, IMT-1}, tile3D), FunctorReadyc55());

  parallel_for ("readyc_6", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc6());
#endif // CANUTO

  parallel_for ("readyc_7_upwell_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc7());

  parallel_for ("readyc_8_upwell_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 1}, koArr3D{KM, JMT-1, IMT}, tile3D), FunctorReadyc8());

  parallel_for ("readyc_9_upwell_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc9());

  parallel_for ("readyc_10_upwell_4", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc10());

  parallel_for ("readyc_11", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc11());

  // advection_momentum(u, v, wka, dlu, dlv, iblock)
  const std::string str_adv_momentum(adv_momentum);
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    parallel_for ("readyc_12_advection_momentum_centered_1", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), 
            FuncAdvMomCen1());
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    parallel_for("readyc_12_advection_momentum_flux_1", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), 
            FuncAdvMomFlu1());
  } else {
    if (mytid == 0) {
      printf ("%s, %d\n", __FILE__, __LINE__);
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    parallel_for ("readyc_13_advection_momentum_centered_2", MDRangePolicy<Kokkos::Rank<3>>(
        koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), 
            FuncAdvMomCen2());
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    parallel_for ("readyc_13_advection_momentum_flux_2", MDRangePolicy<Kokkos::Rank<3>> (
        koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), 
            FuncAdvMomFlu2());
  } else {
    if (mytid == 0) {
      printf ("%s, %d\n", __FILE__, __LINE__);
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  // End advection_momentum(u, v, wka, dlu, dlv, iblock)
#ifdef SMAG
  call smag2(k);
//------------------
#ifdef SMAG_FZ
#else  // SMAG_FZ
#endif // SMAG_FZ
//-----------------
#else // SMAG
  // const int iblock = 0;
  // int block_id = blocks_clinic[iblock];
  // int local_id = iblock + 1;
  // const struct block this_block = CppBlocks::get_block(&block_id, &local_id);
  // const int ib = this_block.ib;
  // const int ie = this_block.ie;
  // const int jb = this_block.jb;
  // const int je = this_block.je;

#ifdef BIHAR
  parallel_for ("readyc_14_hdiffu_del4_1", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, NY_BLOCK, NX_BLOCK}, tile3D), FunctorReadyc14());

  parallel_for ("readyc_15_hdiffu_del4_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, NY_BLOCK, NX_BLOCK}, tile3D), FunctorReadyc15());

  parallel_for ("readyc_16_hdiffu_del4_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, NY_BLOCK, NX_BLOCK}, tile3D), FunctorReadyc16());

  parallel_for ("readyc_17_hdiffu_del4_4", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 2, 2}, koArr3D{KM, JMT-2, IMT-2}, tile3D), FunctorReadyc17());
#else // BIHAR
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      // hdiffu_del2(k, hduk, hdvk, up[iblock][k], vp[iblock][k], this_block)
      parallel_for("readyc_hdiffu_del2_1",
          MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
                  functor_readyc_hdiffu_del2_1());
      parallel_for("readyc_hdiffu_del2_2",
          MDRangePolicy<Kokkos::Rank<2>>
              ({ib-1, jb-1}, {ie, je}), 
                  functor_readyc_hdiffu_del2_2(k, iblock));
      parallel_for("readyc_hdiffu_del2_3",
          MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
                  functor_readyc_hdiffu_del2_3(k));
      // End hdiffu_del2
      parallel_for("readyc_15",
          MDRangePolicy<Kokkos::Rank<2>>
              ({2, 2}, {IMT-2, JMT-2}), 
                  functor_readyc_15(k, iblock));
    }
  }
#endif // BIHAR
#endif // SMAG
  parallel_for ("readyc_18", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc18());

  parallel_for("readyc_19", MDRangePolicy<Kokkos::Rank<3>>
          ({0, 0, 0}, {KM, JMT, IMT}, tile3D), FunctorReadyc19());

  parallel_for ("readyc_20", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorReadyc20());

  parallel_for("readyc_21", MDRangePolicy<Kokkos::Rank<3>>
          ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), FunctorReadyc21());

  parallel_for("readyc_22", MDRangePolicy<Kokkos::Rank<3>>
          ({0, 1, 1}, {KM, JMT-1, IMT-1}, tile3D), FunctorReadyc22());
#ifdef CANUTO
#else  // CANUTO
  if (mytid == 0) {
    printf ("%s, %d\n", __FILE__, __LINE__);
    printf("The false mixing option\n");
  }
  exit(0);
#endif // CANUTO

  return ;
}
//--------------------
//  END READYC
//--------------------
#endif // LICOM_ENABLE_KOKKOS
