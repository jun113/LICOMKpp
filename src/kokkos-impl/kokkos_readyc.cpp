#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_readyc.hpp"

extern "C" void readyc_debug_(double*wk1,double*wk2,double*wk3,double*wk4,double*wk5,double*wk6,double*wk7,
double*wk8,double*wk9,double*wk10,double*wk11,double*wk12,double*wk13);

#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_pconst_mod.h"

void kokkos_readyc () {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  using CppPconstMod::adv_momentum;


#ifdef BCKMEX
  ViewDouble3D v_diff_back("view_diff_back", IMT, JMT, MAX_BLOCKS_CLINIC);
  ViewDouble3D v_diff_back_sh("view_diff_back_sh", IMT, JMT, MAX_BLOCKS_CLINIC);
  ViewDouble3D v_diff_back_nn("view_diff_back_nn", IMT, JMT, MAX_BLOCKS_CLINIC);
#endif // BCKMEX

// using namespace CppDynMod;
// using namespace CppPconstMod;
// using namespace CppTracerMod;
// using namespace CppWorkMod;

// using Kokkos::deep_copy;
// using Kokkos::create_mirror_view;

//   UnManagedViewDouble3D h_v_h0bf (&h0bf[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_h0bl (&h0bl[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_h0 (&h0[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_amld (&amld[0][0][0],
//       1, JMT, IMT); 

//   UnManagedViewDouble4D h_v_up (&up[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_vp (&vp[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_u (&u[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_v (&v[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_akmt (&akmt[0][0][0][0],
//       MAX_BLOCKS_CLINIC, KM, JMT, IMT);
//   UnManagedViewDouble4D h_v_akmu (&akmu[0][0][0][0],
//       MAX_BLOCKS_CLINIC, KM, JMT, IMT);
//   UnManagedViewDouble5D h_v_akt (&akt[0][0][0][0][0],
//       MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
//   UnManagedViewDouble4D h_v_wka (&wka[0][0][0][0],
//       MAX_BLOCKS_CLINIC, KM, JMT, IMT);

//   ViewDouble4D::HostMirror h_v_wp12 = create_mirror_view(*p_v_wp12);
//   ViewDouble4D::HostMirror h_v_wp13 = create_mirror_view(*p_v_wp13);
//   ViewDouble4D::HostMirror h_v_rit = create_mirror_view(*p_v_rit);
//   ViewDouble4D::HostMirror h_v_s2t = create_mirror_view(*p_v_s2t);
//   ViewDouble4D::HostMirror h_v_ws = create_mirror_view(*p_v_ws);
//   ViewDouble3D::HostMirror h_v_work = create_mirror_view(*p_v_work);
//   ViewDouble4D::HostMirror h_v_uk = create_mirror_view(*p_v_uk);
//   ViewDouble4D::HostMirror h_v_vk = create_mirror_view(*p_v_vk);
//   ViewDouble4D::HostMirror h_v_dlu = create_mirror_view(*p_v_dlu);
//   ViewDouble4D::HostMirror h_v_dlv = create_mirror_view(*p_v_dlv);
//   ViewDouble3D::HostMirror h_v_dlub = create_mirror_view(*p_v_dlub);
//   ViewDouble3D::HostMirror h_v_dlvb = create_mirror_view(*p_v_dlvb);
//   ViewDouble3D::HostMirror h_v_su = create_mirror_view(*p_v_su);
//   ViewDouble3D::HostMirror h_v_sv = create_mirror_view(*p_v_sv);

//   for (int j = 0; j < JMT; ++j) {
//     for (int i = 0; i < IMT; ++i) {
//       h_v_su(0,j,i) = h_v_wka(0,0,j,i);
//       h_v_sv(0,j,i) = h_v_wka(0,1,j,i);
//     }
//   }

//   deep_copy(*p_v_h0bf, h_v_h0bf);
//   deep_copy(*p_v_h0, h_v_h0);
//   deep_copy(*p_v_up, h_v_up);
//   deep_copy(*p_v_vp, h_v_vp);
//   deep_copy(*p_v_u, h_v_u);
//   deep_copy(*p_v_v, h_v_v);
//   deep_copy(*p_v_su, h_v_su);
//   deep_copy(*p_v_sv, h_v_sv);

  parallel_for ("readyc_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc1());

  // parallel_for ("readyc_2", MDRangePolicy<Kokkos::Rank<3>> (
  //     koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc2());

  parallel_for ("readyc_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc3());

  parallel_for ("readyc_4", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KMM1, JMT, IMT}, tile3D), FunctorReadyc4());

#ifdef BCKMEX
  // bug no run
  parallel_for("readyc_7", MDRangePolicy<Kokkos::Rank<2>>
      ({0, 0}, {IMT, JMT}), functor_readyc_7(v_diff_back, v_diff_back_sh, v_diff_nh));
  parallel_for("readyc_8", MDRangePolicy<Kokkos::Rank<2>>
      ({0, 0}, {IMT, JMT}), functor_readyc_8(v_diff_back, v_diff_back_sh, v_diff_nh));
#endif // BCKMEX

#if (defined CANUTO) || (defined CANUTO2010)
  parallel_for ("readyc_5", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{1, 1}, koArr2D{JMT-1, IMT-1}, tile2D), FunctorReadyc5());

  parallel_for ("readyc_6", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc6());
#endif // CANUTO

  parallel_for ("readyc_7_upwell_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyc7());

  parallel_for ("readyc_8_upwell_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyc8());

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

#ifdef CANUTO2010
#else  // CANUTO
  if (mytid == 0) {
    printf ("%s, %d\n", __FILE__, __LINE__);
    printf("The false mixing option\n");
  }
  exit(0);
#endif // CANUTO

  // deep_copy(h_v_amld, *p_v_amld);
  // deep_copy(h_v_h0bl, *p_v_h0bl);
  // deep_copy(h_v_h0bf, *p_v_h0bf);
  // deep_copy(h_v_wp12, *p_v_wp12);
  // deep_copy(h_v_wp13, *p_v_wp13);
  // deep_copy(h_v_rit, *p_v_rit);
  // deep_copy(h_v_s2t, *p_v_s2t);
  // deep_copy(h_v_akt, *p_v_akt);
  // deep_copy(h_v_akmt, *p_v_akmt);
  // deep_copy(h_v_akmu, *p_v_akmu);
  // deep_copy(h_v_wka, *p_v_wka);
  // deep_copy(h_v_ws, *p_v_ws);
  // deep_copy(h_v_work, *p_v_work);
  // deep_copy(h_v_uk, *p_v_uk);
  // deep_copy(h_v_vk, *p_v_vk);
  // deep_copy(h_v_dlu, *p_v_dlu);
  // deep_copy(h_v_dlv, *p_v_dlv);
  // deep_copy(h_v_dlub, *p_v_dlub);
  // deep_copy(h_v_dlvb, *p_v_dlvb);

  // readyc_debug_(
  //   h_v_wp12.data(),
  //   h_v_wp13.data(),
  //   h_v_rit.data(),
  //   h_v_s2t.data(),
  //   h_v_wka.data(),
  //   h_v_ws.data(),
  //   h_v_work.data(),
  //   h_v_uk.data(),
  //   h_v_vk.data(),
  //   h_v_dlu.data(),
  //   h_v_dlv.data(),
  //   h_v_dlub.data(),
  //   h_v_dlvb.data()
  // );
  return ;
}
//--------------------
//  END READYC
//--------------------
#endif // LICOM_ENABLE_KOKKOS
