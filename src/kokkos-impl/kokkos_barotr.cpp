#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_barotr.hpp"

void kokkos_barotr() {
  using CppParamMod::KM;
  using CppPconstMod::isb;
  using CppPconstMod::nbb;
  using CppDomain::blocks_clinic;
  using CppDomain::nblocks_clinic;

  // using CppBlocks::ib;
  // using CppBlocks::ie;
  // using CppBlocks::jb;
  // using CppBlocks::je;

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BAROTR
#define LICOM_ENABLE_TEST_BAROTR
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_BAROTR
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_BAROTR

  // const int iiblock = 0;
  // int block_id = blocks_clinic[iiblock];
  // int local_id = iiblock + 1;
  // const struct block this_block = CppBlocks::get_block(&block_id, &local_id);

  // const int ib = this_block.ib;
  // const int ie = this_block.ie;
  // const int jb = this_block.jb;
  // const int je = this_block.je;

  parallel_for ("barotr_1", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorBarotr1());

  for (int nc = 1; nc <= nbb; ++nc) {

    parallel_for ("barotr_2", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr2());

    parallel_for ("barotr_3", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{1, 0}, koArr2D{JMT, IMT-1}, tile2D), FunctorBarotr3());

    parallel_for ("barotr_4", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr4());

    //---------------------------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_start("barotr haloupdate work");
#endif // LICOM_ENABLE_TEST_BAROTR
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    pop_haloupdate_barotr_1(1, 2);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_work).data(), IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr haloupdate work");
#endif // LICOM_ENABLE_TEST_BAROTR
    //-----------------------------------

    parallel_for ("barotr_5", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr5());

#ifdef SMAG1
#else // SMAG1
#ifdef BIHAR

    parallel_for ("barotr_6_hdiffu_del4_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr6());

    parallel_for ("barotr_7_hdiffu_del4_2", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr7());

    parallel_for ("barotr_8_hdiffu_del4_3", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr8());

    parallel_for ("readyc_9_hdiffu_del4_4", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBarotr9());

#else // BIHAR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // hdiffu_del2(0, hduk, hdvk, ubp[iblock], vbp[iblock], this_block);
      parallel_for("barotr_hdiffu_del2_1",
          MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
                  functor_barotr_hdiffu_del2_1());
      const int iiblock = 0;
      int block_id = blocks_clinic[iiblock];
      int local_id = iiblock + 1;
      const struct block this_block = CppBlocks::get_block(&block_id, &local_id);
      const int ib = this_block.ib;
      const int ie = this_block.ie;
      const int jb = this_block.jb;
      const int je = this_block.je;
      parallel_for("barotr_hdiffu_del2_2",
          MDRangePolicy<Kokkos::Rank<2>>
              ({ib-1, jb-1}, {ie, je}), 
                  functor_barotr_hdiffu_del2_2(0, iblock));
      parallel_for("barotr_hdiffu_del2_3",
          MDRangePolicy<Kokkos::Rank<2>>
              ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
                  functor_barotr_hdiffu_del2_3(0));
      // End hdiffu_del2
      parallel_for("barotr_6", MDRangePolicy<Kokkos::Rank<2>>
          ({2, 2}, {IMT-2, JMT-2}), 
              functor_barotr_6(iblock));
    }
#endif // BIHAR
#endif // SMAG1

    parallel_for ("barotr_11", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr11());

    parallel_for ("barotr_12", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr12());

    parallel_for ("barotr_13", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBarotr13());

    parallel_for ("barotr_14", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBarotr14());

    parallel_for ("barotr_15", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{2, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBarotr15());

    //-----------------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_start("barotr haloupdate wka");
#endif // LICOM_ENABLE_TEST_BAROTR
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    pop_haloupdate_barotr_2(2, 2);
#else
    CppPOPHaloMod::pop_halo_update(&((*p_v_wka).data())[2*JMT*IMT], 2, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#endif
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr haloupdate wka");
#endif // LICOM_ENABLE_TEST_BAROTR
    //-----------------------------

    parallel_for ("barotr_16", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr16());

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // ===========================
      // merge and original versions need
      // div(0, div_out, wka[iblock][0], wka[iblock][1]);
      parallel_for ("barotr_17_div", MDRangePolicy<Kokkos::Rank<2>> (
          koArr2D{0, 0}, koArr2D{NY_BLOCK, NX_BLOCK}, tile2D), FunctorBarotr17(iblock));
      // End div
      // ============================
      if (nc % 4 == 0) {
#ifdef BIHAR
        parallel_for ("barotr_18_hdifft_del4_1", MDRangePolicy<Kokkos::Rank<2>> (
            koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr18());

        parallel_for ("barotr_19_hdifft_del4_2", MDRangePolicy<Kokkos::Rank<2>> (
            koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr19());
#else  // BIHAR
        // hdifft_del2(0, dt2k, hdtk, h0p[iblock], this_block)
        // parallel_for("barotr_hdifft_del2_1", MDRangePolicy<Kokkos::Rank<2>>
        //     ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
        //         functor_barotr_hdifft_del2_1(0));

        // parallel_for("barotr_hdifft_del2_2", MDRangePolicy<Kokkos::Rank<2>>
        //     ({0, 0}, {NX_BLOCK, NY_BLOCK}), 
        //         functor_barotr_hdifft_del2_2());

        // parallel_for("barotr_hdifft_del2_3", MDRangePolicy<Kokkos::Rank<2>>
        //     ({ib-1, jb-1}, {ie, je}), 
        //         functor_barotr_hdifft_del2_3());
        // End hdifft_del2
        parallel_for("barotr_13", MDRangePolicy<Kokkos::Rank<2>>
            ({2, 1}, {IMT-2, JMT-2}), 
                functor_barotr_13(iblock));
#endif // BIHAR
      } else {
        parallel_for("barotr_20", MDRangePolicy<Kokkos::Rank<2>> (
            koArr2D{1, 2}, koArr2D{JMT-2, IMT-2}, tile2D), FunctorBarotr20());
      }
    } // End for iblock
    //---------------------------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_start("barotr haloupdate work");
#endif // LICOM_ENABLE_TEST_BAROTR
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    // CUDA HIP memcpy
    pop_haloupdate_barotr_1(1, 2);
#else
    CppPOPHaloMod::pop_halo_update((*p_v_work).data(), IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr haloupdate work");
#endif // LICOM_ENABLE_TEST_BAROTR
    //-----------------------------------

    parallel_for ("barotr_21", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorBarotr21());

    ++isb;
  }
  // if (CppParamMod::mytid == 0) {
  //   std::cout<<"FunctorBarotr1:  "<<sizeof(FunctorBarotr1)<<std::endl;
  //   std::cout<<"FunctorBarotr2:  "<<sizeof(FunctorBarotr2)<<std::endl;
  //   std::cout<<"FunctorBarotr3:  "<<sizeof(FunctorBarotr3)<<std::endl;
  //   std::cout<<"FunctorBarotr4:  "<<sizeof(FunctorBarotr4)<<std::endl;
  //   std::cout<<"FunctorBarotr5:  "<<sizeof(FunctorBarotr5)<<std::endl;
  //   std::cout<<"FunctorBarotr6:  "<<sizeof(FunctorBarotr6)<<std::endl;
  //   std::cout<<"FunctorBarotr7:  "<<sizeof(FunctorBarotr7)<<std::endl;
  //   std::cout<<"FunctorBarotr8:  "<<sizeof(FunctorBarotr8)<<std::endl;
  //   std::cout<<"FunctorBarotr9:  "<<sizeof(FunctorBarotr9)<<std::endl;
  //   std::cout<<"FunctorBarotr10: "<<sizeof(FunctorBarotr10)<<std::endl;
  //   std::cout<<"FunctorBarotr11: "<<sizeof(FunctorBarotr11)<<std::endl;
  //   std::cout<<"FunctorBarotr12: "<<sizeof(FunctorBarotr12)<<std::endl;
  //   std::cout<<"FunctorBarotr13: "<<sizeof(FunctorBarotr13)<<std::endl;
  //   std::cout<<"FunctorBarotr14: "<<sizeof(FunctorBarotr14)<<std::endl;
  //   std::cout<<"FunctorBarotr15: "<<sizeof(FunctorBarotr15)<<std::endl;
  //   std::cout<<"FunctorBarotr16: "<<sizeof(FunctorBarotr16)<<std::endl;
  //   std::cout<<"FunctorBarotr17: "<<sizeof(FunctorBarotr17)<<std::endl;
  //   std::cout<<"FunctorBarotr18: "<<sizeof(FunctorBarotr18)<<std::endl;
  //   std::cout<<"FunctorBarotr19: "<<sizeof(FunctorBarotr19)<<std::endl;
  //   std::cout<<"FunctorBarotr20: "<<sizeof(FunctorBarotr20)<<std::endl;
  //   std::cout<<"FunctorBarotr21: "<<sizeof(FunctorBarotr21)<<std::endl;
  // }
  return ;
}
// End BAROTR
//---------------------------------------
#endif // LICOM_ENABLE_KOKKOS
