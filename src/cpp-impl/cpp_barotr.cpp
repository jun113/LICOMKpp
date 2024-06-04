#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_grid.h"
#ifndef BIHAR
#include "../head/cpp_hmix_del2.h"
#else  // BIHAR
#include "../head/cpp_hmix_del4.h"
#endif // BIHAR
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_work_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"


// BAROTR
void cpp_barotr() {
  // using CppBlocks::all_blocks;
  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::JST;
  using CppParamMod::JET;
  using CppDomain::nblocks_clinic;
  using CppParamMod::KM;
  using CppConstantMod::G;
  using CppConstantMod::C0;
  using CppConstantMod::P25;
  using CppPconstMod::dtb;
  using CppPconstMod::ebea;
  using CppPconstMod::ebeb;
  using CppPconstMod::isb;
  using CppPconstMod::vit;
  using CppPconstMod::viv;
  using CppPconstMod::nbb;
  using CppPconstMod::dzph;
  using CppDynMod::dlub;
  using CppDynMod::dlvb;
  using CppDynMod::h0;
  using CppDynMod::h0f;
  using CppDynMod::h0p;
  using CppDynMod::h0bf;
  using CppDynMod::ub;
  using CppDynMod::vb;
  using CppDynMod::ubp;
  using CppDynMod::vbp;
  // using CppDomain::blocks_clinic;

  using CppGrid::fcor;
  using CppWorkMod::pax;
  using CppWorkMod::pay;
  using CppWorkMod::pxb;
  using CppWorkMod::pyb;
  using CppWorkMod::whx;
  using CppWorkMod::why;
  using CppWorkMod::wka;
  using CppWorkMod::wgp;
  using CppWorkMod::work;

#ifdef BIHAR
using CppHmixDel4::hdifft_del4;
using CppHmixDel4::hdiffu_del4;
#else //  BIHAR
using CppHmixDel2::hdifft_del2;
using CppHmixDel2::hdiffu_del2;
#endif // BIHAR

#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BAROTR
#define LICOM_ENABLE_TEST_BAROTR
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_BAROTR
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_BAROTR

#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_start("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wka[iblock][k][j][i] = 0.0;
        }
      }
    }
  }
  // int errorcode;
#ifdef BIHAR
  double dt2k[JMT][IMT];
#endif // BIHAR
  double div_out[JMT][IMT];
  double gradx[JMT][IMT], grady[JMT][IMT];
  double hduk[JMT][IMT], hdvk[JMT][IMT], hdtk[JMT][IMT];

#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR

  for (int nc = 1; nc <= nbb; ++nc) {

#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_start("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      tgrid_to_ugrid(work[iblock], h0[iblock], iblock);
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 1; j < JMT; ++j) {
        for (int i = 0; i < IMT-1; ++i) {
          wka[iblock][0][j][i] = ub[iblock][j][i]
              * (dzph[iblock][j][i] + work[iblock][j][i]);
          wka[iblock][1][j][i] = vb[iblock][j][i]
              * (dzph[iblock][j][i] + work[iblock][j][i]);
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      div(0, div_out, wka[iblock][0], wka[iblock][1]);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          work[iblock][j][i] = vit[iblock][0][j][i]
              * (-1) * div_out[j][i] * P25;
        }
      }
    }
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr calc");
    my_time.testTime_start("barotr halo work");
#endif // LICOM_ENABLE_TEST_BAROTR
    //-----------------------
    CppPOPHaloMod::pop_halo_update (&(work[0][0][0]), IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
    //----------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr halo work");
    my_time.testTime_start("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          h0[iblock][j][i] = h0p[iblock][j][i] 
              + work[iblock][j][i] * dtb;
        }
      }
    }
#ifdef SMAG1
#else // SMAG1
#ifdef BIHAR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // const struct block this_block = 
      //     all_blocks[blocks_clinic[0] - 1];
      hdiffu_del4(0, hduk, hdvk, ubp[iblock], vbp[iblock]);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          wka[iblock][4][j][i] = hduk[j][i];
          wka[iblock][5][j][i] = hdvk[j][i];
        }
      }
    }
#else // BIHAR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // const struct block this_block = 
      //     all_blocks[blocks_clinic[0] - 1];
      hdiffu_del2 (0, hduk, hdvk, ubp[iblock], vbp[iblock]);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          wka[iblock][4][j][i] = hduk[j][i];
          wka[iblock][5][j][i] = hdvk[j][i];
        }
      }
    }
#endif // BIHAR
#endif // SMAG1
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      grad(0, gradx, grady, h0[iblock]);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          const double gstar = (wgp[iblock][j][i] - 1.0) * G;
          wka[iblock][0][j][i] = wka[iblock][4][j][i]
              + gstar * gradx[j][i];
          wka[iblock][1][j][i] = wka[iblock][5][j][i]
              + gstar * grady[j][i];
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      tgrid_to_ugrid(work[iblock], h0[iblock], iblock);
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          wka[iblock][0][j][i] = viv[iblock][0][j][i]
              * (wka[iblock][0][j][i] + dlub[iblock][j][i]
                  - fcor[iblock][j][i] * vbp[iblock][j][i]
                      + pax[iblock][j][i] + pxb[iblock][j][i]
                          - work[iblock][j][i] * whx[iblock][j][i]);

          wka[iblock][1][j][i] = viv[iblock][0][j][i]
              * (wka[iblock][1][j][i] + dlvb[iblock][j][i]
                  + fcor[iblock][j][i] * ubp[iblock][j][i]
                      + pay[iblock][j][i] + pyb[iblock][j][i]
                          - work[iblock][j][i] * why[iblock][j][i]);
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          wka[iblock][2][j][i] = ebea[iblock][j][i] * wka[iblock][0][j][i] 
                               - ebeb[iblock][j][i] * wka[iblock][1][j][i];
          wka[iblock][3][j][i] = ebea[iblock][j][i] * wka[iblock][1][j][i] 
                               + ebeb[iblock][j][i] * wka[iblock][0][j][i];
        }
      }
    }
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr calc");
    my_time.testTime_start("barotr halo wka");
#endif // LICOM_ENABLE_TEST_BAROTR
    CppPOPHaloMod::pop_halo_update(&(wka[0][2][0][0]), 2, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr halo wka");
    my_time.testTime_start("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          ub[iblock][j][i] = ubp[iblock][j][i]
              + wka[iblock][2][j][i] * dtb;
          vb[iblock][j][i] = vbp[iblock][j][i]
              + wka[iblock][3][j][i] * dtb;
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 1; j < JMT; ++j) {
        for (int i = 0; i < IMT-1; ++i) {
          wka[iblock][0][j][i] = ub[iblock][j][i]
              * (dzph[iblock][j][i] + work[iblock][j][i]);
          wka[iblock][1][j][i] = vb[iblock][j][i]
              * (dzph[iblock][j][i] + work[iblock][j][i]);
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // const struct block this_block = 
      //     all_blocks[blocks_clinic[0] - 1];
      div(0, div_out, wka[iblock][0], wka[iblock][1]);
      if (nc % 4 == 0) {
#ifdef BIHAR
        hdifft_del4(0, dt2k, hdtk, h0p[iblock]);
#else  // BIHAR
        hdifft_del2(0, hdtk, h0p[iblock]);
#endif // BIHAR
      } else {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            hdtk[j][i] = C0;
          }
        }
      }
      for (int j = 1; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          work[iblock][j][i] = vit[iblock][0][j][i]
              * (hdtk[j][i] * 1.0 - div_out[j][i]);
        }
      }
    }
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr calc");
    my_time.testTime_start("barotr halo work");
#endif // LICOM_ENABLE_TEST_BAROTR
    //-----------------------
    CppPOPHaloMod::pop_halo_update(&(work[0][0][0]), IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
    //----------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr halo work");
    my_time.testTime_start("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          h0[iblock][j][i] = h0p[iblock][j][i]
              + work[iblock][j][i] * dtb;
        }
      }
    }
    ++isb;
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          ubp[iblock][j][i] = ub[iblock][j][i];
          vbp[iblock][j][i] = vb[iblock][j][i];
          h0p[iblock][j][i] = h0[iblock][j][i];
        }
      }
    }
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = JST-1; j < JET; ++j) {
        for (int i = 0; i < IMT; ++i) {
          h0f[iblock][j][i]  += h0[iblock][j][i];
          h0bf[iblock][j][i] += h0[iblock][j][i];
        }
      }
    }
#ifdef LICOM_ENABLE_TEST_BAROTR
    my_time.testTime_stop("barotr calc");
#endif // LICOM_ENABLE_TEST_BAROTR
  }
  return ;
}
// End BAROTR
//---------------------------------------
#endif // LICOM_ENABLE_FORTRAN
