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

#include <array>
#include <vector>

using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::JST;
using CppParamMod::JET;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppDomain::nblocks_clinic;

#ifndef BIHAR
static void hdiffu_del2(
    const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const block &);
static void hdifft_del2(const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const block &);
#else  // BIHAR
// static void hdiffu_del4(const int &,
//     double (&)[NY_BLOCK][NX_BLOCK],
//     double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const block &);
static void hdiffu_del4(const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);
// static void hdifft_del4(const int &, 
//     double (&)[NY_BLOCK][NX_BLOCK],
//     double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const block &);
static void hdifft_del4(const int &, 
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);
#endif // BIHAR

// BAROTR
void cpp_barotr() {
  using CppBlocks::all_blocks;
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
  using CppDomain::blocks_clinic;

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
  double dt2k[JMT][IMT];
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
    // pop_haloupdate_barotr1_(&errorcode);
    // CppPOPHaloMod::pop_halo_update_2dr8(work[0],
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
    CppPOPHaloMod::pop_halo_update_2dr8(&(work[0][0][0]), IMT,
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
      // hdiffu_del4(0, hduk, hdvk, 
      //     ubp[iblock], vbp[iblock], this_block);
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
      const struct block this_block = 
          all_blocks[blocks_clinic[0] - 1];
      hdiffu_del2(0, hduk, hdvk, 
          ubp[iblock], vbp[iblock], this_block);
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
    //-----------------------------
    // pop_haloupdate_barotr2_(&errorcode, &errorcode);
    // CppPOPHaloMod::pop_halo_update_2dr8(wka[0][2],
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    // CppPOPHaloMod::pop_halo_update_2dr8(wka[0][3],
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    CppPOPHaloMod::pop_halo_update_3dr8(&(wka[0][2][0][0]), 2, JMT, IMT,
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    //-----------------------------
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
        // hdifft_del4(1, dt2k, hdtk, h0p[iblock], this_block);
        hdifft_del4(1, dt2k, hdtk, h0p[iblock]);
#else  // BIHAR
        hdifft_del2(1, hdtk, h0p[iblock], this_block);
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
    // pop_haloupdate_barotr1_(&errorcode);
    // CppPOPHaloMod::pop_halo_update_2dr8(work[0],
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
    CppPOPHaloMod::pop_halo_update_2dr8(&(work[0][0][0]), IMT,
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

#ifdef BIHAR
// static void hdiffu_del4(const int &k,
//     double (&hduk)[NY_BLOCK][NX_BLOCK], 
//     double (&hdvk)[NY_BLOCK][NX_BLOCK], 
//     const double (&umixk)[NY_BLOCK][NX_BLOCK], 
//     const double (&vmixk)[NY_BLOCK][NX_BLOCK], 
//     const block &this_block) {

//   using CppConstantMod::C0;
//   using CppGrid::kmu;
//   using CppGrid::uarea;

//   using CppHmixDel4::am;
//   using CppHmixDel4::amf;
//   using CppHmixDel4::duc;
//   using CppHmixDel4::due;
//   using CppHmixDel4::dum;
//   using CppHmixDel4::dun;
//   using CppHmixDel4::dus;
//   using CppHmixDel4::duw;
//   using CppHmixDel4::dmc;
//   using CppHmixDel4::dme;
//   using CppHmixDel4::dmn;
//   using CppHmixDel4::dms;
//   using CppHmixDel4::dmw;

//   std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
//       am_factor(nblocks_clinic);
//   for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
//     for (int j = 0; j < NY_BLOCK; ++j) {
//       for (int i = 0; i < NX_BLOCK; ++i) {
//         am_factor[iblock][j][i] = 1.0;
//       }
//     }
//   }

//   double div_out[NY_BLOCK][NX_BLOCK];
//   div(k, div_out, umixk, vmixk);  

//   double curl[NY_BLOCK][NX_BLOCK];
//   zcurl(k, curl, umixk, vmixk);

//   double gradx1[NY_BLOCK][NX_BLOCK];
//   double grady1[NY_BLOCK][NX_BLOCK];
//   grad(k, gradx1, grady1, curl);

//   double gradx2[NY_BLOCK][NX_BLOCK];
//   double grady2[NY_BLOCK][NX_BLOCK];
//   grad(k, gradx2, grady2, curl);

//   const int bid = 0;
//   const int ib = this_block.ib;
//   const int ie = this_block.ie;
//   const int jb = this_block.jb;
//   const int je = this_block.je;

//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       const double dxdy = pow(sqrt(uarea[bid][j][i]), 5) * 45.0;
//       am_factor[bid][j][i] = sqrt(
//           pow(gradx1[j][i], 2) + pow(gradx2[j][i], 2) 
//           + pow(grady1[j][i], 2) + pow(grady2[j][i], 2)) 
//               * dxdy / fabs(am * amf[bid][j][i]);
//     }
//   }
//   for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
//     for (int j = 0; j < NY_BLOCK; ++j) {
//       for (int i = 0; i < NX_BLOCK; ++i) {
//         am_factor[iblock][j][i] = fmin(40.0, am_factor[iblock][j][i]);
//         am_factor[iblock][j][i] = fmax(1.0,  am_factor[iblock][j][i]);
//       }
//     }
//   }
//   double cc[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       cc[j][i] = duc[bid][j][i] + dum[bid][j][i];
//     }
//   }
//   double d2uk[NY_BLOCK][NX_BLOCK];
//   double d2vk[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       d2uk[j][i] = C0;
//       d2vk[j][i] = C0;
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       d2uk[j][i] = (cc[j][i] * umixk[j  ][i  ] +
//               dun[bid][j][i] * umixk[j-1][i  ] +
//               dus[bid][j][i] * umixk[j+1][i  ] +
//               due[bid][j][i] * umixk[j  ][i+1] +
//               duw[bid][j][i] * umixk[j  ][i-1]) +
//              (dmc[bid][j][i] * vmixk[j  ][i  ] +
//               dmn[bid][j][i] * vmixk[j-1][i  ] +
//               dms[bid][j][i] * vmixk[j+1][i  ] +
//               dme[bid][j][i] * vmixk[j  ][i+1] +
//               dmw[bid][j][i] * vmixk[j  ][i-1]);
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       d2vk[j][i] = (cc[j][i] * vmixk[j  ][i  ] +
//               dun[bid][j][i] * vmixk[j-1][i  ] +
//               dus[bid][j][i] * vmixk[j+1][i  ] +
//               due[bid][j][i] * vmixk[j  ][i+1] +
//               duw[bid][j][i] * vmixk[j  ][i-1]) +
//              (dmc[bid][j][i] * umixk[j  ][i  ] +
//               dmn[bid][j][i] * umixk[j-1][i  ] +
//               dms[bid][j][i] * umixk[j+1][i  ] +
//               dme[bid][j][i] * umixk[j  ][i+1] +
//               dmw[bid][j][i] * umixk[j  ][i-1]);
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       if (k <= kmu[bid][j][i]-1) {
//         d2uk[j][i] = am_factor[bid][j][i] *
//             amf[bid][j][i] * d2uk[j][i];
//         d2vk[j][i] = am_factor[bid][j][i] *
//             amf[bid][j][i] * d2vk[j][i];
//       } else {
//         d2uk[j][i] = C0;
//         d2vk[j][i] = C0;
//       }
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       hduk[j][i] = C0;
//       hdvk[j][i] = C0;
//     }
//   }
//   for (int j = jb-1; j < je; ++j) {
//     for (int i = ib-1; i < ie; ++i) {
//       hduk[j][i] = am * ((cc[j][i] * d2uk[j  ][i  ] +
//                     dun[bid][j][i] * d2uk[j-1][i  ] +
//                     dus[bid][j][i] * d2uk[j+1][i  ] +
//                     due[bid][j][i] * d2uk[j  ][i+1] +
//                     duw[bid][j][i] * d2uk[j  ][i-1]) +
//                    (dmc[bid][j][i] * d2vk[j  ][i  ] +
//                     dmn[bid][j][i] * d2vk[j-1][i  ] +
//                     dms[bid][j][i] * d2vk[j+1][i  ] +
//                     dme[bid][j][i] * d2vk[j  ][i+1] +
//                     dmw[bid][j][i] * d2vk[j  ][i-1]));
//     }
//   }
//   for (int j = jb-1; j < je; ++j) {
//     for (int i = ib-1; i < ie; ++i) {
//       hdvk[j][i] = am * ((cc[j][i] * d2vk[j  ][i  ] +
//                     dun[bid][j][i] * d2vk[j-1][i  ] +
//                     dus[bid][j][i] * d2vk[j+1][i  ] +
//                     due[bid][j][i] * d2vk[j  ][i+1] +
//                     duw[bid][j][i] * d2vk[j  ][i-1]) +
//                    (dmc[bid][j][i] * d2uk[j  ][i  ] +
//                     dmn[bid][j][i] * d2uk[j-1][i  ] +
//                     dms[bid][j][i] * d2uk[j+1][i  ] +
//                     dme[bid][j][i] * d2uk[j  ][i+1] +
//                     dmw[bid][j][i] * d2uk[j  ][i-1]));
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       if (k > kmu[bid][j][i]-1) {
//         hduk[j][i] = C0;
//         hdvk[j][i] = C0;
//       }
//     }
//   }
//   return ;
// }

static void hdiffu_del4(const int &k,
    double (&hduk)[NY_BLOCK][NX_BLOCK], 
    double (&hdvk)[NY_BLOCK][NX_BLOCK], 
    const double (&umixk)[NY_BLOCK][NX_BLOCK], 
    const double (&vmixk)[NY_BLOCK][NX_BLOCK]) {

  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;

  using CppConstantMod::C0;
  using CppGrid::kmu;
  using CppGrid::uarea;

  using CppHmixDel4::am;
  using CppHmixDel4::amf;
  using CppHmixDel4::duc;
  using CppHmixDel4::due;
  using CppHmixDel4::dum;
  using CppHmixDel4::dun;
  using CppHmixDel4::dus;
  using CppHmixDel4::duw;
  using CppHmixDel4::dmc;
  using CppHmixDel4::dme;
  using CppHmixDel4::dmn;
  using CppHmixDel4::dms;
  using CppHmixDel4::dmw;

  std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
      am_factor(nblocks_clinic);
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        am_factor[iblock][j][i] = 1.0;
      }
    }
  }

  double div_out[NY_BLOCK][NX_BLOCK];
  div(k, div_out, umixk, vmixk);  

  double curl[NY_BLOCK][NX_BLOCK];
  zcurl(k, curl, umixk, vmixk);

  double gradx1[NY_BLOCK][NX_BLOCK];
  double grady1[NY_BLOCK][NX_BLOCK];
  grad(k, gradx1, grady1, curl);

  double gradx2[NY_BLOCK][NX_BLOCK];
  double grady2[NY_BLOCK][NX_BLOCK];
  grad(k, gradx2, grady2, curl);

  const int bid = 0;

  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      const double dxdy = pow(sqrt(uarea[bid][j][i]), 5) * 45.0;
      am_factor[bid][j][i] = sqrt(
          pow(gradx1[j][i], 2) + pow(gradx2[j][i], 2) 
          + pow(grady1[j][i], 2) + pow(grady2[j][i], 2)) 
              * dxdy / fabs(am * amf[bid][j][i]);
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        am_factor[iblock][j][i] = fmin(40.0, am_factor[iblock][j][i]);
        am_factor[iblock][j][i] = fmax(1.0,  am_factor[iblock][j][i]);
      }
    }
  }
  double cc[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cc[j][i] = duc[bid][j][i] + dum[bid][j][i];
    }
  }
  double d2uk[NY_BLOCK][NX_BLOCK];
  double d2vk[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      d2uk[j][i] = C0;
      d2vk[j][i] = C0;
    }
  }
  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2uk[j][i] = (cc[j][i] * umixk[j  ][i  ] +
              dun[bid][j][i] * umixk[j-1][i  ] +
              dus[bid][j][i] * umixk[j+1][i  ] +
              due[bid][j][i] * umixk[j  ][i+1] +
              duw[bid][j][i] * umixk[j  ][i-1]) +
             (dmc[bid][j][i] * vmixk[j  ][i  ] +
              dmn[bid][j][i] * vmixk[j-1][i  ] +
              dms[bid][j][i] * vmixk[j+1][i  ] +
              dme[bid][j][i] * vmixk[j  ][i+1] +
              dmw[bid][j][i] * vmixk[j  ][i-1]);
    }
  }
  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2vk[j][i] = (cc[j][i] * vmixk[j  ][i  ] +
              dun[bid][j][i] * vmixk[j-1][i  ] +
              dus[bid][j][i] * vmixk[j+1][i  ] +
              due[bid][j][i] * vmixk[j  ][i+1] +
              duw[bid][j][i] * vmixk[j  ][i-1]) +
             (dmc[bid][j][i] * umixk[j  ][i  ] +
              dmn[bid][j][i] * umixk[j-1][i  ] +
              dms[bid][j][i] * umixk[j+1][i  ] +
              dme[bid][j][i] * umixk[j  ][i+1] +
              dmw[bid][j][i] * umixk[j  ][i-1]);
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k <= kmu[bid][j][i]-1) {
        d2uk[j][i] = am_factor[bid][j][i] *
            amf[bid][j][i] * d2uk[j][i];
        d2vk[j][i] = am_factor[bid][j][i] *
            amf[bid][j][i] * d2vk[j][i];
      } else {
        d2uk[j][i] = C0;
        d2vk[j][i] = C0;
      }
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hduk[j][i] = C0;
      hdvk[j][i] = C0;
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hduk[j][i] = am * ((cc[j][i] * d2uk[j  ][i  ] +
                    dun[bid][j][i] * d2uk[j-1][i  ] +
                    dus[bid][j][i] * d2uk[j+1][i  ] +
                    due[bid][j][i] * d2uk[j  ][i+1] +
                    duw[bid][j][i] * d2uk[j  ][i-1]) +
                   (dmc[bid][j][i] * d2vk[j  ][i  ] +
                    dmn[bid][j][i] * d2vk[j-1][i  ] +
                    dms[bid][j][i] * d2vk[j+1][i  ] +
                    dme[bid][j][i] * d2vk[j  ][i+1] +
                    dmw[bid][j][i] * d2vk[j  ][i-1]));
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdvk[j][i] = am * ((cc[j][i] * d2vk[j  ][i  ] +
                    dun[bid][j][i] * d2vk[j-1][i  ] +
                    dus[bid][j][i] * d2vk[j+1][i  ] +
                    due[bid][j][i] * d2vk[j  ][i+1] +
                    duw[bid][j][i] * d2vk[j  ][i-1]) +
                   (dmc[bid][j][i] * d2uk[j  ][i  ] +
                    dmn[bid][j][i] * d2uk[j-1][i  ] +
                    dms[bid][j][i] * d2uk[j+1][i  ] +
                    dme[bid][j][i] * d2uk[j  ][i+1] +
                    dmw[bid][j][i] * d2uk[j  ][i-1]));
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k > kmu[bid][j][i]-1) {
        hduk[j][i] = C0;
        hdvk[j][i] = C0;
      }
    }
  }
  return ;
}
// static void hdifft_del4(const int &k, 
//     double (&d2tk)[NY_BLOCK][NX_BLOCK],
//     double (&hdtk)[NY_BLOCK][NX_BLOCK],
//     const double (&tmix)[NY_BLOCK][NX_BLOCK],
//     const block &this_block) {

//   using CppGrid::kmt;
//   using CppGrid::kmtn;
//   using CppGrid::kmts;
//   using CppGrid::kmte;
//   using CppGrid::kmtw;
//   using CppHmixDel4::ah;
//   using CppHmixDel4::ahf;
//   using CppHmixDel4::dtn;
//   using CppHmixDel4::dts;
//   using CppHmixDel4::dte;
//   using CppHmixDel4::dtw;
//   using CppConstantMod::C0;

//   const int bid = 0;
//   double cc[NY_BLOCK][NX_BLOCK];
//   double cn[NY_BLOCK][NX_BLOCK];
//   double cs[NY_BLOCK][NX_BLOCK];
//   double ce[NY_BLOCK][NX_BLOCK];
//   double cw[NY_BLOCK][NX_BLOCK];

//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       cn[j][i] = 
//           (k <= kmtn[bid][j][i] && k <= kmt[bid][j][i]) 
//               ? dtn[bid][j][i] : C0;
//       cs[j][i] = 
//           (k <= kmts[bid][j][i] && k <= kmt[bid][j][i]) 
//               ? dts[bid][j][i] : C0;
//       ce[j][i] = 
//           (k <= kmte[bid][j][i] && k <= kmt[bid][j][i]) 
//               ? dte[bid][j][i] : C0;
//       cw[j][i] = 
//           (k <= kmtw[bid][j][i] && k <= kmt[bid][j][i]) 
//               ? dtw[bid][j][i] : C0;
//       cc[j][i] = -(cn[j][i] + cs[j][i] + ce[j][i] + cw[j][i]);
//     }
//   }

//   const int ib = this_block.ib;
//   const int ie = this_block.ie;
//   const int jb = this_block.jb;
//   const int je = this_block.je;
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       d2tk[j][i] = ahf[bid][j][i] * (
//           cc[j][i] * tmix[j  ][i  ]
//         + cn[j][i] * tmix[j-1][i  ]
//         + cs[j][i] * tmix[j+1][i  ]
//         + ce[j][i] * tmix[j  ][i+1]
//         + cw[j][i] * tmix[j  ][i-1]);
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       hdtk[j][i] = C0;
//     }
//   }
//   for (int j = jb-1; j < je; ++j) {
//     for (int i = ib-1; i < ie; ++i) {
//       hdtk[j][i] = ah * (cc[j][i] * d2tk[j  ][i  ]
//                        + cn[j][i] * d2tk[j-1][i  ]
//                        + cs[j][i] * d2tk[j+1][i  ]
//                        + ce[j][i] * d2tk[j  ][i+1]
//                        + cw[j][i] * d2tk[j  ][i-1]);
//     }
//   }
//   return ;
// }
static void hdifft_del4(const int &k, 
    double (&d2tk)[NY_BLOCK][NX_BLOCK],
    double (&hdtk)[NY_BLOCK][NX_BLOCK],
    const double (&tmix)[NY_BLOCK][NX_BLOCK]) {
  
  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;

  using CppGrid::kmt;
  using CppGrid::kmtn;
  using CppGrid::kmts;
  using CppGrid::kmte;
  using CppGrid::kmtw;
  using CppHmixDel4::ah;
  using CppHmixDel4::ahf;
  using CppHmixDel4::dtn;
  using CppHmixDel4::dts;
  using CppHmixDel4::dte;
  using CppHmixDel4::dtw;
  using CppConstantMod::C0;

  const int bid = 0;
  double cc[NY_BLOCK][NX_BLOCK];
  double cn[NY_BLOCK][NX_BLOCK];
  double cs[NY_BLOCK][NX_BLOCK];
  double ce[NY_BLOCK][NX_BLOCK];
  double cw[NY_BLOCK][NX_BLOCK];

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cn[j][i] = 
          (k <= kmtn[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dtn[bid][j][i] : C0;
      cs[j][i] = 
          (k <= kmts[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dts[bid][j][i] : C0;
      ce[j][i] = 
          (k <= kmte[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dte[bid][j][i] : C0;
      cw[j][i] = 
          (k <= kmtw[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dtw[bid][j][i] : C0;
      cc[j][i] = -(cn[j][i] + cs[j][i] + ce[j][i] + cw[j][i]);
    }
  }

  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2tk[j][i] = ahf[bid][j][i] * (
          cc[j][i] * tmix[j  ][i  ]
        + cn[j][i] * tmix[j-1][i  ]
        + cs[j][i] * tmix[j+1][i  ]
        + ce[j][i] * tmix[j  ][i+1]
        + cw[j][i] * tmix[j  ][i-1]);
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hdtk[j][i] = C0;
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdtk[j][i] = ah * (cc[j][i] * d2tk[j  ][i  ]
                       + cn[j][i] * d2tk[j-1][i  ]
                       + cs[j][i] * d2tk[j+1][i  ]
                       + ce[j][i] * d2tk[j  ][i+1]
                       + cw[j][i] * d2tk[j  ][i-1]);
    }
  }
  return ;
}
#else // BIHAR 
static void hdiffu_del2(
    const int &k,
    double (&hduk)[NY_BLOCK][NX_BLOCK],
    double (&hdvk)[NY_BLOCK][NX_BLOCK],
    const double (&umixk)[NY_BLOCK][NX_BLOCK],
    const double (&vmixk)[NY_BLOCK][NX_BLOCK],
    const block &this_block) {

  using CppConstantMod::C0;

  using CppHmixDel2::am;
  using CppHmixDel2::dmc;
  using CppHmixDel2::dme;
  using CppHmixDel2::dmn;
  using CppHmixDel2::dms;
  using CppHmixDel2::dmw;
  using CppHmixDel2::duc;
  using CppHmixDel2::due;
  using CppHmixDel2::dum;
  using CppHmixDel2::dun;
  using CppHmixDel2::dus;
  using CppHmixDel2::duw;

  using CppGrid::kmu;
  using CppPconstMod::viv;

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hduk[j][i] = C0;
      hdvk[j][i] = C0;
    }
  }
  const int bid = 0;
  const int ib = this_block.ib;
  const int ie = this_block.ie;
  const int jb = this_block.jb;
  const int je = this_block.je;
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      double cc = duc[bid][j][i] + dum[bid][j][i];
    
      hduk[j][i] = am * ((cc * umixk[j  ][i  ] +
              dun[bid][j][i] * umixk[j-1][i  ] +
              dus[bid][j][i] * umixk[j+1][i  ] +
              due[bid][j][i] * umixk[j  ][i+1] +
              duw[bid][j][i] * umixk[j  ][i-1]) +
             (dmc[bid][j][i] * vmixk[j  ][i  ] +
              dmn[bid][j][i] * vmixk[j-1][i  ] +
              dms[bid][j][i] * vmixk[j+1][i  ] +
              dme[bid][j][i] * vmixk[j  ][i+1] +
              dmw[bid][j][i] * vmixk[j  ][i-1])) *
                  viv[bid][k][j][i];

      hdvk[j][i] = am * ((cc * vmixk[j  ][i  ] +
              dun[bid][j][i] * vmixk[j-1][i  ] +
              dus[bid][j][i] * vmixk[j+1][i  ] +
              due[bid][j][i] * vmixk[j  ][i+1] +
              due[bid][j][i] * vmixk[j  ][i-1]) +
             (dmc[bid][j][i] * umixk[j  ][i  ] +
              dmn[bid][j][i] * umixk[j-1][i  ] +
              dms[bid][j][i] * umixk[j+1][i  ] +
              dme[bid][j][i] * umixk[j  ][i+1] +
              dmw[bid][j][i] * umixk[j  ][i-1])) *
                  viv[bid][k][j][i];
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k > kmu[bid][j][i]-1) {
        hduk[j][i] = C0;
        hdvk[j][i] = C0;
      }
    }
  }
  return ;
}
static void hdifft_del2(const int &k,
    double (&hdtk)[NY_BLOCK][NX_BLOCK],
    const double (&tmix)[NY_BLOCK][NX_BLOCK],
    const block &this_block) {

  using CppGrid::kmt;
  using CppGrid::kmtn;
  using CppGrid::kmts;
  using CppGrid::kmte;
  using CppGrid::kmtw;
  using CppHmixDel2::ah;
  using CppHmixDel2::ahf;
  using CppHmixDel2::dtn;
  using CppHmixDel2::dts;
  using CppHmixDel2::dte;
  using CppHmixDel2::dtw;
  using CppConstantMod::C0;

  const int bid = 0;

  double cc[NY_BLOCK][NX_BLOCK];
  double cn[NY_BLOCK][NX_BLOCK];
  double cs[NY_BLOCK][NX_BLOCK];
  double ce[NY_BLOCK][NX_BLOCK];
  double cw[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cn[j][i] = 
          (k <= kmtn[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dtn[bid][j][i] : C0;
      cs[j][i] = 
          (k <= kmts[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dts[bid][j][i] : C0;
      ce[j][i] = 
          (k <= kmte[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dte[bid][j][i] : C0;
      cw[j][i] = 
          (k <= kmtw[bid][j][i] && k <= kmt[bid][j][i]) 
              ? dtw[bid][j][i] : C0;
      cc[j][i] = -(cn[j][i] + cs[j][i] + ce[j][i] + cw[j][i]);
    }
  }

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hdtk[j][i] = C0;
    }
  }

  const int ib = this_block.ib;
  const int ie = this_block.ie;
  const int jb = this_block.jb;
  const int je = this_block.je;

  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdtk[j][i] = ah * (cc[j][i] * tmix[j  ][i  ]
                       + cn[j][i] * tmix[j-1][i  ]
                       + cs[j][i] * tmix[j+1][i  ]
                       + ce[j][i] * tmix[j  ][i+1]
                       + cw[j][i] * tmix[j  ][i-1]);
    }
  }
  return ;
}
#endif // BIHAR
#endif // LICOM_ENABLE_FORTRAN
