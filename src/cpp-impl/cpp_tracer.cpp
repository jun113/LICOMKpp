#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_buf_mod.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_grid.h"
#ifndef BIHAR
#include "../head/cpp_hmix_del2.h"
#else  // BIHAR
#include "../head/cpp_hmix_del4.h"
#endif // BIHAR
#include "../head/cpp_isopyc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pmix_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_tracer_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include <string>

using CppConstantMod::C0;
using CppConstantMod::P5;
using CppConstantMod::P25;
using CppDomain::nblocks_clinic;
using CppParamMod::mytid;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JST;
using CppParamMod::JET;
using CppParamMod::NTRA;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;

#ifndef ISO
#ifndef SMAG
#ifndef BIHAR
static void hdifft_del2(const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const block &);
#else  // BIHAR
static void hdifft_del4(const int &, 
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);
// static void hdifft_del4(const int &, 
//     double (&)[NY_BLOCK][NX_BLOCK],
//     double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const block &);
#endif // BIHAR
#endif // SMAG
#endif // ISO

static void advection_tracer(const double (&)[KM][JMT][IMT],
    const double (&)[KM][JMT][IMT], const double (&)[KM+1][JMT][IMT],
        const double (&)[KM][JMT][IMT], double (&adv_tt)[KM][JMT][IMT],
            const int &iblock, const int &n);

static void smts(double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], const double &); 

template<typename T> T myMax(const T &val1, const T &val2) {
  return val1 > val2 ? val1 : val2;
}

template<typename T, typename... Args> T myMax(const T &val, const Args &... arg) {
  T result = myMax(arg...);
  return myMax(val, result);
}

template<typename T> T myMin(const T &val1, const T &val2) {
    return val1 < val2 ? val1 : val2;
}

template<typename T, typename... Args> T myMin(const T &val, const Args &... arg) {
  T result = myMin(arg...);
  return myMin(val, result);
}
inline double DEGtoRAD(const double &degree) {
  return (degree * (PI / 180));
}

// TRACER
void cpp_tracer() {

#ifndef ISO
  using CppBlocks   ::all_blocks;
#endif // ISO
#ifdef COUP
  using CppBufMod   ::ifrac;
#endif // COUP
#ifndef ISO
  using CppDomain   ::blocks_clinic;
#endif // ISO
  using CppDynMod   ::ws;
  using CppDynMod   ::h0f;
  using CppDynMod   ::h0l;
  using CppDynMod   ::utf;
  using CppDynMod   ::vtf;
  using CppDynMod   ::utl;
  using CppDynMod   ::vtl;
#ifdef COUP
  using CppForcMod  ::ssf;
#endif // COUP
  using CppForcMod  ::sss;
  using CppForcMod  ::sst;
  using CppForcMod  ::swv;
#ifndef FRC_CORE
  using CppForcMod  ::dqdt;
#endif // FRC_CORE
  using CppForcMod  ::nswv;
  using CppForcMod  ::fresh;
  using CppForcMod  ::seaice;
  using CppForcMod  ::restore;
  using CppGrid     ::kmt;
  using CppGrid     ::tarea;
  using CppGrid     ::area_t;
  using CppGrid     ::horiz_grid_opt;
#ifdef ISO
  // using CppIsopycMod::k3;
  using CppIsopycMod::ahisop;
#endif // ISO
  using CppPconstMod::akt;
  using CppPconstMod::ahv;
  using CppPconstMod::dts;
  using CppPconstMod::ist;
  using CppPconstMod::vit;
  using CppPconstMod::od0;
  using CppPconstMod::aft1;
  using CppPconstMod::aft2;
  using CppPconstMod::od0cp;
  using CppPconstMod::odzp;
  using CppPconstMod::odzt;
  using CppPconstMod::onbc;
  using CppPconstMod::oncc;
  using CppPconstMod::gamma;
  using CppPconstMod::dwndmix;
  using CppPconstMod::adv_tracer;
  using CppPconstMod::simple_assm;
  using CppPconstMod::boundary_restore;
  using CppTracerMod::at;
#ifndef ISO
  using CppTracerMod::dx;
#endif // ISO
  using CppTracerMod::dz;
  using CppTracerMod::atb;
  using CppTracerMod::net;
  using CppTracerMod::tend;
  using CppTracerMod::dt_diff;
  using CppTracerMod::fw_norm2;
  using CppTracerMod::penetrate;
  // using CppWorkMod  ::wkb;
  // using CppWorkMod  ::wkc;
  // using CppWorkMod  ::wkd;
#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_TRACER
#define LICOM_ENABLE_TEST_TRACER
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_TRACER
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_TRACER
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

  double wkb[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double wkc[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double wkd[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#ifdef SOLAR
  using CppPmixMod  ::pen;
#endif // SOLAR
#ifdef SOLARCHLORO
  using CppPmixMod  ::pen_chl;
#endif // SOLARCHLORO

  double stf[MAX_BLOCKS_CLINIC][JMT][IMT];
  double tf[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#ifdef ISO
  double k1[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT];
  double k2[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT];
  double k3[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT];
  double adv_vetiso[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
#endif // ISO

  double aa;
  double c2dtts;
  const std::string str_adv_tracer(adv_tracer);
  if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
    if (ist >= 1) {
      aa = 0.5;
      c2dtts = dts * 2.0;
    } else {
      aa = 0.0;
      c2dtts = dts;
    }
  } else if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
    aa = 0.5;
    c2dtts = dts;
  } else {
    if (mytid == 0) {
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h0f[iblock][j][i] *= onbc;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        stf[iblock][j][i] = aa * h0f[iblock][j][i]
            + (1.0 - aa) * h0l[iblock][j][i];
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          utf[iblock][k][j][i] *= oncc;
          vtf[iblock][k][j][i] *= oncc;
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wkd[iblock][k][j][i] = aa * utf[iblock][k][j][i]
              + (1.0 - aa) * utl[iblock][k][j][i];
          wkb[iblock][k][j][i] = aa * vtf[iblock][k][j][i]
              + (1.0 - aa) * vtl[iblock][k][j][i];
        }
      }
    }
  }
  upwell(wkd, wkb, stf);

#ifdef NODIAG

#ifdef ISO
  isopyc(k1, k2, k3, adv_vetiso);
#endif // ISO

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
#ifdef ISO
          wkc[iblock][k][j][i] = ahv + ahisop[iblock][j][i]
              * k3[iblock][2][j][k][i];
#else  // ISO
          wkc[iblock][k][j][i] = ahv;
#endif // ISO

#ifndef CANUTO
          akt[iblock][0][k][j][i] = wkc[iblock][k][j][i];
          akt[iblock][1][k][j][i] = wkc[iblock][k][j][i];
#endif // CANUTO
        }
      }
    }
  }

  // int errorCode;
  const double aidif = 0.5;
  double vtl_ori[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  const std::string str_horiz_grid_opt(horiz_grid_opt);
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

  for (int n = 0; n < NTRA; ++n) {
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            tf[iblock][k][j][i] = 0.0;
          }
        }
      }
    }

    double adv_tt[KM][JMT][IMT];
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            adv_tt[k][j][i] = 0.0;
          }
        }
      }

      advection_tracer(wkd[iblock], wkb[iblock], ws[iblock], at[iblock][n], 
          adv_tt, iblock, n);

      for (int k = 0; k < KM; ++k) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            tf[iblock][k][j][i] = adv_tt[k][j][i] * vit[iblock][k][j][i];
          }
        }
      }
    }

#ifdef CANUTO
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 1; j < JMT-1; ++j) {
          for (int i = 1; i < IMT-1; ++i) {
            if (akt[iblock][n][0][j][i] < dwndmix) {
              akt[iblock][n][0][j][i] = dwndmix;
            }
#ifdef ISO
            // Note: allocate(k3(imt, 0:km, jmt, 1:3, max_blocks_clinic))
            wkc[iblock][k][j][i] = akt[iblock][n][k][j][i] 
                + ahisop[iblock][j][i] * k3[iblock][2][j][k+1][i];
#else  // ISO
            wkc[iblock][k][j][i] = akt[iblock][n][k][j][i];
#endif // ISO
          }
        }
      }
    }
#endif // CANUTO

#ifdef LPFDIAG
    if (mytid == 0) {
      printf("before ISOFLUX\n");
    }
    call oustsumglobal(tf, outsum);
    if (mytid == 0) {
      printf("outsum: %f\n", outsum);
    }
#endif // LPFDIAG

#ifdef ISO
    isoflux(n, k1, k2, k3, adv_vetiso, tf);
#ifdef LPFDIAG
    if (mytid == 0) {
      printf("after ISOFLUX\n");
    }
    call oustsumglobal(tf, outsum);
    if (mytid == 0) {
      printf("outsum: %f\n", outsum);
    }
#endif // LPFDIAG

#else  // ISO 

#ifdef SMAG
    call SMAG3
#else  // SMAG

#ifdef BIHAR

    double hdtk[JMT][IMT], dt2k[JMT][IMT];
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      // const struct block this_block = 
      //     all_blocks[blocks_clinic[0] - 1];
      for (int k = 0; k < KM; ++k) {
        // hdifft_del4(k, dt2k, hdtk, atb[iblock][n][k+1], this_block);
        hdifft_del4(k, dt2k, hdtk, atb[iblock][n][k+1]);
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            tf[iblock][k][j][i] += hdtk[j][i];
            dx[iblock][n][k][j][i] = hdtk[j][i];
          }
        }
      }
    }

#else  // BIHAR

    double hdtk[JMT][IMT];
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      const struct block this_block = 
          all_blocks[blocks_clinic[0] - 1];
      for (int k = 0; k < KM; ++k) {
        hdifft_del2(k, hdtk, atb[iblock][n][k+1], this_block);
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            tf[iblock][k][j][i] += hdtk[j][i];
            dx[iblock][n][k][j][i] = hdtk[j][i];
          }
        }
      }
    }

#endif // BIHAR
#endif // SMAG
#endif // ISO

#ifdef LPFDIAG
    if (mytid == 0) {
      printf("after1 ISOFLUX\n");
    }
    call oustsumglobal(tf, outsum);
    if (mytid == 0) {
      printf("outsum: %f\n", outsum);
    }
#endif // LPFDIAG

    if (n == 0) {

#ifdef SOLAR
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 1; k < KM-1; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              const double wt1 = swv[iblock][j][i] 
                  * pen[k-1] * vit[iblock][k  ][j][i];
              const double wt2 = swv[iblock][j][i] 
                  * pen[k  ] * vit[iblock][k+1][j][i];

              tf[iblock][k][j][i] += (wt1 - wt2) * odzp[k];

              penetrate[iblock][k][j][i] = (wt1 - wt2) * odzp[k];
            }
          }
        }
      }

      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            const double wt1 = swv[iblock][j][i] 
                * pen[0   ] * vit[iblock][1   ][j][i];
            const double wt2 = swv[iblock][j][i] 
                * pen[KM-2] * vit[iblock][KM-1][j][i];

            tf[iblock][0   ][j][i] -= odzp[0   ] * wt1;
            tf[iblock][KM-1][j][i] += odzp[KM-1] * wt2;

            penetrate[iblock][0   ][j][i] = - wt1 * odzp[0];
            penetrate[iblock][KM-1][j][i] =   wt2 * odzp[KM-1];
          }
        }
      }
#endif // SOLAR

#ifdef SOLARCHLORO
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 1; k < KM-1; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              const double wt1 = swv[iblock][j][i] 
                  * pen_chl[iblock][k-1][j][i] * vit[iblock][k  ][j][i];
              const double wt2 = swv[iblock][j][i] 
                  * pen_chl[iblock][k  ][j][i] * vit[iblock][k+1][j][i];

              tf[iblock][k][j][i] += (wt1 - wt2) * odzp[k];

              penetrate[iblock][k][j][i] = (wt1 - wt2) * odzp[k];
            }
          }
        }
      }

      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            const double wt1 = swv[iblock][j][i] 
                * pen_chl[iblock][0   ][j][i] * vit[iblock][1   ][j][i];
            const double wt2 = swv[iblock][j][i] 
                * pen_chl[iblock][KM-2][j][i] * vit[iblock][KM-1][j][i];

            tf[iblock][0   ][j][i] -= odzp[0   ] * wt1;
            tf[iblock][KM-1][j][i] += odzp[KM-1] * wt2;

            penetrate[iblock][0   ][j][i] = - wt1 * odzp[0];
            penetrate[iblock][KM-1][j][i] =   wt2 * odzp[KM-1];
          }
        }
      }
#endif // SOLARCHLORO
    }

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 1; k < KM-1; ++k) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            const double wt1 = wkc[iblock][k-1][j][i]
                * (atb[iblock][n][k  ][j][i] - atb[iblock][n][k+1][j][i])
                    * odzt[k  ] * vit[iblock][k  ][j][i];
            const double wt2 = wkc[iblock][k  ][j][i]
                * (atb[iblock][n][k+1][j][i] - atb[iblock][n][k+2][j][i])
                    * odzt[k+1] * vit[iblock][k+1][j][i];
            
            tf[iblock][k][j][i] += odzp[k] * (wt1 - wt2) * (1.0 - aidif);

            dz[iblock][n][k][j][i] = odzp[k] * (wt1 - wt2) * (1.0 - aidif);
          }
        }
      }
    }

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            const double wt1 = wkc[iblock][0   ][j][i]
                * (atb[iblock][n][1   ][j][i] - atb[iblock][n][2 ][j][i])
                    * odzt[1   ] * vit[iblock][1   ][j][i];
            const double wt2 = wkc[iblock][KM-2][j][i]
                * (atb[iblock][n][KM-1][j][i] - atb[iblock][n][KM][j][i])
                    * odzt[KM-1] * vit[iblock][KM-1][j][i];
            
            tf[iblock][0   ][j][i] -= odzp[0   ] * wt1 * (1.0 - aidif);
            tf[iblock][KM-1][j][i] += odzp[KM-1] * wt2 * (1.0 - aidif);

            dz[iblock][n][0   ][j][i] = - odzp[0   ] * wt1 * (1.0 - aidif);
            dz[iblock][n][KM-1][j][i] =   odzp[KM-1] * wt2 * (1.0 - aidif);
          }
        }
    }

#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
    if (n == 1) {
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            if (kmt[iblock][j][i] > 0) {
#ifdef COUP
              if (boundary_restore == 2) {
                stf[iblock][j][i] = ssf[iblock][j][i] 
                    / odzp[0] + (gamma * 50.0 * odzp[0]) 
                        * (sss[iblock][j][i] - atb[iblock][1][1][j][i]) 
                            * ifrac[iblock][j][i] 
                    / odzp[0] + (gamma * 30.0 / 365.0 / 4.0 * 50.0 * odzp[0])
                        * (sss[iblock][j][i] - atb[iblock][1][1][j][i]) 
                            / odzp[0];
              } else {
                stf[iblock][j][i] = ssf[iblock][j][i] / odzp[0];
              }
#else  // COUP
#ifdef FRC_CORE
              stf[iblock][j][i] = (fresh[iblock][j][i] * 34.7 * od0 * 1.0e-3
                  + gamma 
                      * (sss[iblock][j][i] - atb[iblock][1][1][j][i]) 
                          * seaice[iblock][j][i] / odzp[0]
                  + gamma * 30.0 / 365.0 / 4.0 * 50.0
                      * (sss[iblock][j][i] - atb[iblock][1][1][j][i])
                          / odzp[0] * (1.0 - seaice[iblock][j][i]));

#else  // FRC_CORE
              stf[iblock][j][i] = gamma 
                  * (sss[iblock][j][i] - atb[iblock][1][1][j][i])
                      / odzp[0];
#endif // FRC_CORE
#endif // COUP
              tf[iblock][0][j][i] += stf[iblock][j][i] * (1.0 - aidif) * odzp[0];

              net[iblock][1][j][i] = stf[iblock][j][i] * odzp[0];
            }
          }
        }
      }
#ifdef SSSNORM
      double err_norm1;
      double err_norm2 = 0.0;

      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 2; j < JMT-2; ++j) {
          for (int i = 2; i < IMT-2; ++i) {
            err_norm2 += tarea[iblock][j][i] * net[iblock][1][j][i]
                * vit[iblock][0][j][i];
          }
        }
      }
#ifdef SPMD
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
    my_time.testTime_start("tracer reduce");
#endif // LICOM_ENABLE_TEST_TRACER
      mpi_reduce_tracer_(&err_norm2, &err_norm1);
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer reduce");
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

      if (0) {
        // TODO
      }

      mpi_bcast_tracer_(&err_norm1);

      err_norm2 = - err_norm1 / area_t;

#else  // SPMD
#endif // SPMD
      fw_norm2 = err_norm2;

      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = JST-1; j < JET; ++j) {
          for (int i = 0; i < IMT; ++i) {
            tf[iblock][0][j][i] += err_norm2 * vit[iblock][0][j][i];
          }
        }
      }
#endif // SSSNORM

#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
    } else {
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            if (kmt[iblock][j][i] > 0) {
#ifdef COUP
              stf[iblock][j][i] = tsf[iblock][j][i];
#else  // COUP
#ifdef FRC_CORE
              stf[iblock][j][i] = ((swv[iblock][j][i] + nswv[iblock][j][i]) 
                  * od0cp + seaice[iblock][j][i] * gamma 
                      * (sst[iblock][j][i] - atb[iblock][0][1][j][i]) / odzp[0]);
#else  // FRC_CORE
              stf[iblock][j][i] = (swv[iblock][j][i] + nswv[iblock][j][i]
                  - dqdt[iblock][j][i] * (sst[iblock][j][i] - atb[iblock][0][1][j][i]))
                      * od0cp;
#endif // FRC_CORE
#endif // COUP
              tf[iblock][0][j][i] += stf[iblock][j][i] * odzp[0] * (1.0 - aidif);

              net[iblock][0][j][i] = stf[iblock][j][i] * odzp[0];
            }
          }
        }
      }
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
    }

#ifdef LICOM_ENABLE_TEST_TRACER
    // my_time.testTime_start("tracer halo net");
#endif // LICOM_ENABLE_TEST_TRACER
    // pop_haloupdate_tracer1_(&errorCode);
#ifdef LICOM_ENABLE_TEST_TRACER
    // my_time.testTime_stop("tracer halo net");
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

    if (simple_assm) {
      if (mytid == 0) {
        printf("into restoring, n = %d\n", n+1);
      }
      /*
      // TODO: 20211112 wjl, restore_at
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 0; k < kvt; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              tf[iblock][k][j][i] += vit[iblock][k][j][i]
                  * (restore_at[iblock][n][k][j][i] - atb[iblock][n][k+1][j][i])
                      * lamda1[k];
              dt_restore[iblock][n][k][j][i] = vit[iblock][k][j][i]
                  * (restore_at[iblock][n][k][j][i] - atb[iblock][n][k+1][j][i])
                      * lamda1[k];

            }
          }
        }
      }
      */
    }
    const double lamda = 1.0 / (15.0 * 86400.0);
    if (boundary_restore == 1) {
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 1; k < KM; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              tf[iblock][k][j][i] += vit[iblock][k][j][i]
                  * (restore[iblock][n][k][j][i] - atb[iblock][n][k+1][j][i])
                      * lamda;
            }
          }
        }
      }

      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 1; k < KM; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              tf[iblock][k][j][i] += vit[iblock][k][j][i]
                  * (restore[iblock][n][k][j][i] - atb[iblock][n][k+1][j][i])
                      * lamda;
            }
          }
        }
      }
    }
    if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 0; k < KM; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              vtl[iblock][k][j][i] = at[iblock][n][k][j][i]
                  + dts * tf[iblock][k][j][i];
            }
          }
        }
      }
    } else if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 0; k < KM; ++k) {
          for (int j = 2; j < JMT-2; ++j) {
            for (int i = 2; i < IMT-2; ++i) {
              vtl[iblock][k][j][i] = at[iblock][n][k][j][i]
                  + c2dtts * tf[iblock][k][j][i];
            }
          }
        }
      }
    } else {
      if (mytid == 0) {
        printf("The false advection option for tracer\n");
      }
      exit(0);
    }

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            vtl_ori[iblock][k][j][i] = vtl[iblock][k][j][i];
          }
        }
      }
    }

    invtrit(vtl, stf, wkc, aidif, c2dtts);

#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
    my_time.testTime_start("tracer halo vtl");
#endif // LICOM_ENABLE_TEST_TRACER
    // pop_haloupdate_tracer2_(&errorCode);
    // CppPOPHaloMod::pop_halo_update_3dr8(vtl[0],
    //     CppDomain::POP_haloClinic_C, 
    //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
    //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
    CppPOPHaloMod::pop_halo_update(&vtl[0][0][0][0],
        KM, JMT, IMT, 
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer halo vtl");
    my_time.testTime_start("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER

    for (int k = 0; k < KM; ++k) {
      if (k == 0) {
        for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              dt_diff[iblock][n][k][j][i] = 
                  (vtl[iblock][k][j][i] - vtl_ori[iblock][k][j][i])
                      / c2dtts * vit[iblock][k][j][i] 
                          - stf[iblock][j][i] * aidif * odzp[k];
            }
          }
        }
      } else {
        for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              dt_diff[iblock][n][k][j][i] = 
                  (vtl[iblock][k][j][i] - vtl_ori[iblock][k][j][i])
                      / c2dtts * vit[iblock][k][j][i];
            }
          }
        }
      }
    }

    const double fil_lat1 = 63.0;
    const double fil_lat2 = 63.0;
    if (str_horiz_grid_opt.find("lat_lon") != str_horiz_grid_opt.npos) {
      if (ist % 180 == 1) {
        smts(vtl, vit, fil_lat2);
      } else {
        if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int j = 0; j < JMT; ++j) {
              for (int i = 0; i < IMT; ++i) {
                vtl[iblock][0][j][i] = vtl[iblock][0][j][i]
                    - at[iblock][n][0][j][i] - net[iblock][n][j][i] * dts;
              }
            }
          }
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int k = 1; k < KM; ++k) {
              for (int j = 0; j < JMT; ++j) {
                for (int i = 0; i < IMT; ++i) {
                  vtl[iblock][k][j][i] -= at[iblock][n][k][j][i];
                }
              }
            }
          }
        } else if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int j = 0; j < JMT; ++j) {
              for (int i = 0; i < IMT; ++i) {
                vtl[iblock][0][j][i] = vtl[iblock][0][j][i]
                    - atb[iblock][n][1][j][i] - net[iblock][n][j][i] * c2dtts;
              }
            }
          }
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int k = 1; k < KM; ++k) {
              for (int j = 0; j < JMT; ++j) {
                for (int i = 0; i < IMT; ++i) {
                  vtl[iblock][k][j][i] -= atb[iblock][n][k+1][j][i];
                }
              }
            }
          }
        } else {
          if (mytid == 0) {
            printf("The false advection option for tracer\n");
          }
          exit(0);
        }

        smts(vtl, vit, fil_lat1);

        if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int j = 0; j < JMT; ++j) {
              for (int i = 0; i < IMT; ++i) {
                vtl[iblock][0][j][i] += at[iblock][n][0][j][i] 
                    + net[iblock][n][j][i] * dts;
              }
            }
          }
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int k = 1; k < KM; ++k) {
              for (int j = 0; j < JMT; ++j) {
                for (int i = 0; i < IMT; ++i) {
                  vtl[iblock][k][j][i] += at[iblock][n][k][j][i];
                }
              }
            }
          }
        } else if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int j = 0; j < JMT; ++j) {
              for (int i = 0; i < IMT; ++i) {
                vtl[iblock][0][j][i] += atb[iblock][n][1][j][i] 
                    + net[iblock][n][j][i] * c2dtts;
              }
            }
          }
          for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
            for (int k = 1; k < KM; ++k) {
              for (int j = 0; j < JMT; ++j) {
                for (int i = 0; i < IMT; ++i) {
                  vtl[iblock][k][j][i] += atb[iblock][n][k+1][j][i];
                }
              }
            }
          }
        } else {
          if (mytid == 0) {
            printf("The false advection option for tracer\n");
          }
          exit(0);
        }
      }
    }

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            tend[iblock][n][k][j][i] = (vtl[iblock][k][j][i] - at[iblock][n][k][j][i])
                / c2dtts * vit[iblock][k][j][i];
          }
        }
      }
    }

    if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 0; k < KM; ++k) {
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              at[iblock][n][k][j][i] = vtl[iblock][k][j][i];
              atb[iblock][n][k+1][j][i] = vtl[iblock][k][j][i];
            }
          }
        }
      }
    } else if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
      if (ist >= 1) {
        for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
          for (int k = 0; k < KM; ++k) {
            for (int j = 0; j < JMT; ++j) {
              for (int i = 0; i < IMT; ++i) {
                atb[iblock][n][k+1][j][i] = aft2 * at[iblock][n][k][j][i]
                    + aft1 * (atb[iblock][n][k+1][j][i] + vtl[iblock][k][j][i]);
              }
            }
          }
        }
      }
      for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
        for (int k = 0; k < KM; ++k) {
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              at[iblock][n][k][j][i] = vtl[iblock][k][j][i];
            }
          }
        }
      }
    } else {
      if (mytid == 0) {
        printf("The false advection option for tracer\n");
      }
      exit(0);
    }
#ifdef LICOM_ENABLE_TEST_TRACER
    my_time.testTime_stop("tracer calc");
#endif // LICOM_ENABLE_TEST_TRACER
  } // End loop: N
#else  // NODIAG
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM+1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          atb[iblock][0][k][j][i] = 12.0;
          atb[iblock][1][k][j][i] = C0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int n = 0; n < NTRA; ++n) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            at[iblock][n][k][j][i] = atb[iblock][n][k+1][j][i];
          }
        }
      }
    }
  }
#endif // NODIAG

  ++ist;
  
  return ;
}

// End TRACER

#ifndef ISO
#ifndef SMAG
#ifdef BIHAR

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

  // const int ib = this_block.ib;
  // const int ie = this_block.ie;
  // const int jb = this_block.jb;
  // const int je = this_block.je;
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

#else // BIHAR 
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
#endif // SMAG
#endif // ISO


static void smts(double (&x)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&z)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double &fil_lat) {

  using CppGrid::tlat;

  const int max_nn = 10;

  int nn[JMT];
  for (int j = 0; j < JMT; ++j) {
    nn[j] = 0;
  }

  for (int j = 2; j < JMT-2; ++j) {
    if (cos(tlat[0][j][0]) <= cos(DEGtoRAD(fil_lat))) {
      nn[j] = std::min(max_nn, 
          abs(static_cast<int>(cos(DEGtoRAD(fil_lat)) / cos(tlat[0][j][0]) * 1.2)));
    }
    if (nn[j] < 0) {
      nn[j] = 0;
    }
  }

  double xs[IMT];
  for (int ncy = 0; ncy < max_nn; ++ncy) {
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int j = 2; j < JMT-2; ++j) {
        if (nn[j] >= ncy) {
          for (int k = 0; k < KM; ++k) {
            for (int i = 0; i < IMT; ++i) {
              xs[i] = x[iblock][k][j][i] * z[iblock][k][j][i];
            }
            for (int i = 2; i < IMT-2; ++i) {
              x[iblock][k][j][i] = (xs[i] * (1.0 - 0.25 * z[iblock][k][j][i-1]
                  - 0.25 * z[iblock][k][j][i+1]) + 0.25 
                      * (xs[i-1] + xs[i+1])) * z[iblock][k][j][i];
            }
          }
        }
      }
    }
  }
  int errorCode;
  pop_haloupdate_smts_(&errorCode);
  return ;
}

static void advection_tracer(const double (&uuu)[KM][JMT][IMT],
    const double (&vvv)[KM][JMT][IMT], const double (&www)[KM+1][JMT][IMT],
        const double (&ttt)[KM][JMT][IMT], double (&adv_tt)[KM][JMT][IMT],
            const int &iblock, const int &mtracer) {

  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::hue;
  using CppGrid::hun;
  using CppGrid::htw;
  using CppGrid::hts;
  using CppGrid::tarea_r;
  using CppPconstMod::dts;
  using CppPconstMod::vit;
  using CppPconstMod::nss;
  using CppPconstMod::odzp;
  using CppPconstMod::odzt;
  using CppPconstMod::adv_tracer;
  using CppTracerMod::ax;
  using CppTracerMod::ay;
  using CppTracerMod::az;

  double u_wface[KM][JMT][IMT], v_sface[KM][JMT][IMT];

  const std::string str_adv_tracer(adv_tracer);
  if (str_adv_tracer.find("centered") != str_adv_tracer.npos ||
      str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          v_sface[k][j][i] = (vvv[k][j][i] + vvv[k][j][i+1]) 
              * hts[iblock][j][i] * P25;
        }
      }
    }
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 0; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j-1][i] + uuu[k][j][i]) 
              * htw[iblock][j][i] * P25;
        }
      }
    }
  }

  if (str_adv_tracer.find("flux") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          v_sface[k][j][i] = (vvv[k][j][i  ] * dxu[iblock][j][i  ]
                            + vvv[k][j][i+1] * dxu[iblock][j][i+1]) * P25;
        }
      }
    }
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 0; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j-1][i] * dyu[iblock][j-1][i]
                            + uuu[k][j  ][i] * dyu[iblock][j  ][i]) * P25; 
        }
      }
    }
  }


  if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
    // k = 0
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {
        adv_tt[0][j][i] = (
            - u_wface[0][j  ][i  ] * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
            - u_wface[0][j  ][i+1] * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ])
            - v_sface[0][j  ][i  ] * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
            - v_sface[0][j-1][i  ] * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]))
                * tarea_r[iblock][j][i];
        
        ax[iblock][mtracer][0][j][i] += (
            - u_wface[0][j  ][i  ] * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
            - u_wface[0][j  ][i+1] * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ]))
                * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        ay[iblock][mtracer][0][j][i] += (
            - v_sface[0][j  ][i  ] * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
            - v_sface[0][j-1][i  ] * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]))
                * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        const double adv_z2 = www[1][j][i] * (ttt[0][j][i] - ttt[1][j][i]);

        adv_tt[0][j][i] -= P5 * odzp[0] * adv_z2;

        az[iblock][mtracer][0][j][i] -= P5 * odzp[0] * adv_z2 
            / static_cast<double>(nss); 
      }
    }

    for (int k = 1; k < KM-1; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_tt[k][j][i] = (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]))
                  * tarea_r[iblock][j][i];
      
          ax[iblock][mtracer][k][j][i] += (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]))
                  * tarea_r[iblock][j][i] / static_cast<double>(nss);

          ay[iblock][mtracer][k][j][i] += (
              - v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]))
                  * tarea_r[iblock][j][i] / static_cast<double>(nss);
               
          const double adv_z1 = www[k][j][i] 
              * (ttt[k-1][j][i] - ttt[k  ][j][i]);
               
          const double adv_z2 = www[k+1][j][i] 
              * (ttt[k  ][j][i] - ttt[k+1][j][i]);

          adv_tt[k][j][i] -= P5 * odzp[k] * (adv_z1 + adv_z2);

          az[iblock][mtracer][k][j][i] -= P5 * odzp[k]
              * (adv_z1 + adv_z2) / static_cast<double>(nss);
        }
      }
    }
    // k = KM - 1
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {
        adv_tt[KM-1][j][i] = (
            - u_wface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
            - u_wface[KM-1][j  ][i+1] 
                * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j-1][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]))
                    * tarea_r[iblock][j][i];
        
        ax[iblock][mtracer][KM-1][j][i] += (
            - u_wface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
            - u_wface[KM-1][j  ][i+1] 
                * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ]))
                    * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        ay[iblock][mtracer][KM-1][j][i] += (
            - v_sface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j-1][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]))
                    * tarea_r[iblock][j][i] / static_cast<double>(nss);

        const double adv_z1 = www[KM-1][j][i] 
            * (ttt[KM-2][j][i] - ttt[KM-1][j][i]);

        adv_tt[KM-1][j][i] -= P5 * odzp[KM-1] * adv_z1;

        az[iblock][mtracer][KM-1][j][i] -= P5 * odzp[KM-1] * adv_z1
            / static_cast<double>(nss); 
      }
    }
  } else if (str_adv_tracer.find("flux") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_tt[k][j][i] = (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] + ttt[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] + ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] + ttt[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] + ttt[k][j  ][i  ]))
                  * tarea_r[iblock][j][i];
          double adv_z1, adv_z2; 
          if (k == 0) {
            adv_z1 = 0.0;
          } else {
            adv_z1 = www[k  ][j][i] * (ttt[k-1][j][i] + ttt[k  ][j][i]) * P5;
          }
          if (k == KM-1) {
            adv_z2 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (ttt[k  ][j][i] + ttt[k+1][j][i]) * P5;
          }
          adv_tt[k][j][i] -= odzp[k] * (adv_z2 - adv_z1);
        }
      }
    }
  } else if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
    double at00[KM][JMT][IMT], atmax[KM][JMT][IMT], atmin[KM][JMT][IMT];

    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          const double adv_x0 = (
              ttt[k][j][i+1] + ttt[k][j][i  ]) 
                  * u_wface[k][j][i+1] * tarea_r[iblock][j][i]
           - (ttt[k][j][i  ] + ttt[k][j][i-1]) 
                  * u_wface[k][j][i  ] * tarea_r[iblock][j][i];


          const double temp_1 = (ttt[k][j+1][i] + ttt[k][j  ][i]) 
              * v_sface[k][j  ][i];
          const double temp_2 = (ttt[k][j  ][i] + ttt[k][j-1][i]) 
              * v_sface[k][j-1][i];
          const double adv_y0 = (temp_1 - temp_2) * tarea_r[iblock][j][i];
          // TODO wjl 20211116
          /*
          const double adv_y0 = (
              (ttt[k][j+1][i] + ttt[k][j  ][i]) * v_sface[k][j  ][i]
            - (ttt[k][j  ][i] + ttt[j][j-1][i]) * v_sface[k][j-1][i])
                * tarea_r[iblock][j][i];
          */

          const double adv_xy1 = - dts * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);

          const double adv_xy2 =   dts * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);

          const double adv_xy3 = - dts * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(v_sface[k][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);

          const double adv_xy4 =   dts * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(v_sface[k][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
             
          const double adv_c1 = - ttt[k][j][i] 
              * (u_wface[k][j  ][i+1] - u_wface[k][j  ][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                   
          const double adv_c2 = - ttt[k][j][i] 
              * (v_sface[k][j  ][i  ] - v_sface[k][j-1][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                    
          double adv_za, adv_zc;
          double adv_zb1, adv_zb2;
          if (k == 0) {
            adv_za = - 0.5 * odzp[0] * www[1][j][i] 
                * (ttt[1][j][i] + ttt[0][j][i]);
                 
            adv_zb1 = 0.0;

            adv_zb2 = 0.5 * odzp[0] * pow(www[1][j][i], 2) * odzt[1]
                * (ttt[0][j][i] - ttt[1][j][i]) * dts;

            adv_zc = odzp[0] * ttt[0][j][i] * www[1][j][i];
          } else if (k == KM-1) {

            adv_za = 0.5 * odzp[KM-1] * www[KM-1][j][i]
                * (ttt[KM-1][j][i] + ttt[KM-2][j][i]);

            adv_zb1 = - 0.5 * odzp[KM-1] * pow(www[KM-1][j][i], 2) * odzt[KM-1]
                * (ttt[KM-2][j][i] - ttt[KM-1][j][i]) * dts;

            adv_zb2 = 0.0;

            adv_zc = - odzp[KM-1] * ttt[KM-1][j][i] * www[KM-1][j][i];
          } else {
            adv_za = 0.5 * odzp[k] * www[k  ][j][i] 
                       * (ttt[k][j][i] + ttt[k-1][j][i])
                   - 0.5 * odzp[k] * www[k+1][j][i] 
                       * (ttt[k][j][i] + ttt[k+1][j][i]);

            adv_zb1 = - 0.5 * odzp[k] * pow(www[k  ][j][i], 2) * odzt[k  ]
                * (ttt[k-1][j][i] - ttt[k  ][j][i]) * dts;

            adv_zb2 =   0.5 * odzp[k] * pow(www[k+1][j][i], 2) * odzt[k+1]
                * (ttt[k  ][j][i] - ttt[k+1][j][i]) * dts;
                
            adv_zc = - odzp[k] * ttt[k][j][i] 
                * (www[k  ][j][i] - www[k+1][j][i]);

          }
          const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
          const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
          const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

          at00[k][j][i] = ttt[k][j][i] + (adv_xx + adv_yy + adv_zz) * dts;
        }
      }
    }

    const double wt1 = - 1.0e10;
    const double wt2 = + 1.0e10;

    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          if (k == 0) {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt2);
          } else if (k == KM-1) {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt2);
          } else {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt1, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt2, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt2);
          }
        }
      }
    }
    // k = 0
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        double adv_xy1;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j  ][i+1] > atmax[0][j  ][i+1] || 
            at00[0][j  ][i+1] < atmin[0][j  ][i+1]) {

          adv_xy1 = - (ttt[0][j  ][i+1] - ttt[0][j  ][i  ])
              * fabs(u_wface[0][j  ][i+1]) * tarea_r[iblock][j][i];
        } else {
          adv_xy1 = - dts * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[0][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
        }

        double adv_xy2;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j  ][i-1] > atmax[0][j  ][i-1] || 
            at00[0][j  ][i-1] < atmin[0][j  ][i-1]) {

          adv_xy2 =   (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
              * fabs(u_wface[0][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy2 =   dts * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[0][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
        }

        double adv_xy3;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j+1][i  ] > atmax[0][j+1][i  ] || 
            at00[0][j+1][i  ] < atmin[0][j+1][i  ]) {

          adv_xy3 = - (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
              * fabs(v_sface[0][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy3 = - dts * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[0][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
        }

        double adv_xy4;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j-1][i  ] > atmax[0][j-1][i  ] || 
            at00[0][j-1][i  ] < atmin[0][j-1][i  ]) {
 
          adv_xy4 =   (ttt[0][j  ][i  ] - ttt[0][j-1][i  ])
              * fabs(v_sface[0][j-1][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy4 =   dts * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[0][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
        }

        const double adv_zb1 = 0.0;
        double adv_zb2;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[1][j  ][i  ] > atmax[1][j  ][i  ] || 
            at00[1][j  ][i  ] < atmin[1][j  ][i  ]) {
 
          adv_zb2 = 0.5 * fabs(www[1][j][i]) * odzp[0]
              * (ttt[0][j][i] - ttt[1][j][i]);
        } else {
          adv_zb2 = 0.5 * odzp[0] * pow(www[1][j][i], 2)
              * odzt[1] * (ttt[0][j][i] - ttt[1][j][i]) * dts;
        }

        const double adv_za = - 0.5 * odzp[0] * www[1][j][i]
            * (ttt[1][j][i] + ttt[0][j][i]);
        const double adv_zc = odzp[0] * ttt[0][j][i] * www[1][j][i];

        const double adv_c1 = - ttt[0][j][i] 
            * (u_wface[0][j  ][i+1] - u_wface[0][j  ][i  ])
                * tarea_r[iblock][j][i] * 2.0; 
        const double adv_c2 = - ttt[0][j][i] 
            * (v_sface[0][j  ][i  ] - v_sface[0][j-1][i  ])
                * tarea_r[iblock][j][i] * 2.0; 

        const double adv_x0 = 
            (ttt[0][j  ][i+1] + ttt[0][j  ][i  ]) * u_wface[0][j  ][i+1] 
                * tarea_r[iblock][j][i]
          - (ttt[0][j  ][i  ] + ttt[0][j  ][i-1]) * u_wface[0][j  ][i  ] 
                * tarea_r[iblock][j][i];
                 
        const double adv_y0 = (
            (ttt[0][j+1][i  ] + ttt[0][j  ][i  ]) * v_sface[0][j  ][i  ]
          - (ttt[0][j  ][i  ] + ttt[0][j-1][i  ]) * v_sface[0][j-1][i  ])
                * tarea_r[iblock][j][i];

        const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

        adv_tt[0][j][i] = adv_xx + adv_yy + adv_zz;
 
        ax[iblock][mtracer][0][j][i] = adv_xx;
        ay[iblock][mtracer][0][j][i] = adv_yy;
        az[iblock][mtracer][0][j][i] = adv_zz;
      }
    }

    for (int k = 1; k < KM-1; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          double adv_xy1;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j  ][i+1] > atmax[k  ][j  ][i+1] || 
              at00[k  ][j  ][i+1] < atmin[k  ][j  ][i+1]) {
         
            adv_xy1 = - (ttt[k][j  ][i+1] - ttt[k][j  ][i  ])
                * fabs(u_wface[k][j  ][i+1]) * tarea_r[iblock][j][i];
          } else {
            adv_xy1 = - dts * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i+1], 2)
                    / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
          }

          double adv_xy2;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j  ][i-1] > atmax[k  ][j  ][i-1] || 
              at00[k  ][j  ][i-1] < atmin[k  ][j  ][i-1]) {
         
            adv_xy2 =   (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
                * fabs(u_wface[k][j  ][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy2 =   dts * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1]) * 2.0
                * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i  ], 2)
                    / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
          }

          double adv_xy3;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j+1][i  ] > atmax[k  ][j+1][i  ] || 
              at00[k  ][j+1][i  ] < atmin[k  ][j+1][i  ]) {
         
            adv_xy3 = - (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
                * fabs(v_sface[k][j  ][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy3 = - dts * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(v_sface[k][j  ][i  ], 2)
                    / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
          }

          double adv_xy4;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j-1][i  ] > atmax[k  ][j-1][i  ] || 
              at00[k  ][j-1][i  ] < atmin[k  ][j-1][i  ]) {
         
            adv_xy4 =   (ttt[k][j  ][i  ] - ttt[k][j-1][i  ])
                * fabs(v_sface[k][j-1][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy4 =   dts * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(v_sface[k][j-1][i  ], 2)
                    / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
          }

          double adv_zb1;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k-1][j  ][i  ] > atmax[k-1][j  ][i  ] || 
              at00[k-1][j  ][i  ] < atmin[k-1][j  ][i  ]) {
  
            adv_zb1 = - 0.5 * fabs(www[k  ][j][i]) * odzp[k]
                * (ttt[k-1][j][i] - ttt[k][j][i]);
          } else {
            adv_zb1 = - 0.5 * odzp[k] * pow(www[k  ][j][i], 2)
                * odzt[k  ] * (ttt[k-1][j][i] - ttt[k  ][j][i]) * dts;
          }

          double adv_zb2;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k+1][j  ][i  ] > atmax[k+1][j  ][i  ] || 
              at00[k+1][j  ][i  ] < atmin[k+1][j  ][i  ]) {
  
            adv_zb2 =   0.5 * fabs(www[k+1][j][i]) * odzp[k]
                * (ttt[k  ][j][i] - ttt[k+1][j][i]);
          } else {
            adv_zb2 =   0.5 * odzp[k] * pow(www[k+1][j][i], 2)
                * odzt[k+1] * (ttt[k  ][j][i] - ttt[k+1][j][i]) * dts;
          }

          const double adv_c1 = - ttt[k][j][i] 
              * (u_wface[k][j  ][i+1] - u_wface[k][j  ][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                   
          const double adv_c2 = - ttt[k][j][i] 
              * (v_sface[k][j  ][i  ] - v_sface[k][j-1][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;

          const double adv_za = 
              0.5 * odzp[k] * www[k  ][j][i] * (ttt[k][j][i] + ttt[k-1][j][i])
            - 0.5 * odzp[k] * www[k+1][j][i] * (ttt[k][j][i] + ttt[k+1][j][i]);

          const double adv_zc = - odzp[k] * ttt[k][j][i]
              * (www[k  ][j][i] - www[k+1][j][i]);

          const double adv_x0 = 
              (ttt[k][j  ][i+1] + ttt[k][j  ][i  ]) * u_wface[k][j  ][i+1] 
                  * tarea_r[iblock][j][i]
            - (ttt[k][j  ][i  ] + ttt[k][j  ][i-1]) * u_wface[k][j  ][i  ] 
                  * tarea_r[iblock][j][i];
                   
          const double adv_y0 = (
              (ttt[k][j+1][i  ] + ttt[k][j  ][i  ]) * v_sface[k][j  ][i  ]
            - (ttt[k][j  ][i  ] + ttt[k][j-1][i  ]) * v_sface[k][j-1][i  ])
                  * tarea_r[iblock][j][i];

          const double adv_xx = - (adv_x0 + adv_xy1 +adv_xy2 + adv_c1);
          const double adv_yy = - (adv_y0 + adv_xy3 +adv_xy4 + adv_c2);
          const double adv_zz = - (adv_za + adv_zb1 +adv_zb2 + adv_zc);

          adv_tt[k][j][i] = adv_xx + adv_yy + adv_zz;

          ax[iblock][mtracer][k][j][i] = adv_xx;
          ay[iblock][mtracer][k][j][i] = adv_yy;
          az[iblock][mtracer][k][j][i] = adv_zz;
        }
      }
    }

    // k = KM - 1
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        double adv_xy1;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i+1] > atmax[KM-1][j  ][i+1] || 
            at00[KM-1][j  ][i+1] < atmin[KM-1][j  ][i+1]) {

          adv_xy1 = - (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ])
              * fabs(u_wface[KM-1][j  ][i+1]) * tarea_r[iblock][j][i];
        } else {
          adv_xy1 = - dts * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[KM-1][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
        }

        double adv_xy2;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i-1] > atmax[KM-1][j  ][i-1] || 
            at00[KM-1][j  ][i-1] < atmin[KM-1][j  ][i-1]) {

          adv_xy2 =   (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
              * fabs(u_wface[KM-1][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy2 =   dts * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[KM-1][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
        }

        double adv_xy3;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j+1][i  ] > atmax[KM-1][j+1][i  ] || 
            at00[KM-1][j+1][i  ] < atmin[KM-1][j+1][i  ]) {

          adv_xy3 = - (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
              * fabs(v_sface[KM-1][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy3 = - dts * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[KM-1][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
        }

        double adv_xy4;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j-1][i  ] > atmax[KM-1][j-1][i  ] || 
            at00[KM-1][j-1][i  ] < atmin[KM-1][j-1][i  ]) {

          adv_xy4 =   (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ])
              * fabs(v_sface[KM-1][j-1][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy4 =   dts * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[KM-1][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
        }
        double adv_zb1;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-2][j  ][i  ] > atmax[KM-2][j  ][i  ] || 
            at00[KM-2][j  ][i  ] < atmin[KM-2][j  ][i  ]) {
 
          adv_zb1 = - 0.5 * fabs(www[KM-1][j][i]) * odzp[KM-1]
              * (ttt[KM-2][j][i] - ttt[KM-1][j][i]);
        } else {
          adv_zb1 = - 0.5 * odzp[KM-1] * pow(www[KM-1][j][i], 2)
              * odzt[KM-1] * (ttt[KM-2][j][i] - ttt[KM-1][j][i]) * dts;
        }

        const double adv_zb2 = 0.0;

        const double adv_c1 = - ttt[KM-1][j][i] * 
            (u_wface[KM-1][j  ][i+1] - u_wface[KM-1][j  ][i  ])
                * tarea_r[iblock][j][i] * 2.0;
                 
        const double adv_c2 = - ttt[KM-1][j][i] * 
            (v_sface[KM-1][j  ][i  ] - v_sface[KM-1][j-1][i  ])
                * tarea_r[iblock][j][i] * 2.0;

        const double adv_za = 0.5 * odzp[KM-1] * www[KM-1][j][i]
            * (ttt[KM-1][j][i] + ttt[KM-2][j][i]);

        const double adv_zc = - odzp[KM-1] * ttt[KM-1][j][i] * www[KM-1][j][i];

        const double adv_x0 = 
            (ttt[KM-1][j  ][i+1] + ttt[KM-1][j  ][i  ]) 
                * u_wface[KM-1][j  ][i+1] * tarea_r[iblock][j][i]
          - (ttt[KM-1][j  ][i  ] + ttt[KM-1][j  ][i-1])
                * u_wface[KM-1][j  ][i  ] * tarea_r[iblock][j][i];
                 
        const double adv_y0 = (
            (ttt[KM-1][j+1][i  ] + ttt[KM-1][j  ][i  ]) 
                * v_sface[KM-1][j  ][i  ]
          - (ttt[KM-1][j  ][i  ] + ttt[KM-1][j-1][i  ])
                * v_sface[KM-1][j-1][i  ]) * tarea_r[iblock][j][i];

        const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

        adv_tt[KM-1][j][i] = adv_xx + adv_yy + adv_zz;

        ax[iblock][mtracer][KM-1][j][i] = adv_xx;
        ay[iblock][mtracer][KM-1][j][i] = adv_yy;
        az[iblock][mtracer][KM-1][j][i] = adv_zz;
      }
    }
  }
  return ;
}
#endif // LICOM_ENABLE_FORTRAN
