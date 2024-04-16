#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_canuto_mod.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_hmix_del2.h"
#include "../head/cpp_hmix_del4.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pmix_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/fortran_blocks.h"
#include "../head/fortran_canuto_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using CppDomain  ::nblocks_clinic;
using CppParamMod::mytid;
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::KMM1;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JST;
using CppParamMod::JET;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

static void advection_momentum(
    const double (&)[KM][JMT][IMT],
    const double (&)[KM][JMT][IMT],
    const double (&)[KM][JMT][IMT],
    double (&)[KM][JMT][IMT],
    double (&)[KM][JMT][IMT],
    const int &);
#ifndef BIHAR
static void hdiffu_del2(
    const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const block &);
#else  // BIHAR
static void hdiffu_del4(
    const int &,
    double (&)[NY_BLOCK][NX_BLOCK],
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);
// static void hdiffu_del4(
//     const int &,
//     double (&)[NY_BLOCK][NX_BLOCK],
//     double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const double (&)[NY_BLOCK][NX_BLOCK],
//     const block &);
#endif // BIHAR

void cpp_readyc() {

  using CppConstantMod::C0;

  using CppDomain::    blocks_clinic;
  using CppDynMod::    dlu;
  using CppDynMod::    dlv;
  using CppDynMod::    dlub;
  using CppDynMod::    dlvb;
  using CppDynMod::    h0;
  using CppDynMod::    h0bl;
  using CppDynMod::    h0bf;
  using CppDynMod::    u;
  using CppDynMod::    v;
  using CppDynMod::    up;
  using CppDynMod::    vp;
  using CppDynMod::    ws;
  using CppDynMod::    bbcx;
  using CppDynMod::    bbcy;
  using CppDynMod::    sbcx;
  using CppDynMod::    sbcy;
  using CppForcMod::   buoysol;
  using CppForcMod::   buoytur;
  using CppForcMod::   su;
  using CppForcMod::   sv;
  using CppForcMod::   ustar;
  using CppForcMod::   wave_dis;
  // using CppGrid::      au0;
  // using CppGrid::      aus;
  // using CppGrid::      auw;
  // using CppGrid::      ausw;
  using CppGrid::      kmt;
  using CppGrid::      kmu;
  using CppGrid::      fcort;
  using CppPconstMod:: akmu;
  using CppPconstMod:: akt;
  using CppPconstMod:: akmt;
  using CppPconstMod:: ak_tide;
  using CppPconstMod:: c0f;
  using CppPconstMod:: cag;
  using CppPconstMod:: sag;
  using CppPconstMod:: dzp;
  using CppPconstMod:: dfricmx;
  using CppPconstMod:: dwndmix;
  using CppPconstMod:: fztidal;
  using CppPconstMod:: fz_tide;
  using CppPconstMod:: ncc;
  using CppPconstMod:: richardson;
  using CppPconstMod:: viv;
  using CppPconstMod:: vit;
  using CppPconstMod:: od0;
  using CppPconstMod:: odzt;
  using CppPconstMod:: odzp;
  using CppPconstMod:: ohbu;
  using CppPconstMod:: snlat;
  using CppPconstMod:: wp3_tidal;
  using CppPconstMod:: zkp;
  using CppPmixMod::   ric;
  // using CppPmixMod::   riu;
  using CppPmixMod::   rit;
  using CppPmixMod::   rict;
  // using CppPmixMod::   rict_ref;
  using CppPmixMod::   ridt;
  using CppPmixMod::   ricdttms;
  using CppPmixMod::   s2t;
  using CppPmixMod::   s2u;
  using CppTracerMod:: amld;
  using CppTracerMod:: at;
  using CppTracerMod:: pdensity;
  using CppWorkMod::   wka;

  using CppPconstMod::BACK_TIDALMIXING;
  using CppPconstMod::MAX_TIDALMIXING;
  using CppPconstMod::MIXING_EF;
  using CppPconstMod::LOCAL_MIXING_FRACTION;

  double riv1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double riv2[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = JST - 1; j < JET; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h0bl[iblock][j][i] = h0bf[iblock][j][i];
        h0bf[iblock][j][i] = h0  [iblock][j][i];
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          s2t[iblock][k][j][i]  = C0;
          ridt[iblock][k][j][i] = C0;
        }
      }
    }
  }
  // for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
  //   for (int k = 0; k < KM+1; ++k) {
  //     for (int j = 0; j < JMT; ++j) {
  //       for (int i = 0; i < IMT; ++i) {
  //         riu[iblock][k][j][i] = C0;
  //       }
  //     }
  //   }
  // }
  double wp12[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double wp13[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wp12[iblock][k][j][i] = C0;
          wp13[iblock][k][j][i] = C0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      ugrid_to_tgrid(wp12[iblock][k], up[iblock][k], iblock, k);
      ugrid_to_tgrid(wp13[iblock][k], vp[iblock][k], iblock, k);
    }
  }

  const double epsln = 1.0e-25;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 1; j < JMT; ++j) {
        for (int i = 0; i < IMT-1; ++i) {
          riv1[iblock][k][j][i] = 
              wp12[iblock][k  ][j][i] * vit[iblock][k  ][j][i] - 
              wp12[iblock][k+1][j][i] * vit[iblock][k+1][j][i];
          riv2[iblock][k][j][i] = 
              wp13[iblock][k  ][j][i] * vit[iblock][k  ][j][i] - 
              wp13[iblock][k+1][j][i] * vit[iblock][k+1][j][i];

          s2t[iblock][k][j][i] = vit[iblock][k+1][j][i] * 
              (riv1[iblock][k][j][i] * riv1[iblock][k][j][i] + 
               riv2[iblock][k][j][i] * riv2[iblock][k][j][i]) * 
               odzt[k+1] * odzt[k+1];
#ifdef CANUTO
          rit[iblock][k][j][i] = vit[iblock][k+1][j][i] * rict[iblock][k][j][i] / 
              (s2t[iblock][k][j][i] + epsln);
#else
          // TODO wjl 2021/09/30,  rit[iblock][iblock][k][j][i] ?
          rit[iblock][k][j][i] = rit[iblock][iblock][k][j][i] + vit[iblock][iblock][k+1][j][i] * 
              rict[iblock][iblock][k][j][i] / (s2t[iblock][iblock][k][j][i] + epsln);
#endif // CANUTO
        }
      }
    }
  }
  
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          riv1[iblock][k][j][i] = up[iblock][k][j][i] - up[iblock][k+1][j][i];
          riv2[iblock][k][j][i] = vp[iblock][k][j][i] - vp[iblock][k+1][j][i];
          s2u[iblock][k][j][i] = viv[iblock][k+1][j][i] * 
              (riv1[iblock][k][j][i] * riv1[iblock][k][j][i] +
               riv2[iblock][k][j][i] * riv2[iblock][k][j][i]) *
               odzt[k+1] * odzt[k+1];
          // riu[iblock][k+1][j][i] = viv[iblock][k+1][j][i] * 
          //     ric[iblock][k][j][i]/(s2u[iblock][k][j][i] + epsln);
        }
      }
    }
  }
#ifdef BCKMEX
  double diff_back[MAX_BLOCKS_CLINIC][JMT][IMT];
  double diff_back_sh[MAX_BLOCKS_CLINIC][JMT][IMT];
  double diff_back_nh[MAX_BLOCKS_CLINIC][JMT][IMT];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT;  ++j) {
      for (int i = 0; i < IMT; ++i) {
        diff_back[iblock][j][i]    = 0.0;
        diff_back_sh[iblock][j][i] = 0.0;
        diff_back_nh[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        if (tlat[iblock][j][i] < 0.0) {
          diff_back_sh[iblock][j][i] = diff_back_coef_max * 
              exp(-pow(0.4e0 * (tlat[iblock][j][i] / DEGtoRAD + 28.9), 2); 
        } else {
          diff_back_nh[iblock][j][i] = diff_back_coef_max * 
              exp(-pow(0.4e0 * (tlat[iblock][j][i] / DEGtoRAD - 28.9), 2);
        }
        diff_back[iblock][j][i] = diff_back_eq + 
            diff_back_sh[iblock][j][i] + diff_back_nh[iblock][j][i];

        if (tlat[iblock][j][i] < -10.0 * DEGtoRAD) {
          diff_back[iblock][j][i] += diff_back_coef;
        } else if (fabs(tlat[iblock][j][i]) <= 10.0 * DEGtoRAD) {
          diff_back[iblock][j][i] += diff_back_coef *
              pow((fabs(tlat[iblock][j][i] / DEGtoRAD) / 10.0), 2);
        } else {
          diff_back[iblock][j][i] += diff_back_coef;
        }
      }
    }
  }
#endif // BCKMEX

#ifdef CANUTO
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT;  ++j) {
      for (int i = 0; i < IMT; ++i) {
        amld[iblock][j][i] = C0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT;  ++j) {
        for (int i = 0; i < IMT; ++i) {
          akmt[iblock][k][j][i] = C0;
          akmu[iblock][k][j][i] = C0;
        }
      }
    }
  }
  double akm_back[KM - 1];
  double akt_back[KM - 1];
  double aks_back[KM - 1];
  for (int k = 0; k < KM-1; ++k) {
    akm_back[k] = C0;
    akt_back[k] = C0;
    aks_back[k] = C0;
  }
  double wk1[KM-1], wk2[KM-1], wk3[KM-1];
  double wp1[KM], wp2[KM], wp3[KM];
  double wp4[KM], wp5[KM], wp6[KM];
  double wp7[KM], wp8[KM];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        if (vit[iblock][0][j][i] > 0.5) {
          for (int k = 0; k < KM; ++k) {
            wp1[k] = 0.0;
            wp2[k] = 0.0;
            wp3[k] = 0.0;
            wp4[k] = 0.0;
            wp5[k] = 0.0;
            wp6[k] = 0.0;
            wp7[k] = 0.0;
            wp8[k] = 0.0;
          }
          for (int k = 0; k < kmt[iblock][j][i] - 1; ++k) {
            wp8[k] = - vit[iblock][k+1][j][i] * zkp[k+1] * 1.0e+2;
          }
          for (int k = 0; k < kmt[iblock][j][i] - 1; ++k) {
            wp1[k] = vit[iblock][k+1][j][i] * 
                (at[iblock][0][k][j][i] - (at[iblock][0][k][j][i] - at[iblock][0][k+1][j][i]) * 
                 dzp[k] / (dzp[k] + dzp[k+1]));
            wp2[k] = vit[iblock][k+1][j][i] *
                (at[iblock][1][k][j][i] - (at[iblock][1][k][j][i] - at[iblock][1][k+1][j][i]) *
                dzp[k] / (dzp[k] + dzp[k+1]));
            wp4[k] =  vit[iblock][k+1][j][i] * rit[iblock][k][j][i];
            wp5[k] =  vit[iblock][k+1][j][i] * ricdttms[iblock][k][j][i];
            wp6[k] =  vit[iblock][k+1][j][i] * s2t[iblock][k][j][i];
            wp7[k] =  vit[iblock][k+1][j][i] * rict[iblock][k][j][i];
            /*
            // wjl 20211020
            double tq =  wp1[k] - to[0];
            double sq = (wp2[k] - 35.0) * 1.0e-3 - so[0];
            */

            wp3[k] =  vit[iblock][k+1][j][i] * (pdensity[iblock][k][j][i] + 
                (pdensity[iblock][k+1][j][i] -  pdensity[iblock][k][j][i]) * 
                 dzp[k] / (dzp[k] + dzp[k+1])) * 1.e-3;
          }
          double wp9  = vit[iblock][0][j][i] * ustar  [iblock][j][i] * 1.0e+2;
          double wp10 = vit[iblock][0][j][i] * buoytur[iblock][j][i] * 1.0e+4;
          double wp11 = vit[iblock][0][j][i] * buoysol[iblock][j][i] * 1.0e+4;


#ifdef CANUTOMIXOUT
          wp10_canuto[iblock][j][i] = wp10;
          wp11_canuto[iblock][j][i] = wp11;
#endif // CANUTOMIXOUT
          int iwk  = kmt[iblock][j][i] - 1;
          int iwk1 = std::max(1, kmt[iblock][j][i]-1);
#ifdef CANUTO2010
          double tau_mag = ustar[iblock][j][i] * ustar[iblock][j][i] / od0;
          if(mytid == 0) write(*,*) "into canuto2010"
          // Fortran
          call canuto_2010_interface(wk1, wk2, wk3, wk4, 
              amld(i,j,iblock) ,tau_mag, wp1, wp2, wp3, wp4, wp5, 
              wp7, wp6, ulat(i,j,iblock)/degtorad , wp8, 
              kmt(i,j,iblock), i, j, iblock, isc)
#endif // CANUTO2010
//---------------------------------
          double mldtmp;
          double dfricmx_temp = dfricmx * 1.0e+4;
          double dwndmix_temp = dwndmix * 1.0e+4;
          int nmax        = KM - 1;
          int isurfuse    = 1;
          int ifextermld  = 0;
          int ifoutput    = 0;

          turb_2(wp8, wp1, wp2, wp3, wp4, wp5, wp6, 
              dfricmx_temp, dwndmix_temp, akm_back, akt_back, aks_back,
              wp7, wp9, wp10, wp11, fcort[iblock][j][i], mldtmp, 
              wk1, wk2, wk3, iwk, iwk1, nmax, isurfuse, 
              ifextermld, ifoutput, i, j);

//----------------------------------
          amld[iblock][j][i] = mldtmp * 1.0e-2;
#ifdef TIDEMIX
          for (int k = 0; k < KM; ++k) {
            ak_tide[iblock][k][j][i] = 0.0;
          }
          for (int k = 0; k < kmt[iblock][j][i] - 1; ++k) {
            ak_tide[iblock][k][j][i] = BACK_TIDALMIXING + MIXING_EF * 
                LOCAL_MIXING_FRACTION * wave_dis[iblock][j][i] * 
                fz_tide[iblock][k][j][i] /
                (fmax(rict[iblock][k][j][i], 1.0e-8) * wp3[k] * 1000.0);

            ak_tide[iblock][k][j][i] = fmin(ak_tide[iblock][k][j][i],
                MAX_TIDALMIXING);

            richardson[iblock][k][j][i] = rict[iblock][k][j][i];
            fztidal[iblock][k][j][i]    = fz_tide[iblock][k][j][i];
            wp3_tidal[iblock][k][j][i]  = wp3[k];
          }
#ifdef CANUTOMIXOUT
          for (int k = 0; k < kmt[iblock][j][i] - 1; ++k) {
            wp1_canuto[iblock][k][j][i]  = wp1[k];
            wp2_canuto[iblock][k][j][i]  = wp2[k];
            wp3_canuto[iblock][k][j][i]  = wp3[k];
            wp4_canuto[iblock][k][j][i]  = wp4[k];
            wp5_canuto[iblock][k][j][i]  = wp5[k];
            wp6_canuto[iblock][k][j][i]  = wp7[k];
            wp8_canuto[iblock][k][j][i]  = wp8[k];
            wp12_canuto[iblock][k][j][i] = wp12[k];
            wp13_canuto[iblock][k][j][i] = wp13[k];
            wk4_canuto[iblock][k][j][i]  = wk4[k];
          }
#endif // CANUTOMIXOUT
          for (int k = kmt[iblock][j][i] - 3; k >= 0; --k) {
            ak_tide[iblock][k][j][i] = fmin(ak_tide[iblock][k][j][i],
                ak_tide[iblock][k+1][j][i]);
          }
#endif // TIDEMIX
          for (int k = 0; k < KM - 1; ++k) {
#ifdef TIDEMIX
            akmt[iblock][k][j][i] = wk1[k] * 1.0e-4 +
                ak_tide[iblock][k][j][i] * 5.0;
            akt[iblock][0][k][j][i] += (wk2[k] * 1.0e-4 + 
                ak_tide[iblock][k][j][i]) / 
                    static_cast<float>(ncc);
            akt[iblock][1][k][j][i] += (wk3[k] * 1.0e-4 + 
                ak_tide[iblock][k][j][i]) / 
                    static_cast<float>(ncc);
#ifdef CANUTOMIXOUT
            wk1_canuto[iblock][k][j][i] = wk1[k];
            wk2_canuto[iblock][k][j][i] = wk2[k];
            wk3_canuto[iblock][k][j][i] = wk3[k];
#endif // CANUTOMIXOUT
#ifdef BCKMEX
            akmt[iblock][k][j][i] += diff_back[iblock][j][i] * 
                10.0 * 1.0e-4;
            akt[iblock][0][k][j][i] += diff_back[iblock][j][i] /
                    static_cast<float>(ncc) * 1.0e-4;
            akt[iblock][1][k][j][i] += diff_back[iblock][j][i] /
                    static_cast<float>(ncc) * 1.0e-4;
#endif // BCKMEX
#else // TIDEMIX
            akmt[iblock][k][j][i] = wk1[k] * 1.0e-4;
            akt[iblock][0][k][j][i] += wk2[k] * 1.0e-4 / 
                static_cast<float>(ncc);
            akt[iblock][1][k][j][i] += wk3[k] * 1.0e-4 / 
                static_cast<float>(ncc);
#ifdef BCKMEX
            akmt[iblock][k][j][i]   += diff_back[iblock][j][i] * 
                10.0 * 1.0e-4;
            akt[iblock][0][k][j][i] += diff_back[iblock][j][i] / 
                static_cast<float>(ncc) * 1.0e-4;
            akt[iblock][1][k][j][i] += diff_back[iblock][j][i] / 
                static_cast<float>(ncc) * 1.0e-4;
#endif // BCKMEX
#endif // TIDEMIX
          }
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < NY_BLOCK-1; ++j) {
        for (int i = 1; i < NX_BLOCK; ++i) {
          // akmu[iblock][k][j][i] = 
          //     au0 [iblock][j][i] * akmt[iblock][k][j  ][i  ] +
          //     aus [iblock][j][i] * akmt[iblock][k][j+1][i  ] +
          //     auw [iblock][j][i] * akmt[iblock][k][j  ][i-1] +
          //     ausw[iblock][j][i] * akmt[iblock][k][j+1][i-1];
          akmu[iblock][k][j][i] = 
              CppConstantMod::P25 * akmt[iblock][k][j  ][i  ] +
              CppConstantMod::P25 * akmt[iblock][k][j+1][i  ] +
              CppConstantMod::P25 * akmt[iblock][k][j  ][i-1] +
              CppConstantMod::P25 * akmt[iblock][k][j+1][i-1];
        }
      }
      for (int i = 0; i < NX_BLOCK; ++i) {
        akmu[iblock][k][NY_BLOCK-1][i] = C0;
      }
      for (int j = 0; j < NY_BLOCK; ++j) {
        akmu[iblock][k][j][0] = C0;
      }
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          akmu[iblock][k][j][i] *= viv[iblock][k+1][j][i];
        }
      }
    }
  }

#endif // CANUTO

  int errorcode;
  pop_haloupdate_readyc_(&errorcode);

  // wjl 20231123
  // for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
  //   for (int j = 0; j < JMT; ++j) {
  //     for (int i = 0; i < IMT; ++i) {
  //       rict_ref[iblock][j][i] = rict[iblock][13][j][i];
  //       for (int k = 0; k < KMM1; ++k) {
  //         if (fabs(zkp[k+1]) > amld[iblock][j][i]) {
  //           rict_ref[iblock][j][i] = rict[iblock][k][j][i];
  //           break;
  //         }
  //       }
  //     }
  //   }
  // }
  upwell(u, v, h0);

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          dlu[iblock][k][j][i] = C0;
          dlv[iblock][k][j][i] = C0;
          wka[iblock][k][j][i] = C0;
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < NY_BLOCK-1; ++j) {
        for (int i = 1; i < NX_BLOCK; ++i) {
          // wka[iblock][k][j][i] = 
          //     au0 [iblock][j][i] * ws[iblock][k][j  ][i  ] +
          //     aus [iblock][j][i] * ws[iblock][k][j+1][i  ] +
          //     auw [iblock][j][i] * ws[iblock][k][j  ][i-1] +
          //     ausw[iblock][j][i] * ws[iblock][k][j+1][i-1];
          wka[iblock][k][j][i] = 
              CppConstantMod::P25 * ws[iblock][k][j  ][i  ] +
              CppConstantMod::P25 * ws[iblock][k][j+1][i  ] +
              CppConstantMod::P25 * ws[iblock][k][j  ][i-1] +
              CppConstantMod::P25 * ws[iblock][k][j+1][i-1];
        }
      }
      for (int i = 0; i < NX_BLOCK; ++i) {
        wka[iblock][k][NY_BLOCK-1][i] = C0;
      }
      for (int j = 0; j < NY_BLOCK; ++j) {
        wka[iblock][k][j][0] = C0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    advection_momentum(u[iblock], v[iblock], wka[iblock],
        dlu[iblock], dlv[iblock], iblock);
  }
#ifdef SMAG
  // TODO
  call smag2(k);

//------------------
#ifdef SMAG_FZ
#else  // SMAG_FZ
#endif // SMAG_FZ
//-----------------

#else // SMAG

  double hduk[JMT][IMT];
  double hdvk[JMT][IMT];

  // const int iiblock = 0;
  // int block_id = blocks_clinic[iiblock];
  // int local_id = iiblock + 1;
  // const struct block this_block = CppBlocks::get_block(&block_id, &local_id);

#ifdef BIHAR
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      // hdiffu_del4(k, hduk, hdvk, up[iblock][k], vp[iblock][k],
      //     this_block);
      hdiffu_del4(k, hduk, hdvk, up[iblock][k], vp[iblock][k]);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          dlv[iblock][k][j][i] += hdvk[j][i];
          dlu[iblock][k][j][i] += hduk[j][i];
        }
      }
    }
  }
#else  // BIHAR
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      hdiffu_del2(k, hduk, hdvk, up[iblock][k], vp[iblock][k],
          this_block);
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          dlv[iblock][k][j][i] += hdvk[j][i];
          dlu[iblock][k][j][i] += hduk[j][i];
        }
      }
    }
  }
#endif // BIHAR
#endif // SMAG

  vinteg(dlu, dlub);
  vinteg(dlv, dlvb);

  int kmb;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        kmb = kmu[iblock][j][i];
        if (kmb >= 0) {
          sbcx[iblock][j][i] = su[iblock][j][i] * od0;
          sbcy[iblock][j][i] = sv[iblock][j][i] * od0;
   
          bbcx[iblock][j][i] = c0f * sqrt(
              up[iblock][kmb-1][j][i] * up[iblock][kmb-1][j][i] + 
              vp[iblock][kmb-1][j][i] * vp[iblock][kmb-1][j][i]) * (
              up[iblock][kmb-1][j][i] * cag + snlat[iblock][j][i] *
              vp[iblock][kmb-1][j][i] * sag);

          bbcy[iblock][j][i] = c0f * sqrt(
              up[iblock][kmb-1][j][i] * up[iblock][kmb-1][j][i] + 
              vp[iblock][kmb-1][j][i] * vp[iblock][kmb-1][j][i]) * (
            - snlat[iblock][j][i] * up[iblock][kmb-1][j][i] * 
              sag + vp[iblock][kmb-1][j][i] * cag);
        } else {
          sbcx[iblock][j][i] = 0.0;
          sbcy[iblock][j][i] = 0.0;
          bbcx[iblock][j][i] = 0.0;
          bbcy[iblock][j][i] = 0.0;
        }
        dlub[iblock][j][i] += (sbcx[iblock][j][i] - bbcx[iblock][j][i]) 
            * ohbu[iblock][j][i];
        dlvb[iblock][j][i] += (sbcy[iblock][j][i] - bbcy[iblock][j][i]) 
            * ohbu[iblock][j][i];
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wka[iblock][k][j][i] = c0f * 
              sqrt(up[iblock][k][j][i] * up[iblock][k][j][i] + 
                   vp[iblock][k][j][i] * vp[iblock][k][j][i]);
        }
      }
    }
  }
#ifdef CANUTO
  const double aidif = 0.5;
#endif // CANUTO
  double diff_u1, diff_u2;
  double diff_v1, diff_v2;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
#ifdef CANUTO
          if (k == 0) {
            diff_v1 = sv[iblock][j][i] * od0 * (1 - aidif);
            diff_u1 = su[iblock][j][i] * od0 * (1 - aidif);
          } else {
            diff_v1 = akmu[iblock][k-1][j][i] * (1 - aidif) *
                (vp[iblock][k-1][j][i] - vp[iblock][k][j][i]) *
                odzt[k] * viv[iblock][k][j][i] + 
                (1.0 - viv[iblock][k][j][i]) * 
                wka[iblock][k-1][j][i] * (1 - aidif) * 
                (-snlat[iblock][j][i] * up[iblock][k-1][j][i] * sag + 
                                        vp[iblock][k-1][j][i] * cag);

            diff_u1 = akmu[iblock][k-1][j][i] * (1 - aidif) *
                (up[iblock][k-1][j][i] - up[iblock][k][j][i]) *
                odzt[k] * viv[iblock][k][j][i] + 
                (1.0 - viv[iblock][k][j][i]) * 
                wka[iblock][k-1][j][i] * (1 - aidif) * 
                (up[iblock][k-1][j][i] * cag + snlat[iblock][j][i] *
                 vp[iblock][k-1][j][i] * sag);
          }
          if (k == KM-1) {
            diff_v2 = wka[iblock][k][j][i] * (-snlat[iblock][j][i] *
                up[iblock][k][j][i] * sag + 
                vp[iblock][k][j][i] * cag) * (1 - aidif);
            diff_u2 = wka[iblock][k][j][i] * (up[iblock][k][j][i] * cag +
                snlat[iblock][j][i] * vp[iblock][k][j][i] * sag) * 
                (1 - aidif);
          } else {
            diff_v2 = akmu[iblock][k][j][i] * (1 - aidif) * 
                (vp[iblock][k][j][i] - vp[iblock][k+1][j][i]) *
                odzt[k+1] * viv[iblock][k+1][j][i] +
                (1.0 - viv[iblock][k+1][j][i]) * wka[iblock][k][j][i] *
                (1 - aidif) * (-snlat[iblock][j][i] * 
                up[iblock][k][j][i] * sag + vp[iblock][k][j][i] * cag);
            diff_u2 = akmu[iblock][k][j][i] * (1 - aidif) *
                (up[iblock][k][j][i] - up[iblock][k+1][j][i]) *
                odzt[k+1] * viv[iblock][k+1][j][i] +
                (1.0 - viv[iblock][k+1][j][i]) * wka[iblock][k][j][i] *
                (1 - aidif) * (up[iblock][k][j][i] * cag +
                snlat[iblock][j][i] * vp[iblock][k][j][i] * sag);
          }
#else  // CANUTO
          if (mytid == 0) {
            printf("The false mixing option\n");
          }
          exit(0);
#endif // CANUTO
          dlv[iblock][k][j][i] += odzp[k] * (diff_v1 - diff_v2);
          dlu[iblock][k][j][i] += odzp[k] * (diff_u1 - diff_u2);
        }
      }
    }
  }
  return ;
}
//--------------------
//  END READYC
//--------------------
static void advection_momentum(
    const double (&uuu)[KM][JMT][IMT],
    const double (&vvv)[KM][JMT][IMT],
    const double (&www)[KM][JMT][IMT],
    double (&adv_uu)[KM][JMT][IMT],
    double (&adv_vv)[KM][JMT][IMT],
    const int &iblock){

  using CppConstantMod::C0;
  using CppConstantMod::P5;
  using CppConstantMod::P25;
  using CppPconstMod::odzp;
  using CppPconstMod::adv_momentum;
  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::hue;
  using CppGrid::hun;
  using CppGrid::uarea_r;

  for (int k = 0; k < KM; ++k) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        adv_uu[k][j][i] = C0;
        adv_vv[k][j][i] = C0;
      }
    }
  }
  std::string str_adv_momentum(adv_momentum);
  double u_wface[KM][JMT][IMT];
  double v_sface[KM][JMT][IMT];
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j][i-1] + uuu[k][j  ][i]) *
              P25 * hue[iblock][j  ][i-1];
          v_sface[k][j][i] = (vvv[k][j][i  ] + vvv[k][j+1][i]) *
              P25 * hun[iblock][j+1][i  ];
        }
      }
    }
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          u_wface[k][j][i] = 
              (uuu[k][j  ][i-1] * dyu[iblock][j  ][i-1] +
               uuu[k][j  ][i  ] * dyu[iblock][j  ][i  ]) * P25;
          v_sface[k][j][i] = 
              (vvv[k][j  ][i  ] * dxu[iblock][j  ][i  ] +
               vvv[k][j+1][i  ] * dxu[iblock][j+1][i  ]) * P25;
        }
      }
    }
  } else {
    if (mytid == 0) {
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  double adv_z1, adv_z2, adv_z3, adv_z4;
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_uu[k][j][i] = (
              - u_wface[k][j  ][i  ] * (uuu[k][j  ][i  ] - uuu[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (uuu[k][j  ][i+1] - uuu[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (uuu[k][j+1][i  ] - uuu[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (uuu[k][j  ][i  ] - uuu[k][j-1][i  ])) 
                  * uarea_r[iblock][j][i];

          adv_vv[k][j][i] = (
              - u_wface[k][j  ][i  ] * (vvv[k][j  ][i  ] - vvv[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (vvv[k][j  ][i+1] - vvv[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (vvv[k][j+1][i  ] - vvv[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (vvv[k][j  ][i  ] - vvv[k][j-1][i  ])) 
                  * uarea_r[iblock][j][i];

          if (k == 0) {
            adv_z1 = 0.0;
            adv_z3 = 0.0;
          } else {
            adv_z1 = www[k][j][i] * (uuu[k-1][j][i] - uuu[k][j][i]);
            adv_z3 = www[k][j][i] * (vvv[k-1][j][i] - vvv[k][j][i]);
          }

          if (k == KM-1) {
            adv_z2 = 0.0;
            adv_z4 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (uuu[k][j][i] - uuu[k+1][j][i]);
            adv_z4 = www[k+1][j][i] * (vvv[k][j][i] - vvv[k+1][j][i]);
          }

          adv_uu[k][j][i] -= P5 * odzp[k] *
              (adv_z1 + adv_z2);
          adv_vv[k][j][i] -= P5 * odzp[k] *
              (adv_z3 + adv_z4);
        }
      }
    }
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_uu[k][j][i] = (
              - u_wface[k][j  ][i  ] * (uuu[k][j  ][i  ] + uuu[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (uuu[k][j  ][i+1] + uuu[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (uuu[k][j  ][i  ] + uuu[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (uuu[k][j+1][i  ] + uuu[k][j  ][i  ])) 
                  * uarea_r[iblock][j][i];

          adv_vv[k][j][i] = (
              - u_wface[k][j  ][i  ] * (vvv[k][j  ][i  ] + vvv[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (vvv[k][j  ][i+1] + vvv[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (vvv[k][j  ][i  ] + vvv[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (vvv[k][j+1][i  ] + vvv[k][j  ][i  ])) 
                  * uarea_r[iblock][j][i];
   
          if (k == 0) {
            adv_z1 = 0.0;
            adv_z3 = 0.0;
          } else {
            adv_z1 = www[k][j][i] * (
                uuu[k-1][j][i] + uuu[k][j][i]) * P5;
            adv_z3 = www[k][j][i] * (
                vvv[k-1][j][i] + vvv[k][j][i]) * P5;
          }
 
          if (k == KM-1) {
            adv_z2 = 0.0;
            adv_z4 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (
                uuu[k][j][i] + uuu[k+1][j][i]) * P5; 
            adv_z4 = www[k+1][j][i] * (
                vvv[k][j][i] + vvv[k+1][j][i]) * P5; 
          }
          adv_uu[k][j][i] -= odzp[k] *
              (adv_z2 - adv_z1);
          adv_vv[k][j][i] -= odzp[k] *
              (adv_z4 - adv_z3);
        }
      }
    }
  } else {
    if (mytid == 0) {
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  return ;
}

#ifndef BIHAR
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
  const int ib = this_block.ib;
  const int ie = this_block.ie;
  const int jb = this_block.jb;
  const int je = this_block.je;
  const int bid = this_block.local_id;
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
#else // BIHAR
static void hdiffu_del4(
    const int &k,
    double (&hduk)[NY_BLOCK][NX_BLOCK],
    double (&hdvk)[NY_BLOCK][NX_BLOCK],
    const double (&umixk)[NY_BLOCK][NX_BLOCK],
    const double (&vmixk)[NY_BLOCK][NX_BLOCK]) {

  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;
  
  using CppConstantMod::C0;
  using CppConstantMod::P5;
  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::dxur;
  using CppGrid::dyur;
  using CppGrid::htw;
  using CppGrid::hts;
  using CppGrid::kmt;
  using CppGrid::kmu;
  using CppGrid::uarea;
  using CppGrid::tarea_r;

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
  //const int bid = this_block.local_id - 1;
  const int bid = 0;
  // TODO this_block

  double div_out[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      div_out[j][i] = C0;
    }
  }
  for (int j = 1; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK-1; ++i) {
      if (k <= kmt[bid][j][i]-1) {
        div_out[j][i] = P5 * (
            (umixk[j  ][i+1] + umixk[j-1][i+1] * htw[bid][j  ][i+1]) -
            (umixk[j  ][i  ] + umixk[j-1][i  ] * htw[bid][j  ][i  ]) +
            (vmixk[j  ][i+1] + vmixk[j  ][i  ] * hts[bid][j  ][i  ]) -
            (vmixk[j-1][i+1] + vmixk[j-1][i  ] * hts[bid][j-1][i  ])) *
                tarea_r[bid][j][i];
      }
    }
  }
  double curl[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      curl[j][i] = C0;
    }
  }
  for (int j = 1; j < NY_BLOCK; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      if (k <= kmt[bid][j][i]-1) {
        curl[j][i] = P5 * (
            vmixk[j  ][i  ] * dyu[bid][j  ][i  ] +
            vmixk[j-1][i  ] * dyu[bid][j-1][i  ] -
            vmixk[j  ][i-1] * dyu[bid][j  ][i-1] -
            vmixk[j-1][i-1] * dyu[bid][j-1][i-1] -
            umixk[j  ][i  ] * dxu[bid][j  ][i  ] -
            umixk[j  ][i-1] * dxu[bid][j  ][i-1] +
            umixk[j-1][i  ] * dxu[bid][j-1][i  ] +
            umixk[j-1][i-1] * dxu[bid][j-1][i-1]);
      }
    }
  }
  double gradx1[NY_BLOCK][NX_BLOCK];
  double grady1[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      gradx1[j][i] = C0;
      grady1[j][i] = C0;
    }
  }
  for (int j = 0; j < NY_BLOCK-1; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      gradx1[j][i] = dxur[bid][j][i] * P5 * (
          curl[j+1][i  ] - curl[j  ][i-1] -
          curl[j+1][i-1] + curl[j  ][i  ]);
      grady1[j][i] = dyur[bid][j][i] * P5 * (
          curl[j+1][i  ] - curl[j  ][i-1] +
          curl[j+1][i-1] - curl[j  ][i  ]);
    }
  }
  double gradx2[NY_BLOCK][NX_BLOCK];
  double grady2[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      gradx2[j][i] = C0;
      grady2[j][i] = C0;
    }
  }
  for (int j = 0; j < NY_BLOCK-1; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      gradx2[j][i] = dxur[bid][j][i] * P5 * (
          div_out[j+1][i  ] - div_out[j  ][i-1] -
          div_out[j+1][i-1] + div_out[j  ][i  ]);
      grady2[j][i] = dyur[bid][j][i] * P5 * (
          div_out[j+1][i  ] - div_out[j  ][i-1] +
          div_out[j+1][i-1] - div_out[j  ][i  ]);
    }
  }

  std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
      am_factor(nblocks_clinic);

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        am_factor[iblock][j][i] = 1.0;
      }
    }
  }
  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      double dxdy = pow(sqrt(uarea[bid][j][i]), 5) * 45.0;
      am_factor[bid][j][i] = sqrt(
          pow(gradx1[j][i], 2) + pow(gradx2[j][i], 2) + 
          pow(grady1[j][i], 2) + pow(grady2[j][i], 2)) *
          dxdy / fabs(am * amf[bid][j][i]);
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
// static void hdiffu_del4(
//     const int &k,
//     double (&hduk)[NY_BLOCK][NX_BLOCK],
//     double (&hdvk)[NY_BLOCK][NX_BLOCK],
//     const double (&umixk)[NY_BLOCK][NX_BLOCK],
//     const double (&vmixk)[NY_BLOCK][NX_BLOCK],
//     const block &this_block) {
  
//   using CppConstantMod::C0;
//   using CppConstantMod::P5;
//   using CppGrid::dxu;
//   using CppGrid::dyu;
//   using CppGrid::dxur;
//   using CppGrid::dyur;
//   using CppGrid::htw;
//   using CppGrid::hts;
//   using CppGrid::kmt;
//   using CppGrid::kmu;
//   using CppGrid::uarea;
//   using CppGrid::tarea_r;

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
//   //const int bid = this_block.local_id - 1;
//   const int bid = 0;
//   // TODO this_block

//   double div_out[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       div_out[j][i] = C0;
//     }
//   }
//   for (int j = 1; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK-1; ++i) {
//       if (k <= kmt[bid][j][i]-1) {
//         div_out[j][i] = P5 * (
//             (umixk[j  ][i+1] + umixk[j-1][i+1] * htw[bid][j  ][i+1]) -
//             (umixk[j  ][i  ] + umixk[j-1][i  ] * htw[bid][j  ][i  ]) +
//             (vmixk[j  ][i+1] + vmixk[j  ][i  ] * hts[bid][j  ][i  ]) -
//             (vmixk[j-1][i+1] + vmixk[j-1][i  ] * hts[bid][j-1][i  ])) *
//                 tarea_r[bid][j][i];
//       }
//     }
//   }
//   double curl[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       curl[j][i] = C0;
//     }
//   }
//   for (int j = 1; j < NY_BLOCK; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       if (k <= kmt[bid][j][i]-1) {
//         curl[j][i] = P5 * (
//             vmixk[j  ][i  ] * dyu[bid][j  ][i  ] +
//             vmixk[j-1][i  ] * dyu[bid][j-1][i  ] -
//             vmixk[j  ][i-1] * dyu[bid][j  ][i-1] -
//             vmixk[j-1][i-1] * dyu[bid][j-1][i-1] -
//             umixk[j  ][i  ] * dxu[bid][j  ][i  ] -
//             umixk[j  ][i-1] * dxu[bid][j  ][i-1] +
//             umixk[j-1][i  ] * dxu[bid][j-1][i  ] +
//             umixk[j-1][i-1] * dxu[bid][j-1][i-1]);
//       }
//     }
//   }
//   double gradx1[NY_BLOCK][NX_BLOCK];
//   double grady1[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       gradx1[j][i] = C0;
//       grady1[j][i] = C0;
//     }
//   }
//   for (int j = 0; j < NY_BLOCK-1; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       gradx1[j][i] = dxur[bid][j][i] * P5 * (
//           curl[j+1][i  ] - curl[j  ][i-1] -
//           curl[j+1][i-1] + curl[j  ][i  ]);
//       grady1[j][i] = dyur[bid][j][i] * P5 * (
//           curl[j+1][i  ] - curl[j  ][i-1] +
//           curl[j+1][i-1] - curl[j  ][i  ]);
//     }
//   }
//   double gradx2[NY_BLOCK][NX_BLOCK];
//   double grady2[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       gradx2[j][i] = C0;
//       grady2[j][i] = C0;
//     }
//   }
//   for (int j = 0; j < NY_BLOCK-1; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       gradx2[j][i] = dxur[bid][j][i] * P5 * (
//           div_out[j+1][i  ] - div_out[j  ][i-1] -
//           div_out[j+1][i-1] + div_out[j  ][i  ]);
//       grady2[j][i] = dyur[bid][j][i] * P5 * (
//           div_out[j+1][i  ] - div_out[j  ][i-1] +
//           div_out[j+1][i-1] - div_out[j  ][i  ]);
//     }
//   }
//   const int ib = this_block.ib;
//   const int ie = this_block.ie;
//   const int jb = this_block.jb;
//   const int je = this_block.je;

//   std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
//       am_factor(nblocks_clinic);

//   for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
//     for (int j = 0; j < NY_BLOCK; ++j) {
//       for (int i = 0; i < NX_BLOCK; ++i) {
//         am_factor[iblock][j][i] = 1.0;
//       }
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       double dxdy = pow(sqrt(uarea[bid][j][i]), 5) * 45.0;
//       am_factor[bid][j][i] = sqrt(
//           pow(gradx1[j][i], 2) + pow(gradx2[j][i], 2) + 
//           pow(grady1[j][i], 2) + pow(grady2[j][i], 2)) *
//           dxdy / fabs(am * amf[bid][j][i]);
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
#endif // BIHAR

#endif // LICOM_ENABLE_FORTRAN
