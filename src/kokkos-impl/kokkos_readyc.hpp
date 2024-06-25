#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYC_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYC_HPP_

#include "../head/def-undef.h"

#include "../head/cpp_blocks.h"
#ifdef CANUTO
#include "../head/cpp_canuto_mod.h"
#endif // CANUTO
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"

#ifdef BIHAR
#include "../head/cpp_hmix_del4.h"
#else // BIHAR
#include "../head/cpp_hmix_del2.h"
#endif // BIHAR

#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"

#ifdef CANUTO
#include "../head/kokkos_canuto_mod.h"
#endif // CANUTO
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_grid.h"

#ifdef BIHAR
#include "../head/kokkos_hmix_del4.h"
#else // BIHAR
#include "../head/kokkos_hmix_del2.h"
#endif // BIHAR

#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_pmix_mod.h"
#include "../head/kokkos_tracer_mod.h"
#include "../head/kokkos_tmp_var.h"
#include "../head/kokkos_work_mod.h"

#include "../head/kokkos_config.hpp"

#include "../head/fortran_blocks.h"
#ifdef CANUTO
#include "../head/fortran_canuto_mod.h"
#endif // CANUTO
#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include "Kokkos_Core.hpp"

#include <cmath>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <string>

using CppParamMod::mytid;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::KMM1;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

using CppConstantMod::C0;
using CppConstantMod::P5;
using CppConstantMod::P25;

using CppPconstMod  ::ncc;
using CppPconstMod  ::MIXING_EF;
using CppPconstMod  ::MAX_TIDALMIXING;
using CppPconstMod  ::BACK_TIDALMIXING;
using CppPconstMod  ::LOCAL_MIXING_FRACTION;

#ifdef CANUTO
using KokkosCanutoMod:: p_v_rib;
using KokkosCanutoMod:: p_v_ridb;
using KokkosCanutoMod:: p_v_irimax;
using KokkosCanutoMod:: p_v_sma1;
using KokkosCanutoMod:: p_v_sha1;
using KokkosCanutoMod:: p_v_ssa1;
using KokkosCanutoMod:: p_v_smb;
using KokkosCanutoMod:: p_v_shb;
using KokkosCanutoMod:: p_v_ssb;
using KokkosCanutoMod:: p_v_slq2b;
using KokkosCanutoMod:: p_v_sm_r1;
using KokkosCanutoMod:: p_v_sh_r1;
using KokkosCanutoMod:: p_v_ss_r1;
using KokkosCanutoMod:: p_v_slq2_r1;
using KokkosCanutoMod:: p_v_and2on2a1;
using KokkosCanutoMod:: p_v_amtaun2a1;
using KokkosCanutoMod:: p_v_back_ra_r;
#endif // CANUTO
using KokkosDynMod::    p_v_dlu;
using KokkosDynMod::    p_v_dlv;
using KokkosDynMod::    p_v_dlub;
using KokkosDynMod::    p_v_dlvb;
using KokkosDynMod::    p_v_h0;
using KokkosDynMod::    p_v_h0bl;
using KokkosDynMod::    p_v_h0bf;
using KokkosDynMod::    p_v_u;
using KokkosDynMod::    p_v_v;
using KokkosDynMod::    p_v_up;
using KokkosDynMod::    p_v_vp;
using KokkosDynMod::    p_v_ws;
using KokkosForcMod::   p_v_buoysol;
using KokkosForcMod::   p_v_buoytur;
using KokkosForcMod::   p_v_su;
using KokkosForcMod::   p_v_sv;
using KokkosForcMod::   p_v_ustar;
using KokkosForcMod::   p_v_wave_dis;
using KokkosGrid::      p_v_at0;
using KokkosGrid::      p_v_atn;
using KokkosGrid::      p_v_ate;
using KokkosGrid::      p_v_atne;
using KokkosGrid::      p_v_dxu;
using KokkosGrid::      p_v_dyu;
using KokkosGrid::      p_v_dxyur;
using KokkosGrid::      p_v_hue;
using KokkosGrid::      p_v_hun;
using KokkosGrid::      p_v_htw;
using KokkosGrid::      p_v_hts;
using KokkosGrid::      p_v_kmt;
using KokkosGrid::      p_v_kmu;
using KokkosGrid::      p_v_ulat;
using KokkosGrid::      p_v_fcort;
using KokkosGrid::      p_v_uarea;
using KokkosGrid::      p_v_tarea_r;
using KokkosGrid::      p_v_uarea_r;
using KokkosPconstMod:: p_v_akmu;
using KokkosPconstMod:: p_v_akt;
using KokkosPconstMod:: p_v_akmt;
using KokkosPconstMod:: p_v_dzp;
using KokkosPconstMod:: p_v_fztidal;
using KokkosPconstMod:: p_v_fz_tide;
using KokkosPconstMod:: p_v_richardson;
using KokkosPconstMod:: p_v_viv;
using KokkosPconstMod:: p_v_vit;
using KokkosPconstMod:: p_v_odzt;
using KokkosPconstMod:: p_v_odzp;
using KokkosPconstMod:: p_v_odz_pt;
using KokkosPconstMod:: p_v_ohbu;
using KokkosPconstMod:: p_v_ohbt;
using KokkosPconstMod:: p_v_to;
using KokkosPconstMod:: p_v_so;
using KokkosPconstMod:: p_v_snlat;
using KokkosPconstMod:: p_v_wp3_tidal;
using KokkosPconstMod:: p_v_zkp;
using KokkosPmixMod::   p_v_rit;
using KokkosPmixMod::   p_v_rict;
using KokkosPmixMod::   p_v_ricdt;
using KokkosPmixMod::   p_v_ricdttms;
using KokkosPmixMod::   p_v_s2t;
using KokkosTracerMod:: p_v_amld;
using KokkosTracerMod:: p_v_at;
using KokkosTracerMod:: p_v_pdensity;
using KokkosWorkMod::   p_v_uk;
using KokkosWorkMod::   p_v_vk;
using KokkosWorkMod::   p_v_wka;
using KokkosWorkMod::   p_v_wkb;
using KokkosWorkMod::   p_v_work;
#ifdef BIHAR
using KokkosHmixDel4::  p_v_amf;

using KokkosHmixDel4::  p_v_du_cnsewm;
using KokkosHmixDel4::  p_v_dm_cnsew;
using KokkosTmpVar::p_v_curl;
using KokkosTmpVar::p_v_d2uk;
using KokkosTmpVar::p_v_d2vk;
#else // BIHAR
using KokkosHmixDel2::  p_v_duc;
using KokkosHmixDel2::  p_v_dum;
using KokkosHmixDel2::  p_v_dun;
using KokkosHmixDel2::  p_v_dus;
using KokkosHmixDel2::  p_v_due;
using KokkosHmixDel2::  p_v_duw;
using KokkosHmixDel2::  p_v_dmc;
using KokkosHmixDel2::  p_v_dmn;
using KokkosHmixDel2::  p_v_dms;
using KokkosHmixDel2::  p_v_dme;
using KokkosHmixDel2::  p_v_dmw;
#endif // BIHAR

using KokkosTmpVar::p_v_wp12;
using KokkosTmpVar::p_v_wp13;
using KokkosTmpVar::p_v_uv_ws_face;

using KokkosTmpVar::p_v_ak_tide_mixing;
// using KokkosTmpVar::p_v_mld_lev;
using KokkosTmpVar::p_v_zlev;
using KokkosTmpVar::p_v_wp3;

class FunctorReadyc1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_amld_(iblock, j, i) = C0;
    v_h0bl_(iblock, j, i) = v_h0bf_(iblock, j, i);
    v_h0bf_(iblock, j, i) = v_h0_(iblock, j, i);
    return ;
  }
 private:
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_amld_ = *p_v_amld;
  const ViewDouble3D v_h0bf_ = *p_v_h0bf;
  const ViewDouble3D v_h0bl_ = *p_v_h0bl;
};

//---------------------------------------------------------------------
// Calculating Richardson number riu (at U/V-point);
//-------------------------------------------------------------------
// using wka for ...
// FunctorReadyc2: ugrid_to_tgrid
// vit[iblock][k][j][i] / 
//     (viv[iblock][k][j  ][i  ] * at0 [iblock][j][i] +
//      viv[iblock][k][j-1][i  ] * atn [iblock][j][i] +
//      viv[iblock][k][j  ][i+1] * ate [iblock][j][i] +
//      viv[iblock][k][j-1][i+1] * atne[iblock][j][i] + 
//         epsln);
// -------------------------------------------------------------------
// class FunctorReadyc2 {
//  public:
//   KOKKOS_INLINE_FUNCTION void operator () (
//       const int &k, const int &j, const int &i) const {
//     const int iblock = 0;
//     const double epsln = 1.0e-25;
//     if (i < (NX_BLOCK-1) && j >= 1) {
//       v_wka_(iblock, k, j, i) = v_vit_(iblock, k, j, i) / (
//           v_viv_(iblock, k, j  , i  ) * v_at0_ (iblock, j, i) +
//           v_viv_(iblock, k, j-1, i  ) * v_atn_ (iblock, j, i) +
//           v_viv_(iblock, k, j  , i+1) * v_ate_ (iblock, j, i) +
//           v_viv_(iblock, k, j-1, i+1) * v_atne_(iblock, j, i)
//               + epsln);
//     }
//     return ;
//   }
//  private:
//   const ViewDouble3D v_at0_  = *p_v_at0;
//   const ViewDouble3D v_atn_  = *p_v_atn;
//   const ViewDouble3D v_ate_  = *p_v_ate;
//   const ViewDouble3D v_atne_ = *p_v_atne;
//   const ViewDouble4D v_viv_  = *p_v_viv;
//   const ViewDouble4D v_vit_  = *p_v_vit;
//   const ViewDouble4D v_wka_  = *p_v_wka;
// };

class FunctorReadyc3 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_akmt_(iblock, k, j, i) = C0;
    v_akmu_(iblock, k, j, i) = C0;
    ugrid_to_tgrid (iblock, k, j, i, v_wp12_, v_up_);
    ugrid_to_tgrid (iblock, k, j, i, v_wp13_, v_vp_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void ugrid_to_tgrid (const int &iblock, 
			const int &k, const int &j, const int &i,
          const ViewDouble4D &v_tgrid,
              const ViewDouble4D &v_ugrid) const {
    if (i < (NX_BLOCK-1) && j >= 1) {
      const double epsln = 1.0e-25;
      v_tgrid(iblock, k, j, i) =  v_vit_ (iblock, k, j  , i  ) * 
         (v_at0_ (iblock, j, i) * v_ugrid(iblock, k, j  , i  ) + 
          v_atn_ (iblock, j, i) * v_ugrid(iblock, k, j-1, i  ) +
          v_ate_ (iblock, j, i) * v_ugrid(iblock, k, j  , i+1) +
          v_atne_(iblock, j, i) * v_ugrid(iblock, k, j-1, i+1)) / 
         (v_viv_(iblock, k, j  , i  ) * v_at0_ (iblock, j, i) +
          v_viv_(iblock, k, j-1, i  ) * v_atn_ (iblock, j, i) +
          v_viv_(iblock, k, j  , i+1) * v_ate_ (iblock, j, i) +
          v_viv_(iblock, k, j-1, i+1) * v_atne_(iblock, j, i)
              + epsln);
    }
    if (i == (NX_BLOCK-1) || j == 0) {
      v_tgrid(iblock, k, j, i) = C0;
    }
    return ;
  }
 private:
  const ViewDouble3D v_at0_  = *p_v_at0;
  const ViewDouble3D v_atn_  = *p_v_atn;
  const ViewDouble3D v_ate_  = *p_v_ate;
  const ViewDouble3D v_atne_ = *p_v_atne;
  const ViewDouble4D v_up_   = *p_v_up;
  const ViewDouble4D v_vp_   = *p_v_vp;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_wp12_ = *p_v_wp12;
  const ViewDouble4D v_wp13_ = *p_v_wp13;
  const ViewDouble4D v_akmt_ = *p_v_akmt;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
};

class FunctorReadyc4 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    const double epsln = 1.0e-25;
    if (i == (IMT-1) || j == 0) {
      v_s2t_(iblock, k, j, i) = C0;
    }
    if (i < (IMT-1) && j >= 1) {
      const double riv1 =
          v_wp12_(iblock, k  , j, i) * v_vit_(iblock, k  , j, i) -
          v_wp12_(iblock, k+1, j, i) * v_vit_(iblock, k+1, j, i);
      const double riv2 =
          v_wp13_(iblock, k  , j, i) * v_vit_(iblock, k  , j, i) -
          v_wp13_(iblock, k+1, j, i) * v_vit_(iblock, k+1, j, i);
      v_s2t_(iblock, k, j, i) = v_vit_(iblock, k+1, j, i) *
          (riv1 * riv1 + riv2 * riv2) * v_odzt_(k+1) * v_odzt_(k+1);
#if (defined CANUTO) || (defined CANUTO2010)
      v_rit_(iblock, k, j, i) = v_vit_(iblock, k+1, j, i) 
          * v_rict_(iblock, k, j, i) 
              / (v_s2t_(iblock, k, j, i) + epsln);
#else
      v_rit_(iblock, k, j, i) += v_vit_(iblock, k+1, j, i) 
          * v_rict_(iblock, k, j, i)
              / (v_s2t_(iblock, k, j, i) + epsln);
#endif // CANUTO
    }
    return ;
  }
 private:
  const ViewDouble1D v_odzt_ = *p_v_odzt; 
  const ViewDouble4D v_vit_  = *p_v_vit; 
  const ViewDouble4D v_s2t_  = *p_v_s2t; 
  const ViewDouble4D v_rit_  = *p_v_rit; 
  const ViewDouble4D v_rict_ = *p_v_rict; 
  const ViewDouble4D v_wp12_ = *p_v_wp12;
  const ViewDouble4D v_wp13_ = *p_v_wp13;
} ;

#if (defined CANUTO) || (defined CANUTO2010)
class FunctorReadyc5 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;

      for (int k = 0; k < KM; ++k) {
        v_wk1_(j, i, k) = 0.0;
        v_wk2_(j, i, k) = 0.0;
        v_wk3_(j, i, k) = 0.0;
        v_wk4_(j, i, k) = 0.0;
        v_wp1_(j, i, k) = 0.0;
        v_wp2_(j, i, k) = 0.0;
        v_wp3_(j, i, k) = 0.0;
        v_wp4_(j, i, k) = 0.0;
        v_wp5_(j, i, k) = 0.0;
        v_wp6_(j, i, k) = 0.0;
        v_wp7_(j, i, k) = 0.0;
        v_wp8_(j, i, k) = 0.0;
      }
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;
      for (int k = 0; k < KM - 1; ++k) {
        v_wp8_(j, i, k) = - v_vit_(iblock, k+1, j, i) * v_zkp_(k+1);
      }

      for (int k = 0; k < kmt_m1; ++k) {
        if (v_vit_(iblock, k+1, j, i) > 0.0) {
          v_wp1_(j, i, k) = v_at_(iblock, 0, k, j, i) - (
              v_at_(iblock, 0, k, j, i) - v_at_(iblock, 0, k+1, j, i)) 
                  * v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1));
          v_wp2_(j, i, k) = (v_at_(iblock, 1, k, j, i) - (
              v_at_(iblock, 1, k, j, i) - v_at_(iblock, 1, k+1, j, i)) 
                  * v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1))) * 1000.0 + 35.0;
          v_wp4_(j, i, k) = v_rit_(iblock, k, j, i);
          v_wp5_(j, i, k) = v_ricdt_(iblock, k, j, i);
          v_wp6_(j, i, k) = v_s2t_(iblock, k, j, i);
          v_wp7_(j, i, k) = v_rict_(iblock, k, j, i);
        }

      }
      for (int k = 0; k < v_kmt_(iblock, j, i); ++k) {
        const double tq = v_at_(iblock, 0, k, j, i) - v_to_(0);
        const double sq = v_at_(iblock, 1, k, j, i) - v_so_(0);
        v_wp3_(j, i, k) = dens(tq, sq, 0) + 1000.0;
      }
      for (int k = 0; k < KM; ++k) {
        if (v_wp7_(j, i, k) < 0.0) {
          v_wp7_(j, i, k) = 0.0;
        }
      }

#ifdef CANUTOMIXOUT
      // To do something
#endif // CANUTOMIXOUT
#ifdef CANUTO2010
      canuto_2010_interface (j, i,
          v_wk1_, v_wk2_, v_wk3_, v_wk4_, v_amld_(iblock, j, i),
          v_wp1_, v_wp2_, v_wp3_, v_wp4_, v_wp5_,
          v_wp7_, v_wp6_, v_ulat_(iblock, j, i) / DegToRad_, 
          v_wp8_, v_kmt_(iblock, j, i));
#endif // CANUTO2010

#ifdef TIDEMIX
      //double ak_tide_mixing[KM];
      for (int k = 0; k < KM; ++k) {
        v_ak_tide_mixing_(j, i, k) = 0.0;
      }
      for (int k = 0; k < kmt_m1; ++k) {
        v_ak_tide_mixing_(j, i, k) = BACK_TIDALMIXING + MIXING_EF * 
            LOCAL_MIXING_FRACTION * v_wave_dis_(iblock, j, i) * 
            v_fz_tide_(iblock, k, j, i) /
            (fmax(v_rict_(iblock, k, j, i), 1.0e-8) * v_wp3_(j, i, k));

        v_ak_tide_mixing_(j, i, k) = fmin(v_ak_tide_mixing_(j, i, k),
            MAX_TIDALMIXING);

        // v_richardson_(iblock, k, j, i) = v_rict_(iblock, k, j, i);
        // v_fztidal_(iblock, k, j, i)    = v_fz_tide_(iblock, k, j, i);
        // v_wp3_tidal_(iblock, k, j, i)  = v_wp3_(j, i, k);
      }
#ifdef CANUTOMIXOUT
#endif // CANUTOMIXOUT
      for (int k = kmt_m1 - 2; k >= 0; --k) {
        v_ak_tide_mixing_(j, i, k) = fmin(
            v_ak_tide_mixing_(j, i, k), v_ak_tide_mixing_(j, i, k+1));
      }
#endif // TIDEMIX
      for (int k = 0; k < KM; ++k) {
        v_akmt_(iblock, k, j, i)    = v_wk1_(j, i, k);
        v_akt_(iblock, 0, k, j, i) += v_wk2_(j, i, k);
        v_akt_(iblock, 1, k, j, i) += v_wk3_(j, i, k);
#ifdef BCKMEX
        v_akmt_(iblock, k, j, i) += diff_back[iblock][j][i] * 
            10.0 * 1.0e-4;
        v_akt_(iblock, 0, k, j, i) += diff_back[iblock][j][i] /
                static_cast<float>(ncc) * 1.0e-4;
        v_akt_(iblock, 1, k, j, i) += diff_back[iblock][j][i] /
                static_cast<float>(ncc) * 1.0e-4;
#endif // BCKMEX
#ifdef TIDEMIX
        v_akmt_(iblock, k, j, i)   += (v_ak_tide_mixing_(j, i, k) * 5.0);
        v_akt_(iblock, 0, k, j, i) += v_ak_tide_mixing_(j, i, k);
        v_akt_(iblock, 1, k, j, i) += v_ak_tide_mixing_(j, i, k);
#endif // TIDEMIX
        const double tmp_dzp = 0.008 * v_dzp_(k) * v_dzp_(k);
        v_akmt_(iblock, k, j, i)   = std::min(v_akmt_(iblock, k, j, i),
            tmp_dzp);
        v_akt_(iblock, 0, k, j, i) = std::min(v_akt_(iblock, 0, k, j, i), tmp_dzp);
        v_akt_(iblock, 1, k, j, i) = std::min(v_akt_(iblock, 1, k, j, i), tmp_dzp);
      }
    }
    return ;
  };
 private:
  const double OMEGA_      = CppConstantMod::OMEGA;
  const double DegToRad_   = CppConstantMod::DEGTORAD;
  const double KARMAN_     = CppConstantMod::KARMAN;
  const double VERY_SMALL_ = CppConstantMod::VERY_SMALL;
  const ViewInt3D    v_kmt_            = *p_v_kmt;
  const ViewDouble3D v_wk1_            = *KokkosTmpVar::p_v_wk1;
  const ViewDouble3D v_wk2_            = *KokkosTmpVar::p_v_wk2;
  const ViewDouble3D v_wk3_            = *KokkosTmpVar::p_v_wk3;
  const ViewDouble3D v_wk4_            = *KokkosTmpVar::p_v_wk4;
  const ViewDouble3D v_wp1_            = *KokkosTmpVar::p_v_wp1;
  const ViewDouble3D v_wp2_            = *KokkosTmpVar::p_v_wp2;
  const ViewDouble3D v_wp3_            = *KokkosTmpVar::p_v_wp3;
  const ViewDouble3D v_wp4_            = *KokkosTmpVar::p_v_wp4;
  const ViewDouble3D v_wp5_            = *KokkosTmpVar::p_v_wp5;
  const ViewDouble3D v_wp6_            = *KokkosTmpVar::p_v_wp6;
  const ViewDouble3D v_wp7_            = *KokkosTmpVar::p_v_wp7;
  const ViewDouble3D v_wp8_            = *KokkosTmpVar::p_v_wp8;
  const ViewDouble3D v_zlev_           = *KokkosTmpVar::p_v_zlev;
  const ViewDouble3D v_ak_tide_mixing_ = *KokkosTmpVar::p_v_ak_tide_mixing;
  const ViewDouble3D v_Ri_             = *KokkosTmpVar::p_v_Ri;
  const ViewDouble3D v_Rrho_           = *KokkosTmpVar::p_v_Rrho;
  const ViewDouble3D v_Gm_             = *KokkosTmpVar::p_v_Gm;
  const ViewDouble1D v_to_             = *p_v_to;
  const ViewDouble1D v_so_             = *p_v_so;
  const ViewDouble1D v_dzp_            = *p_v_dzp;
  const ViewDouble1D v_zkp_            = *p_v_zkp;
  const ViewDouble2D v_c_              = *KokkosPconstMod::p_v_c;
  const ViewDouble3D v_ulat_           = *p_v_ulat;
  const ViewDouble3D v_amld_           = *p_v_amld;
  const ViewDouble3D v_wave_dis_       = *p_v_wave_dis;
  const ViewDouble4D v_rit_            = *p_v_rit;
  const ViewDouble4D v_rict_           = *p_v_rict;
  const ViewDouble4D v_ricdt_          = *p_v_ricdt;
  const ViewDouble4D v_akmt_           = *p_v_akmt;
  // const ViewDouble4D v_fztidal_        = *p_v_fztidal;
  // const ViewDouble4D v_wp3_tidal_      = *p_v_wp3_tidal;
  // const ViewDouble4D v_richardson_     = *p_v_richardson;
  const ViewDouble4D v_fz_tide_        = *p_v_fz_tide;
  const ViewDouble4D v_s2t_            = *p_v_s2t;
  const ViewDouble4D v_vit_            = *p_v_vit;
  const ViewDouble5D v_at_             = *p_v_at;
  const ViewDouble5D v_akt_            = *p_v_akt;

  KOKKOS_INLINE_FUNCTION double dens (
      const double &tq, const double &sq, const int &kk) const {
    double dens;
    dens = (v_c_(0, kk) + (v_c_(3, kk) + v_c_(6, kk) * sq) * sq +
           (v_c_(2, kk) +  v_c_(7, kk) * sq + v_c_(5, kk) * tq) * tq) * tq +
           (v_c_(1, kk) + (v_c_(4, kk) + v_c_(8, kk) * sq) * sq) * sq;
    return dens;
  }

  /*
  template<typename T> KOKKOS_INLINE_FUNCTION 
  T sign (const T &x, const T &y) const {
    return y >= static_cast<T>(0) ? std::abs(x) : -std::abs(x);
  }
  */
  KOKKOS_INLINE_FUNCTION 
  double sign (const double &x, const double &y) const {
    return y >= 0.0 ? fabs(x) : - fabs(x);
  }

  KOKKOS_INLINE_FUNCTION 
  void Rf_calc (double &Rf_out, const double &Ri_in, const double &Rrho_in, 
      const double &struct_h, const double &struct_m, const double &r_in) const {
    Rf_out = Ri_in * struct_h / struct_m * (1.0 - Rrho_in / r_in);
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void mixed_layer_TKE_calc (double &TKE_out, const double &mld_in, 
      const double &Gm_in, const double &s2_in, const double &lev_in, 
          const double &Rf_in, const double &Rf_inf_in) const {
    const double l0 = 0.17 * mld_in;
    const double lB = KARMAN_ * std::abs(lev_in * l0) / (l0 + KARMAN_ * std::abs(lev_in));
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    const double tmp = std::pow(std::pow(1.0 - Rf_in / Rf_inf_in, 4.0), 1.0 / 3.0); 
    const double l = lB * tmp;
    const double B1 = 21.6;
    TKE_out = B1 * B1 * std::pow(Gm_in, - 3.0 / 2.0)
        * l * l * std::pow(s2_in, 3.0/2.0);
#else
    const double tmp = static_cast<double>(
        std::pow(
        std::pow(static_cast<long double>(1.0 - Rf_in / Rf_inf_in), 4.0), 
            1.0 / 3.0)); 
    const double l = lB * tmp;
    const double B1 = 21.6;
    TKE_out = B1 * B1 * 
        static_cast<double>(std::pow(static_cast<long double>(Gm_in), - 3.0 / 2.0))
            * l * l * static_cast<double>(
                std::pow(static_cast<long double>(s2_in), 3.0/2.0));
#endif

    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void thermocline_mixing_coeff_calc (double &Km_out, 
             double &Kh_out,          double &Ks_out,          double &Kd_out,
       const double &mix_eff_m, const double &mix_eff_h, const double &mix_eff_s, 
       const double &mix_eff_d, const double &n2_in,     const double &lat_in) 
          const {
    const double N0     = 5.24e-3;
    const double F30    = 2.0 * OMEGA_ * std::sin(30.0 * DegToRad_);
    const double Tke_n2 = 0.288e-4;
 
    const double f_lat = std::abs (2.0 * OMEGA_ * std::sin(lat_in * DegToRad_) + VERY_SMALL_);
 
    const double cond1 = std::max(std::sqrt(n2_in) / f_lat, 1.0);
 
    // double L_lat = (f_lat * std::acosh(cond1)) / (F30 * std::acosh(N0 / F30));
    // const double L_lat = (f_lat * std::acosh(cond1)) / (F30 * std::acosh(N0 / F30));
 
    // Km_out = L_lat * Tke_n2 * mix_eff_m;
    // Kh_out = L_lat * Tke_n2 * mix_eff_h;
    // Ks_out = L_lat * Tke_n2 * mix_eff_s;
    // Kd_out = L_lat * Tke_n2 * mix_eff_d;
    const double L_lat = (f_lat * acosh(cond1)) / (F30 * acosh(N0 / F30)) * Tke_n2;
 
    Km_out = L_lat * mix_eff_m;
    Kh_out = L_lat * mix_eff_h;
    Ks_out = L_lat * mix_eff_s;
    Kd_out = L_lat * mix_eff_d;
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void shengjin_calc (const double* array, const int &nmax, 
      double& real_x) const {
    const double epsln1 = 1.0e-6;
 
    const double sa = array[0];
    const double sb = array[1];
    const double sc = array[2];
    const double sd = array[3];
 
    real_x = 1.0;
 
    if (std::abs(sa) <= epsln1 && std::abs(sb) <= epsln1) {
      if ((-sd / sc) < 0.0) {
        // printf ("No real positive root at position 1\n");
        return;
      } else {
        // printf ("Only one real positive root exists 1!\n");
        real_x = -sd / sc;
        return;
      }
    }
 
    double delta;
    double x[3];
    if (std::abs(sa) <= epsln1 && std::abs(sb) >= epsln1) {
      delta = sc * sc - 4.0 * sb * sd;
      if (delta < 0.0) {
        // printf ("No real roots in this secondary order equation 2!\n");
        return;
      } else if (delta >= 0.0) {
        const double tmp = std::sqrt(delta);
        x[0] = (-sc + tmp) / (2.0 * sb);
        x[1] = (-sc - tmp) / (2.0 * sb);
        if (x[0] < 0.0 && x[1] < 0.0) {
          // std::cout << "No real positive root at position 2" << std::endl;
          return;
        } else {
          // std::cout << "Two real with one positive at least 2!" << std::endl;
          if (x[0] > 0.0 && x[1] > 0.0) {
            real_x = std::min(x[0], x[1]);
          } else if (x[0] >  0.0 && x[1] <= 0.0) {
            real_x = x[0];
          } else if (x[0] <= 0.0 && x[1]  > 0.0) {
            real_x = x[1];
          }
          return;
        }
      }
    }
 
    const double A = sb * sb - 3.0 * sa * sc;
    const double B = sb * sc - 9.0 * sa * sd;
    const double C = sc * sc - 3.0 * sb * sd;
 
    if (std::abs(A) <= epsln1 && std::abs(B) <= epsln1) {
      if ((-sc / sb) < 0.0) {
        // std::cout << "No real positive root at position 3" << std::endl;
        return;
      } else {
        // std::cout << "Only one real positive root exists 3!" << std::endl;
        real_x = -sc / sb;
        return;
      }
    }
 
    delta = B * B - 4.0 * A * C;
 
    double y[2];
    if (delta > 0.0) {
      double tmp = std::sqrt(delta);
      y[0] = A * sb + 3.0 * sa * (-B + tmp) / 2.0;
      y[1] = A * sb + 3.0 * sa * (-B - tmp) / 2.0;
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
      real_x = (-sb - sign(std::pow(std::abs(y[0]), 1.0 / 3.0), y[0]) 
                    - sign(std::pow(std::abs(y[1]), 1.0 / 3.0), y[1])) / (3.0 * sa);
#else
      real_x = (-sb - 
          sign(static_cast<double>(
              std::pow(static_cast<long double>(std::abs(y[0])), 1.0 / 3.0)), y[0]) 
        - sign(static_cast<double>(
              std::pow(static_cast<long double>(std::abs(y[1])), 1.0 / 3.0)), y[1])) 
                  / (3.0 * sa);
#endif
      if (real_x < 0.0) {
        real_x = 1.0;
        return;
      } else {
        // std::cout << "Only one real positive root exists 4!" << std::endl;
        // real_x = real_x;
        return;
      }
    }
 
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    const double T = (2.0 * A * sb - 3.0 * sa * B) / (2.0 * std::pow(A, 1.5));
#else
    const double T = (2.0 * A * sb - 3.0 * sa * B) / (2.0 
        * static_cast<double>(std::pow(static_cast<long double>(A), 1.5)));
#endif
    const double the  = std::acos(T);
    const double tmp1 = std::sqrt(A);
    const double tmp2 = std::cos(the / 3.0);
    const double tmp3 = std::sqrt(3.0) * std::sin(the / 3.0);
    x[0] = (-sb - 2.0  *  tmp1 *  tmp2        ) / 3.0 / sa;
    x[1] = (-sb +         tmp1 * (tmp2 + tmp3)) / 3.0 / sa;
    x[2] = (-sb +         tmp1 * (tmp2 - tmp3)) / 3.0 / sa;
 
    real_x = 500.0;
    for (int k = 0; k < 3; ++k) {
      if (x[k] >= 0.0 && x[k] <= 1000.0) {
        if (x[k] < real_x) {
          real_x = x[k];
        }
      }
    }
 
    real_x = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; ++i) {
      if (x[i] > 0.0) {
        real_x = std::min (real_x, x[i]);
      }
    }
 
    real_x = std::min (real_x, 2.0e3);
    real_x = std::max (real_x, 1.0);
 
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void dyn_time_scale_calc (const double &Ri, const double &Rrho, 
      const double* pi_n, double &Gm, double* const cubic_coef) const {
 
    double A1, A2, A3, A4, A5, A6;
    double temp1, temp2, temp3, temp4, temp5;
    double temp6, temp7, temp8, temp9, temp10;
    double temp11;
 
    temp1 = pi_n[0] * pi_n[3] * (pi_n[3] - pi_n[0] * Rrho);
    temp2 = pi_n[1] * (15.0 * pi_n[2] + 7.0) * (Rrho * Rrho + 1.0);
    temp3 = (14.0 * (pi_n[1] - pi_n[2]) - 15.0 * pi_n[2] * pi_n[2]) * Rrho;
    // hwy_need_check
    const double one_minus_Rrho = 1.0 - Rrho;
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    temp4 = 150.0 * one_minus_Rrho * one_minus_Rrho * one_minus_Rrho; 
#else
    long double tmp = static_cast<long double>(one_minus_Rrho);
    temp4 = 150.0 * static_cast<double>(tmp * tmp * tmp); 
#endif
 
    A1 = temp1 * (temp2 + temp3) / temp4;
 
    temp1 = pi_n[0] * pi_n[3];
    temp2 = pi_n[1] * (210.0 * pi_n[0] - 150.0 * pi_n[2] + 7.0) 
        * (Rrho * Rrho + 1.0);
    temp3 = 14.0 * (pi_n[1] - pi_n[2]) * (1.0 + 15.0 * pi_n[0] + 15.0 * pi_n[3]);
    temp4 = 150.0 * pi_n[2] * pi_n[2];
    temp5 = (temp3 + temp4) * Rrho;
    temp6 = 210.0 * pi_n[1] * (pi_n[3] - pi_n[0]);
 
    // hwy_need_check
    temp7 = 9000.0 * (one_minus_Rrho * one_minus_Rrho); 
 
    A2 = temp1 * (temp2 + temp5 + temp6) / temp7;
 
    temp1 = pi_n[0];
    temp2 = 5.0 * pi_n[1] * pi_n[3] * (30.0 * pi_n[2] + 17.0);
    temp3 = pi_n[0] * (15.0 * pi_n[2] + 7.0);
    temp4 = Rrho * Rrho + 1.0;
    temp5 = temp1 * (temp2 + temp3) * temp4;
    temp6 = -(15.0 * pi_n[2] + 7.0) * (pi_n[0] * pi_n[0] - pi_n[3] * pi_n[3]);
    temp7 = 10.0 * pi_n[0] * pi_n[2] * pi_n[3] * (15.0 * pi_n[2] + 17.0);
    temp8 = 15.0 * pi_n[1]          * (pi_n[0] * pi_n[0] + pi_n[3] * pi_n[3]);
    temp9 = 14.0 * pi_n[0] * pi_n[3] * (1.0 - 10.0 * pi_n[1]);
    temp10 = -(temp7 + temp8 + temp9) * Rrho;
    temp11 = 150.0 * (one_minus_Rrho * one_minus_Rrho);  // hwy_need_check
 
    A3 = (temp5 + temp6 + temp10) / temp11;
 
    temp1 = 150.0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
    temp2 = -7.0 * pi_n[0] * (1.0 + 30.0 * pi_n[0]);
    temp3 = (temp1 + temp2) * Rrho;
    temp4 = -150.0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
    temp5 = 7.0 * pi_n[3] * (1.0 + 30.0 * pi_n[3]);
    temp6 = 9000.0 * one_minus_Rrho;  // hwy_need_check
 
    A4 = (temp3 + temp4 + temp5) / temp6;
 
    temp1 = -30.0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
    temp2 = -17.0 * pi_n[0];
    temp3 = (temp1 + temp2) * Rrho;
    temp4 = 30.0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
    temp5 = 17.0 * pi_n[3];
    temp6 = 30.0 * one_minus_Rrho;  // hwy_need_check
 
    A5 = (temp3 + temp4 + temp5) / temp6;
 
    A6 = - 1.0 / 60.0;
 
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    cubic_coef[0] = A1 * (Ri * Ri * Ri) + A2 * (Ri * Ri);
#else
    tmp = static_cast<long double>(Ri);
    cubic_coef[0] = A1 * static_cast<double>(tmp * tmp * tmp) + A2 * (Ri * Ri);
#endif
    cubic_coef[1] = A3 * (Ri * Ri) + A4 * Ri;
    cubic_coef[2] = A5 * Ri + A6;
    cubic_coef[3] = 1.0;
    shengjin_calc(cubic_coef, 4, Gm);
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void prepare_pi (const double &Ri, const double &Rrho, 
      double* const out_pi) const {
    const double a     = 10.0;
    const double Ko    = 1.66;
    const double sig_t = 0.72;
    const double zero  = 1.0e-10;
 
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    const double pi0_1 = 1.0 / (std::sqrt(27.0 / 5.0 * std::pow(Ko, 3)) 
        * (1.0 + 1.0 / sig_t));
#else
    const double pi0_1 = 1.0 / (std::sqrt(27.0 / 5.0 * 
        static_cast<double>(std::pow(static_cast<long double>(Ko), 3))) 
            * (1.0 + 1.0 / sig_t));
#endif
    const double pi0_2 = 1.0 / 3.0;
    const double pi0_3 = sig_t;
    const double pi0_4 = pi0_1;
    const double pi0_5 = sig_t;
 
    if (Ri > zero && Rrho > zero) {
      out_pi[0] = pi0_1 / (1.0 + Ri * Rrho / (a + Rrho));
      out_pi[1] = pi0_2 / (1.0 + Ri) * (1.0 + 2.0 * Ri * Rrho 
          / (1.0 + Rrho * Rrho));
      out_pi[2] = pi0_3;
      out_pi[3] = pi0_4 / (1.0 + Ri / (1.0 + a * Rrho));
      out_pi[4] = pi0_5;
    }
 
    if (Ri > zero && Rrho <= zero) {
      out_pi[0] = pi0_1 / (1.0 + Ri);
      out_pi[1] = pi0_2;
      out_pi[2] = pi0_3;
      out_pi[3] = pi0_4 / (1.0 + Ri);
      out_pi[4] = pi0_5;
    }
 
    if (Ri <= zero) {
      out_pi[0] = pi0_1;
      out_pi[1] = pi0_2;
      out_pi[2] = pi0_3;
      out_pi[3] = pi0_4;
      out_pi[4] = pi0_5;
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void struct_function (const double &Ri, const double &Rrho, const double *pi_n, 
      const double &Gm, double &struct_m, double &struct_h, double &struct_s, 
          double &struct_rho, double &r_out, double &r_out_new) const {
 
    // hwy_need_check
#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
    const double x = Ri * Gm * std::pow(1.0 - Rrho, -1);
#else
    const double x = Ri * Gm * 
        static_cast<double>(std::pow(static_cast<long double>(1.0 - Rrho), -1));
#endif
    const double p = pi_n[3] * pi_n[4] - pi_n[3] * pi_n[1] * (1.0 + Rrho);
    const double q = pi_n[0] * pi_n[1] * (1.0 + Rrho) - pi_n[0] * pi_n[2] * Rrho;
 
    // hwy_need_check
    const double r     = (pi_n[3] / pi_n[0]) / Rrho * (1.0 + q * x) / (1.0 + p * x);  
    // hwy_need_check
    const double r_new = (pi_n[3] / pi_n[0])        * (1.0 + q * x) / (1.0 + p * x);
    r_out     = r;
    r_out_new = r_new;
 
    const double Ah = pi_n[3] / (1.0 + p * x + pi_n[3] * pi_n[1] * x * (1.0 - 1.0 / r));
    const double As = Ah / r_new;
 
    const double Am1 = 4.0 / 5.0 - (pi_n[3] - pi_n[0] + (pi_n[0] - 1.0 / 150.0) 
        * (1.0 - 1.0 / r)) * x * Ah;
    const double Am2 = 10.0 + (pi_n[3] - pi_n[0] * Rrho) * x + 1.0 / 50.0 * Gm;
    const double Am = Am1 / Am2;
 
    const double LX = (1.0 - 1.0 / r) * x * Ah;
 
    const double ratio = 2.0 / 3.0 / (1.0 + 2.0 / 15.0 * LX + 1.0 / 10.0 * Am * Gm);
 
    struct_m = ratio * Am1 / Am2;
    struct_h = ratio * Ah;
    struct_s = ratio * As;
    struct_rho = (struct_h - struct_s * Rrho) / (1.0 - Rrho);
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void find_mix_layer_depth (const int& j, const int &i,
      double &mld_out, int &mld_lev_out, const ViewDouble3D &v_den, 
          const double &delta, const ViewDouble3D &v_zlev_in, const int &n) const {
    int k;
    for (k = 0; k < n; ++k) {
      v_zlev_(j, i, k) = std::abs(v_zlev_in(j, i, k));
    }
    for (k = 1; k < n; ++k) {
      if (std::abs(v_den(j, i, k) - v_den(j, i, 0)) > delta) {
        const double den_m = v_den(j, i, 0) - sign(delta, v_den(j, i, 0) - v_den(j, i, k));
        mld_out = v_zlev_(j, i, k) + (v_zlev_(j, i, k-1) - v_zlev_(j, i, k)) * (den_m - v_den(j, i, k)) 
            / (v_den(j, i, k-1) - v_den(j, i, k) + VERY_SMALL_);
        mld_lev_out = k + 1;
        break;
      }
    }
    if (k == n) {
      mld_out     = v_zlev_(j, i, n - 1);
      mld_lev_out = n;
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void mixing_efficiency (double &mix_eff_var, const double &struct_var, 
      const double &Ri, const double &Gm) const {
    mix_eff_var = 0.5 * Ri * Gm * struct_var;
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void canuto_2010_interface (const int &j, const int &i,
    const ViewDouble3D &v_Km_out,  // wk1
    const ViewDouble3D &v_Kh_out,  // wk2
    const ViewDouble3D &v_Ks_out,  // wk3
    const ViewDouble3D &v_Kd_out,  // wk4
          double       &mld_out,   // amld
    const ViewDouble3D &v_ts_in,   // wp1
    const ViewDouble3D &v_ss_in,   // wp2
    const ViewDouble3D &v_rho_in,  // wp3
    const ViewDouble3D &v_ri_in,   // wp4
    const ViewDouble3D &v_rrho_in, // wp5
    const ViewDouble3D &v_n2_in,   // wp7
    const ViewDouble3D &v_s2_in,   // wp6
    const double       &lat_in,    // ulat[iblock][j][i] / DegToRad
    const ViewDouble3D &v_lev_in,  // wp8
    const int          &num_lev)   /* kmt[iblock][j][i] */ const {

    int mld_lev = 0;
    find_mix_layer_depth (j, i, mld_out, mld_lev, 
        v_rho_in, 0.03, v_lev_in, num_lev);
 
    for (int k = 0; k < KM-1; ++k) {
      v_Ri_(j, i, k)   = v_ri_in(j, i, k);
      v_Rrho_(j, i, k) = v_rrho_in(j, i, k);
    }
 
    const double Ri_low  = - 1.0e+10;
    const double Ri_high =   1.0e+10;
 
    for (int k = 0; k < KM - 1; ++k) {
      if (v_Ri_(j, i, k) > Ri_high) {
        v_Ri_(j, i, k) = Ri_high;
      } else if (v_Ri_(j, i, k) < Ri_low) {
        v_Ri_(j, i, k) = Ri_low;
      }
    }
 
    const double Rrho_bound = 1.0e-3;
    for (int k = 0; k < num_lev - 1; ++k) {
      if (std::abs(v_Rrho_(j, i, k) - 1.0) < Rrho_bound) {
        if (v_Rrho_(j, i, k) >= 1.0) {
          v_Rrho_(j, i, k) = 1.001;
        } else if (v_Rrho_(j, i, k) <= 1.0) {
          v_Rrho_(j, i, k) = 0.999;
        }
      }
    }
    double out_pi[5];
    double cube[4];
    double struct_m_inf(0.0), struct_h_inf(0.0), struct_s_inf(0.0), struct_rho_inf(0.0);
    double R_inf(0.0), Rnew_inf(0.0);
    double Rf_inf(0.0);
    double struct_m(0.0), struct_h(0.0), struct_s(0.0), struct_rho(0.0);
    double R(0.0), Rnew(0.0);
    double mix_eff_m(0.0), mix_eff_h(0.0), mix_eff_s(0.0), mix_eff_rho(0.0);
    double Rf(0.0);
    double TKE_mld(0.0);
    for (int k = 0; k < num_lev - 1; ++k) {
      prepare_pi (1.0e10, v_Rrho_(j, i, k), out_pi);
 
      dyn_time_scale_calc (1.0e10, v_Rrho_(j, i, k), out_pi, v_Gm_(j, i, k), cube);
 
      struct_function(1.0e10, v_Rrho_(j, i, k), out_pi, v_Gm_(j, i, k), 
          struct_m_inf, struct_h_inf, struct_s_inf, struct_rho_inf, R_inf, Rnew_inf);
 
      Rf_calc(Rf_inf, 1.0e10, v_Rrho_(j, i, k), struct_h_inf, struct_m_inf, Rnew_inf);
 
      prepare_pi(v_Ri_(j, i, k), v_Rrho_(j, i, k), out_pi);
 
      dyn_time_scale_calc(v_Ri_(j, i, k), v_Rrho_(j, i, k), out_pi, v_Gm_(j, i, k), cube);
 
      struct_function(v_Ri_(j, i, k), v_Rrho_(j, i, k), out_pi, v_Gm_(j, i, k), 
          struct_m, struct_h, struct_s, struct_rho, R, Rnew);
    
      mixing_efficiency(mix_eff_m,   struct_m,   v_Ri_(j, i, k), v_Gm_(j, i, k));
      mixing_efficiency(mix_eff_h,   struct_h,   v_Ri_(j, i, k), v_Gm_(j, i, k));
      mixing_efficiency(mix_eff_s,   struct_s,   v_Ri_(j, i, k), v_Gm_(j, i, k));
      mixing_efficiency(mix_eff_rho, struct_rho, v_Ri_(j, i, k), v_Gm_(j, i, k));
 
      Rf_calc(Rf, v_Ri_(j, i, k), v_Rrho_(j, i, k), struct_h, struct_m, Rnew);
 
      if (k <= (mld_lev - 1)) {
        mixed_layer_TKE_calc (TKE_mld, mld_out, v_Gm_(j, i, k), 
            v_s2_in(j, i, k), v_lev_in(j, i, k), Rf, Rf_inf);
      
        v_Km_out(j, i, k) = mix_eff_m   * TKE_mld / (v_n2_in(j, i, k) + VERY_SMALL_);
        v_Kh_out(j, i, k) = mix_eff_h   * TKE_mld / (v_n2_in(j, i, k) + VERY_SMALL_);
        v_Ks_out(j, i, k) = mix_eff_s   * TKE_mld / (v_n2_in(j, i, k) + VERY_SMALL_);
        v_Kd_out(j, i, k) = mix_eff_rho * TKE_mld / (v_n2_in(j, i, k) + VERY_SMALL_);
      }
 
      if (k > (mld_lev - 1)) {
        thermocline_mixing_coeff_calc(
            v_Km_out(j, i, k), v_Kh_out(j, i, k), v_Ks_out(j, i, k), v_Kd_out(j, i, k), 
            mix_eff_m, mix_eff_h, mix_eff_s, mix_eff_rho, v_n2_in(j, i, k), lat_in);
      }
 
      if (k < 3) {
          v_Km_out(j, i, k) = std::max(1.0e-3, v_Km_out(j, i, k));
          v_Kh_out(j, i, k) = std::max(1.0e-3, v_Kh_out(j, i, k));
          v_Ks_out(j, i, k) = std::max(1.0e-3, v_Ks_out(j, i, k));
      } else {
          v_Km_out(j, i, k) = std::max(1.0e-4, v_Km_out(j, i, k));
          v_Kh_out(j, i, k) = std::max(1.0e-5, v_Kh_out(j, i, k));
          v_Ks_out(j, i, k) = std::max(1.0e-5, v_Ks_out(j, i, k));
      }
 
      v_Km_out(j, i, k) = std::min(1.2e-1, v_Km_out(j, i, k));
      v_Kh_out(j, i, k) = std::min(1.2e-1, v_Kh_out(j, i, k));
      v_Ks_out(j, i, k) = std::min(1.2e-1, v_Ks_out(j, i, k));
 
      if (v_n2_in(j, i, k) < VERY_SMALL_) {
        const double tmp = 8.0e-3 * v_dzp_(k) * v_dzp_(k);
        v_Kh_out(j, i, k) = std::min(tmp, 8.0);
        v_Ks_out(j, i, k) = std::min(tmp, 8.0);
        v_Km_out(j, i, k) = std::min(tmp, 8.0);
      }
    }
    return;
  }
};
#endif // (defined CANUTO) || (defined CANUTO2010)

// calculate the vertical mixing on U-grid
class FunctorReadyc6 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    if (k < KMM1) {
      tgrid_to_ugrid(iblock, k, j, i, v_akmu_, v_akmt_);
      v_akmu_(iblock, k, j, i) *= v_viv_(iblock, k+1, j, i);
    } else {
      v_akmu_(iblock, k, j, i) = C0;
    }
   return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (const int &iblock, 
			const int &k, const int &j, const int &i,
          const ViewDouble4D &v_ugrid, const ViewDouble4D &v_tgrid) 
              const {
    if (i >= 1 && j <(NY_BLOCK-1)) {
      // v_ugrid(iblock, k, j, i) = 
      //     v_au0_ (iblock, j, i) * v_tgrid(iblock, k, j  , i  ) +
      //     v_aus_ (iblock, j, i) * v_tgrid(iblock, k, j+1, i  ) +
      //     v_auw_ (iblock, j, i) * v_tgrid(iblock, k, j  , i-1) +
      //     v_ausw_(iblock, j, i) * v_tgrid(iblock, k, j+1, i-1);
      // v_ugrid(iblock, k, j, i) = P25 * v_tgrid(iblock, k, j  , i  ) +
      //                            P25 * v_tgrid(iblock, k, j+1, i  ) +
      //                            P25 * v_tgrid(iblock, k, j  , i-1) +
      //                            P25 * v_tgrid(iblock, k, j+1, i-1);
      v_ugrid(iblock, k, j, i) = P25 
          * (v_tgrid(iblock, k, j  , i  ) 
           + v_tgrid(iblock, k, j+1, i  ) 
           + v_tgrid(iblock, k, j  , i-1) 
           + v_tgrid(iblock, k, j+1, i-1));
    }
    if (i == 0 || j == (NY_BLOCK-1)) {
      v_ugrid(iblock, k, j, i) = 0;
    }
    return ;
  }
 private:
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
  const ViewDouble4D v_akmt_ = *p_v_akmt;
};
//-----------------------------------------------------------------
// COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
//-----------------------------------------------------------------
class FunctorReadyc7 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    upwell_1(iblock, j, i, v_h0_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_1 (const int &iblock,
      const int &j, const int &i, const ViewDouble3D &v_h0wk)
          const {
    tgrid_to_ugrid(iblock, j, i, v_work_, v_h0wk);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
      const int &iblock, const int &j, const int &i,
          const ViewDouble3D &v_ugrid, const ViewDouble3D &v_tgrid) 
              const {
    if (i >= 1 && j < (NY_BLOCK-1)) {
      // v_ugrid(iblock, j, i) = 
      //     v_au0_ (iblock, j, i) * v_tgrid(iblock, j  , i  ) +
      //     v_aus_ (iblock, j, i) * v_tgrid(iblock, j+1, i  ) +
      //     v_auw_ (iblock, j, i) * v_tgrid(iblock, j  , i-1) +
      //     v_ausw_(iblock, j, i) * v_tgrid(iblock, j+1, i-1);
      // v_ugrid(iblock, j, i) = P25 * v_tgrid(iblock, j  , i  ) +
      //                         P25 * v_tgrid(iblock, j+1, i  ) +
      //                         P25 * v_tgrid(iblock, j  , i-1) +
      //                         P25 * v_tgrid(iblock, j+1, i-1);
      v_ugrid(iblock, j, i) = P25 
          * (v_tgrid(iblock, j  , i  ) 
           + v_tgrid(iblock, j+1, i  ) 
           + v_tgrid(iblock, j  , i-1) 
           + v_tgrid(iblock, j+1, i-1));
    }
    if (i == 0 || j == (NY_BLOCK-1)) {
      v_ugrid(iblock, j, i) = 0;
    }
    return ;
  }
 private:
  const ViewDouble3D v_h0_   = *p_v_h0;
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble3D v_work_ = *p_v_work;
};

class FunctorReadyc8 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    upwell_2(iblock, k, j, i, v_u_, v_v_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_2 (const int &iblock,
      const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uwk, const ViewDouble4D &v_vwk)
              const {
    if (i >= 1 && j < (JMT-1)) {
      v_uk_(iblock, k, j, i) = (1.0 + v_work_(iblock, j, i) 
          * v_ohbu_(iblock, j, i)) * v_uwk(iblock, k, j, i);
      v_vk_(iblock, k, j, i) = (1.0 + v_work_(iblock, j, i) 
          * v_ohbu_(iblock, j, i)) * v_vwk(iblock, k, j, i);
    } else {
      v_uk_(iblock, k, j, i) = 0.0;
      v_vk_(iblock, k, j, i) = 0.0;
    }
    return ;
  }
 private:
  const ViewDouble3D v_ohbu_ = *p_v_ohbu;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_u_    = *p_v_u;
  const ViewDouble4D v_v_    = *p_v_v;
  const ViewDouble4D v_uk_   = *p_v_uk;
  const ViewDouble4D v_vk_   = *p_v_vk;
};
class FunctorReadyc9 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    upwell_3(iblock, k, j, i);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_3 (const int &iblock,
      const int &k, const int &j, const int &i) const {
    div(iblock, k, j, i, v_wka_, v_uk_, v_vk_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void div (const int &iblock, const int &k, const int &j, 
      const int &i, const ViewDouble4D &v_div_out, const ViewDouble4D &v_ux, 
          const ViewDouble4D &v_uy) const {
    v_div_out(iblock, k, j, i) = C0;
    if (i < (NX_BLOCK-1) && j >= 1) {
      const int bid = 0;
      if (k <= (v_kmt_(bid, j, i) - 1)) {
        v_div_out(iblock, k, j, i) = P5 * (
            (v_ux(iblock, k, j  , i+1) + v_ux(iblock, k, j-1, i+1)) * v_htw_(bid, j  , i+1)
          - (v_ux(iblock, k, j  , i  ) + v_ux(iblock, k, j-1, i  )) * v_htw_(bid, j  , i  )
          + (v_uy(iblock, k, j  , i+1) + v_uy(iblock, k, j  , i  )) * v_hts_(bid, j  , i  )
          - (v_uy(iblock, k, j-1, i+1) + v_uy(iblock, k, j-1, i  )) * v_hts_(bid, j-1, i  )) 
              * v_tarea_r_(bid, j, i);
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_uk_      = *p_v_uk;
  const ViewDouble4D v_vk_      = *p_v_vk;
  const ViewDouble4D v_wka_     = *p_v_wka;
};

class FunctorReadyc10 {
 public:
  KOKKOS_INLINE_FUNCTION void operator() (
      const int &j, const int &i) const {
    upwell_4(j, i, v_h0_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_4 (const int &j,
      const int &i, const ViewDouble3D &v_h0wk) const {
    const int iblock = 0;
    double work = C0;

    if (i >= 1 && i < (IMT-1) && j >= 1 && j < (JMT-1)) {
      for (int k = 0; k < KM; ++k) {
        work -= v_dzp_(k) * v_wka_(iblock, k, j, i)
                * v_vit_(iblock, k, j, i);
      }
      double ws = v_ws_(iblock, 0, j, i);
      for (int k = 1; k < KM; ++k) {
        ws = v_vit_(iblock, k, j, i) * (ws + v_dzp_(k-1) 
            * (work * v_ohbt_(iblock, j, i)
                + v_wka_(iblock, k-1, j, i)));
        v_ws_(iblock, k, j, i) = ws;
      }
      
      work = 1.0 / (1.0 + v_h0wk(iblock, j, i) * v_ohbt_(iblock, j, i));

      for (int k = 1; k < KM; ++k) {
        v_ws_(iblock, k, j, i) *= work;
      }
    }
    return ;
  }
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_ohbt_ = *p_v_ohbt;
  const ViewDouble4D v_ws_   = *p_v_ws;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

class FunctorReadyc11 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    tgrid_to_ugrid(iblock, k, j, i, v_wka_, v_ws_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_ugrid, const ViewDouble4D &v_tgrid) 
              const {
    if (i >= 1 && j < (NY_BLOCK-1)) {
      // v_ugrid(iblock, k, j, i) = 
      //     v_au0_ (iblock, j, i) * v_tgrid(iblock, k, j  , i  ) +
      //     v_aus_ (iblock, j, i) * v_tgrid(iblock, k, j+1, i  ) +
      //     v_auw_ (iblock, j, i) * v_tgrid(iblock, k, j  , i-1) +
      //     v_ausw_(iblock, j, i) * v_tgrid(iblock, k, j+1, i-1);
      // v_ugrid(iblock, k, j, i) = P25 * v_tgrid(iblock, k, j  , i  ) +
      //                            P25 * v_tgrid(iblock, k, j+1, i  ) +
      //                            P25 * v_tgrid(iblock, k, j  , i-1) +
      //                            P25 * v_tgrid(iblock, k, j+1, i-1);
      v_ugrid(iblock, k, j, i) = P25 
          * (v_tgrid(iblock, k, j  , i  ) 
           + v_tgrid(iblock, k, j+1, i  ) 
           + v_tgrid(iblock, k, j  , i-1) 
           + v_tgrid(iblock, k, j+1, i-1));
    }
    if (i == 0 || j == (NY_BLOCK-1)) {
      v_ugrid(iblock, k, j, i) = 0;
    }
    return ;
  }
 private:
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble4D v_ws_   = *p_v_ws;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

//-----------------------------------------------------------------
// COMPUTE THE ADVECTIVE TERMS
//-----------------------------------------------------------------
// advection_momentum(u, v, wka, dlu, dlv, iblock)
class FuncAdvMomCen1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advcetion_momentum_centered_1 (iblock, k, j, i,
        v_u_, v_v_, v_dlu_, v_dlv_);
    return ;
  };
  KOKKOS_INLINE_FUNCTION void advcetion_momentum_centered_1 (
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uuu, const ViewDouble4D &v_vvv,
              const ViewDouble4D &v_adv_uu, const ViewDouble4D &v_adv_vv) 
                  const {
    v_adv_uu(iblock, k, j, i) = C0;
    v_adv_vv(iblock, k, j, i) = C0;
    if (i > 0 && j < (JMT-1)) {
      v_uv_ws_face_(k, j, i, 0) = (
          v_uuu(iblock, k, j, i-1) + v_uuu(iblock, k, j  , i))
              * P25 * v_hue_(iblock, j  , i-1);
      v_uv_ws_face_(k, j, i, 1) = (
          v_vvv(iblock, k, j, i  ) + v_vvv(iblock, k, j+1, i))
              * P25 * v_hun_(iblock, j+1, i  );
    }
    return ;
  }
 private:
  const ViewDouble3D v_hun_        = *p_v_hun;
  const ViewDouble3D v_hue_        = *p_v_hue;
  const ViewDouble4D v_u_          = *p_v_u;
  const ViewDouble4D v_v_          = *p_v_v;
  const ViewDouble4D v_dlu_        = *p_v_dlu;
  const ViewDouble4D v_dlv_        = *p_v_dlv;
  const ViewDouble4D v_uv_ws_face_ = *p_v_uv_ws_face;
};
class FuncAdvMomFlu1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator() 
      (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    advcetion_momentum_flux_1(iblock, k, j, i,
        v_u_, v_v_, v_dlu_, v_dlv_);
    return ;
  };
  KOKKOS_INLINE_FUNCTION void advcetion_momentum_flux_1(
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uuu, const ViewDouble4D &v_vvv,
              const ViewDouble4D &v_adv_uu, const ViewDouble4D &v_adv_vv) 
                  const {
    v_adv_uu(iblock, k, j, i) = C0;
    v_adv_vv(iblock, k, j, i) = C0;
    if (i > 0 && j < (JMT-1)) {
      v_uv_ws_face_(k, j, i, 0) = 
          (v_uuu(iblock, k, j, i-1) * v_dyu_(iblock, j, i-1)  
         + v_uuu(iblock, k, j, i  ) * v_dyu_(iblock, j, i  )) * P25;
      v_uv_ws_face_(k, j, i, 1) = 
          (v_vvv(iblock, k, j  , i) * v_dxu_(iblock, j  , i)  
         + v_vvv(iblock, k, j+1, i) * v_dxu_(iblock, j+1, i)) * P25;
    }
    return ;
  }
 private:
  const ViewDouble3D v_dxu_        = *p_v_dxu;
  const ViewDouble3D v_dyu_        = *p_v_dyu;
  const ViewDouble4D v_u_          = *p_v_u;
  const ViewDouble4D v_v_          = *p_v_v;
  const ViewDouble4D v_dlu_        = *p_v_dlu;
  const ViewDouble4D v_dlv_        = *p_v_dlv;
  const ViewDouble4D v_uv_ws_face_ = *p_v_uv_ws_face;
};

class FuncAdvMomCen2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advcetion_momentum_centered_2 (iblock, k, j, i,
        v_u_, v_v_, v_wka_, v_dlu_, v_dlv_);
    return ;
  };
  KOKKOS_INLINE_FUNCTION void advcetion_momentum_centered_2 (
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uuu, const ViewDouble4D &v_vvv,
              const ViewDouble4D &v_www, const ViewDouble4D &v_adv_uu, 
                  const ViewDouble4D &v_adv_vv) const {
    double adv_z1, adv_z2, adv_z3, adv_z4;
    v_adv_uu(iblock, k, j, i) = (
        - v_uv_ws_face_(k, j  , i  , 0) 
            * (v_uuu(iblock, k, j  , i  ) - v_uuu(iblock, k, j  , i-1))
        - v_uv_ws_face_(k, j  , i+1, 0) 
            * (v_uuu(iblock, k, j  , i+1) - v_uuu(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j  , i  , 1) 
            * (v_uuu(iblock, k, j+1, i  ) - v_uuu(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j-1, i  , 1) 
            * (v_uuu(iblock, k, j  , i  ) - v_uuu(iblock, k, j-1, i  ))) 
                * v_uarea_r_(iblock, j, i);

    v_adv_vv(iblock, k, j, i) = (
        - v_uv_ws_face_(k, j  , i  , 0) 
            * (v_vvv(iblock, k, j  , i  ) - v_vvv(iblock, k, j  , i-1))
        - v_uv_ws_face_(k, j  , i+1, 0)               
            * (v_vvv(iblock, k, j  , i+1) - v_vvv(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j  , i  , 1)               
            * (v_vvv(iblock, k, j+1, i  ) - v_vvv(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j-1, i  , 1)               
            * (v_vvv(iblock, k, j  , i  ) - v_vvv(iblock, k, j-1, i  ))) 
                * v_uarea_r_(iblock, j, i);

    if (k == 0) {
      adv_z1 = 0.0;
      adv_z3 = 0.0;
    } else {
      adv_z1 = v_www(iblock, k, j, i) 
          * (v_uuu(iblock, k-1, j, i  ) - v_uuu(iblock, k, j, i));
      adv_z3 = v_www(iblock, k, j, i) 
          * (v_vvv(iblock, k-1, j, i  ) - v_vvv(iblock, k, j, i));
    }

    if (k == KM-1) {
      adv_z2 = 0.0;
      adv_z4 = 0.0;
    } else {
      adv_z2 = v_www(iblock, k+1, j, i) 
          * (v_uuu(iblock, k  , j  , i  ) - v_uuu(iblock, k+1, j, i));
      adv_z4 = v_wka_(iblock, k+1, j, i) 
          * (v_vvv(iblock, k  , j  , i  ) - v_vvv(iblock, k+1, j, i));
    }

    v_adv_uu(iblock, k, j, i) -= P5 * v_odzp_(k)
        * (adv_z1 + adv_z2);
    v_adv_vv(iblock, k, j, i) -= P5 * v_odzp_(k)
        * (adv_z3 + adv_z4);
    return ;
  }
 private:
  const ViewDouble1D v_odzp_    = *p_v_odzp;
  const ViewDouble3D v_uarea_r_ = *p_v_uarea_r;
  const ViewDouble4D v_u_       = *p_v_u;
  const ViewDouble4D v_v_       = *p_v_v;
  const ViewDouble4D v_dlu_     = *p_v_dlu;
  const ViewDouble4D v_dlv_     = *p_v_dlv;
  const ViewDouble4D v_wka_     = *p_v_wka;
  const ViewDouble4D v_uv_ws_face_ = *p_v_uv_ws_face;
};

class FuncAdvMomFlu2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator() 
      (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    advcetion_momentum_flux_2(iblock, k, j, i,
        v_u_, v_v_, v_wka_, v_dlu_, v_dlv_);
    return ;
  };
  KOKKOS_INLINE_FUNCTION void advcetion_momentum_flux_2(
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uuu, const ViewDouble4D &v_vvv,
              const ViewDouble4D &v_www, const ViewDouble4D &v_adv_uu, 
                  const ViewDouble4D &v_adv_vv) const {
    double adv_z1, adv_z2, adv_z3, adv_z4;
    v_adv_uu(iblock, k, j, i) = (
        - v_uv_ws_face_(k, j  , i  , 0) 
            * (v_uuu(iblock, k, j  , i  ) + v_uuu(iblock, k, j  , i-1))
        + v_uv_ws_face_(k, j  , i+1, 0) 
            * (v_uuu(iblock, k, j  , i+1) + v_uuu(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j-1, i  , 1) 
            * (v_uuu(iblock, k, j  , i  ) + v_uuu(iblock, k, j-1, i  ))
        + v_uv_ws_face_(k, j  , i  , 1) 
            * (v_uuu(iblock, k, j+1, i  ) + v_uuu(iblock, k, j  , i  ))) 
                * v_uarea_r_(iblock, j, i);

    v_adv_vv(iblock, k, j, i) = (
        - v_uv_ws_face_(k, j  , i  , 0) 
            * (v_vvv(iblock, k, j  , i  ) + v_vvv(iblock, k, j  , i-1))
        + v_uv_ws_face_(k, j  , i+1, 0)               
            * (v_vvv(iblock, k, j  , i+1) + v_vvv(iblock, k, j  , i  ))
        - v_uv_ws_face_(k, j-1, i  , 1)               
            * (v_vvv(iblock, k, j  , i  ) + v_vvv(iblock, k, j-1, i  ))
        + v_uv_ws_face_(k, j  , i  , 1)               
            * (v_vvv(iblock, k, j+1, i  ) + v_vvv(iblock, k, j  , i  ))) 
                * v_uarea_r_(iblock, j, i);

    if (k == 0) {
      adv_z1 = 0.0;
      adv_z3 = 0.0;
    } else {
      adv_z1 = v_www(iblock, k, j, i) 
          * (v_uuu(iblock, k-1, j, i  ) + v_uuu(iblock, k, j, i)) * P5;
      adv_z3 = v_www(iblock, k, j, i) 
          * (v_vvv(iblock, k-1, j, i  ) + v_vvv(iblock, k, j, i)) * P5;
    }

    if (k == KM-1) {
      adv_z2 = 0.0;
      adv_z4 = 0.0;
    } else {
      adv_z2 = v_www(iblock, k+1, j, i) 
          * (v_uuu(iblock, k, j, i) + v_uuu(iblock, k+1, j, i)) * P5;
      adv_z4 = v_www(iblock, k+1, j, i) 
          * (v_vvv(iblock, k, j, i) - v_vvv(iblock, k+1, j, i)) * P5;
    }

    v_adv_uu(iblock, k, j, i) -= v_odzp_(k) * (adv_z2 - adv_z1);
    v_adv_vv(iblock, k, j, i) -= v_odzp_(k) * (adv_z4 - adv_z3);
    return ;
  }
 private:
  const ViewDouble1D v_odzp_    = *p_v_odzp;
  const ViewDouble3D v_uarea_r_ = *p_v_uarea_r;
  const ViewDouble4D v_u_       = *p_v_u;
  const ViewDouble4D v_v_       = *p_v_v;
  const ViewDouble4D v_dlu_     = *p_v_dlu;
  const ViewDouble4D v_dlv_     = *p_v_dlv;
  const ViewDouble4D v_wka_     = *p_v_wka;
  const ViewDouble4D v_uv_ws_face_ = *p_v_uv_ws_face;
};
// End advection_momentum(u, v, wka, dlu, dlv, iblock)

class FunctorReadyc14 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    hdiffu_del4_1(k, j, i, v_up_, v_vp_);
#endif // BIHAR
    return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_1(
      const int &k, const int &j, const int &i,
          const ViewDouble4D &v_umixk, const ViewDouble4D &v_vmixk)
              const {
    div  (k, j, i, v_wka_, v_umixk, v_vmixk);
    zcurl(k, j, i, v_wkb_, v_umixk, v_vmixk);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void div(const int &k, const int &j, 
      const int &i, const ViewDouble4D &v_div_out, 
          const ViewDouble4D &v_ux, const ViewDouble4D &v_uy) const {
    const int iblock = 0;
    v_div_out(iblock, k, j, i) = C0;
    if (i < (NX_BLOCK-1) && j >= 1) {
      const int bid = 0;
      if (k <= (v_kmt_(bid, j, i) - 1)) {
        v_div_out(iblock, k, j, i) = P5 * (
            (v_ux(iblock, k, j  , i+1) + v_ux(iblock, k, j-1, i+1)) * v_htw_(bid, j  , i+1)
          - (v_ux(iblock, k, j  , i  ) + v_ux(iblock, k, j-1, i  )) * v_htw_(bid, j  , i  )
          + (v_uy(iblock, k, j  , i+1) + v_uy(iblock, k, j  , i  )) * v_hts_(bid, j  , i  )
          - (v_uy(iblock, k, j-1, i+1) + v_uy(iblock, k, j-1, i  )) * v_hts_(bid, j-1, i  )) 
              * v_tarea_r_(bid, j, i);
      }
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void zcurl (const int &k,
      const int &j, const int &i, const ViewDouble4D &v_curl, 
          const ViewDouble4D &v_ux, const ViewDouble4D &v_uy) 
              const {
    const int iblock = 0;
    v_curl(iblock, k, j, i) = C0;
    if (i >= 1 && j >= 1) {
      const int bid = 0;
      if (k <= v_kmt_(bid, j, i) - 1) {
        v_curl(iblock, k, j, i) = P5 * (
            v_uy(iblock, k, j  , i  ) * v_dyu_(bid, j  , i  )
          + v_uy(iblock, k, j-1, i  ) * v_dyu_(bid, j-1, i  )
          - v_uy(iblock, k, j  , i-1) * v_dyu_(bid, j  , i-1)
          - v_uy(iblock, k, j-1, i-1) * v_dyu_(bid, j-1, i-1)
          - v_ux(iblock, k, j  , i  ) * v_dxu_(bid, j  , i  )
          - v_ux(iblock, k, j  , i-1) * v_dxu_(bid, j  , i-1)
          + v_ux(iblock, k, j-1, i  ) * v_dxu_(bid, j-1, i  )
          + v_ux(iblock, k, j-1, i-1) * v_dxu_(bid, j-1, i-1));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble3D v_dxu_     = *p_v_dxu;
  const ViewDouble3D v_dyu_     = *p_v_dyu;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_up_      = *p_v_up;
  const ViewDouble4D v_vp_      = *p_v_vp;
  const ViewDouble4D v_wka_     = *p_v_wka;
  const ViewDouble4D v_wkb_     = *p_v_wkb;
#endif // BIHAR
};

class FunctorReadyc15 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    hdiffu_del4_2 (k, j, i, v_up_, v_vp_);
#endif // BIHAR
    return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_2 (
      const int &k, const int &j, const int &i,
          const ViewDouble4D &v_umixk, const ViewDouble4D &v_vmixk)
              const {

    if (i >= (ib_-2) && i < (ie_+1) && j >= (jb_-2) && j < (je_+1)) {
      const int bid = 0;
      const double cc = v_du_cnsewm_(bid, j, i, 0) 
          + v_du_cnsewm_(bid, j, i, 5);
      v_d2uk_(k, j, i) = (           cc * v_umixk(bid, k, j  , i  )
           + v_du_cnsewm_(bid, j, i, 1) * v_umixk(bid, k, j-1, i  )
           + v_du_cnsewm_(bid, j, i, 2) * v_umixk(bid, k, j+1, i  )
           + v_du_cnsewm_(bid, j, i, 3) * v_umixk(bid, k, j  , i+1)
           + v_du_cnsewm_(bid, j, i, 4) * v_umixk(bid, k, j  , i-1))
          + (v_dm_cnsew_ (bid, j, i, 0) * v_vmixk(bid, k, j  , i  )
           + v_dm_cnsew_ (bid, j, i, 1) * v_vmixk(bid, k, j-1, i  )
           + v_dm_cnsew_ (bid, j, i, 2) * v_vmixk(bid, k, j+1, i  )
           + v_dm_cnsew_ (bid, j, i, 3) * v_vmixk(bid, k, j  , i+1)
           + v_dm_cnsew_ (bid, j, i, 4) * v_vmixk(bid, k, j  , i-1));
      
      v_d2vk_(k, j, i) = (           cc * v_vmixk(bid, k, j  , i  )
           + v_du_cnsewm_(bid, j, i, 1) * v_vmixk(bid, k, j-1, i  )
           + v_du_cnsewm_(bid, j, i, 2) * v_vmixk(bid, k, j+1, i  )
           + v_du_cnsewm_(bid, j, i, 3) * v_vmixk(bid, k, j  , i+1)
           + v_du_cnsewm_(bid, j, i, 4) * v_vmixk(bid, k, j  , i-1))
          - (v_dm_cnsew_ (bid, j, i, 0) * v_umixk(bid, k, j  , i  )
           + v_dm_cnsew_ (bid, j, i, 1) * v_umixk(bid, k, j-1, i  )
           + v_dm_cnsew_ (bid, j, i, 2) * v_umixk(bid, k, j+1, i  )
           + v_dm_cnsew_ (bid, j, i, 3) * v_umixk(bid, k, j  , i+1)
           + v_dm_cnsew_ (bid, j, i, 4) * v_umixk(bid, k, j  , i-1));
    } else {
      v_d2uk_(k, j, i) = 0.0;
      v_d2vk_(k, j, i) = 0.0;
    }
    return ;
  }
#endif // BIHAR
 private:
#ifdef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const ViewDouble3D v_d2uk_      = *p_v_d2uk;
  const ViewDouble3D v_d2vk_      = *p_v_d2vk;
  const ViewDouble4D v_up_        = *p_v_up;
  const ViewDouble4D v_vp_        = *p_v_vp;
  const ViewDouble4D v_dm_cnsew_  = *p_v_dm_cnsew;
  const ViewDouble4D v_du_cnsewm_ = *p_v_du_cnsewm;
#endif // BIHAR
};
class FunctorReadyc16 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    const int iblock = 0;
    hdiffu_del4_3 (iblock, k, j, i);
#endif // BIHAR
    return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_3 (const int &iblock,
      const int &k, const int &j, const int &i) const {
    const int bid = 0;

    double am_factor = 1.0;
    const double amf = v_amf_(bid, j, i);
    if (i >= (ib_-2) && i < (ie_+1) && j >= (jb_-2) && j < (je_+1)) {
      double gradx1, grady1;
      grad (k, j, i, gradx1, grady1, v_wka_);
      gradx1 *= gradx1;
      grady1 *= grady1;
      double gradx2, grady2;
      grad (k, j, i, gradx2, grady2, v_wkb_);
      gradx2 *= gradx2;
      grady2 *= grady2;
      const double sqrt_uarea = sqrt (v_uarea_(bid, j, i));
      const double dxdy = sqrt_uarea * sqrt_uarea * sqrt_uarea * sqrt_uarea * sqrt_uarea * 45.0;
      am_factor = sqrt(gradx1 + gradx2 + grady1 + grady2)  
          * dxdy / fabs(am_ * amf);
    }
    am_factor = fmin(40.0, am_factor);
    am_factor = fmax(1.0,  am_factor);
    if (k <= v_kmu_(bid, j, i) - 1) {
      v_d2uk_(k, j, i) *= am_factor * amf;
      v_d2vk_(k, j, i) *= am_factor * amf;
    } else {
      v_d2uk_(k, j, i) = C0; 
      v_d2vk_(k, j, i) = C0;
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &k, const int &j,
      const int &i, double &gradx, double &grady, const ViewDouble4D &v_f) 
              const {
    const int bid = 0;
    gradx = 0.0;
    grady = 0.0;
    if (i >= 1 && j < (NY_BLOCK-1)) {
      if (k <= v_kmu_(bid, j, i) - 1) {
        gradx = v_dxyur_(bid, j, i, 0) * P5 * 
            (v_f(bid, k, j+1, i  ) - v_f(bid, k, j, i-1) - 
             v_f(bid, k, j+1, i-1) + v_f(bid, k, j, i  ));
      
        grady = v_dxyur_(bid, j, i, 1) * P5 * 
            (v_f(bid, k, j+1, i  ) - v_f(bid, k, j, i-1) + 
             v_f(bid, k, j+1, i-1) - v_f(bid, k, j, i  ));
      }
    }
    return ;
  }
#endif // BIHAR
 private:
#ifdef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double am_ = CppHmixDel4::am;
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble3D v_d2uk_  = *p_v_d2uk;
  const ViewDouble3D v_d2vk_  = *p_v_d2vk;
  const ViewDouble3D v_amf_   = *p_v_amf;
  const ViewDouble3D v_uarea_ = *p_v_uarea;
  const ViewDouble4D v_wka_   = *p_v_wka;
  const ViewDouble4D v_wkb_   = *p_v_wkb;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
#endif // BIHAR
};
class FunctorReadyc17 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    const int iblock = 0;
    double hduk, hdvk;
    hdiffu_del4_4 (iblock, k, j, i, hduk, hdvk);
    v_dlu_(iblock, k, j, i) += hduk;
    v_dlv_(iblock, k, j, i) += hdvk;
#endif // BIHAR
  return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_4 (
      const int &iblock, const int &k, const int &j, const int &i,
          double &hduk, double &hdvk)
                  const {
    const int bid = 0;
    hduk = C0;
    hdvk = C0;
    if (i >= (ib_-1) && i < (ie_) && j >= (jb_-1) && j < (je_)) {
      const double cc = v_du_cnsewm_(bid, j, i, 0) 
          + v_du_cnsewm_(bid, j, i, 5);
      hduk = am_ * ((               cc * v_d2uk_(k, j  , i  )
          + v_du_cnsewm_(bid, j, i, 1) * v_d2uk_(k, j-1, i  )
          + v_du_cnsewm_(bid, j, i, 2) * v_d2uk_(k, j+1, i  )
          + v_du_cnsewm_(bid, j, i, 3) * v_d2uk_(k, j  , i+1)
          + v_du_cnsewm_(bid, j, i, 4) * v_d2uk_(k, j  , i-1))
         + (v_dm_cnsew_( bid, j, i, 0) * v_d2vk_(k, j  , i  )
          + v_dm_cnsew_( bid, j, i, 1) * v_d2vk_(k, j-1, i  )
          + v_dm_cnsew_( bid, j, i, 2) * v_d2vk_(k, j+1, i  )
          + v_dm_cnsew_( bid, j, i, 3) * v_d2vk_(k, j  , i+1)
          + v_dm_cnsew_( bid, j, i, 4) * v_d2vk_(k, j  , i-1)));
   
      hdvk = am_ * ((               cc * v_d2vk_(k, j  , i  )
          + v_du_cnsewm_(bid, j, i, 1) * v_d2vk_(k, j-1, i  )
          + v_du_cnsewm_(bid, j, i, 2) * v_d2vk_(k, j+1, i  )
          + v_du_cnsewm_(bid, j, i, 3) * v_d2vk_(k, j  , i+1)
          + v_du_cnsewm_(bid, j, i, 4) * v_d2vk_(k, j  , i-1))
         - (v_dm_cnsew_( bid, j, i, 0) * v_d2uk_(k, j  , i  )
          + v_dm_cnsew_( bid, j, i, 1) * v_d2uk_(k, j-1, i  )
          + v_dm_cnsew_( bid, j, i, 2) * v_d2uk_(k, j+1, i  )
          + v_dm_cnsew_( bid, j, i, 3) * v_d2uk_(k, j  , i+1)
          + v_dm_cnsew_( bid, j, i, 4) * v_d2uk_(k, j  , i-1)));
    }
    if (k > v_kmu_(bid, j, i) - 1) {
      hduk = C0;
      hdvk = C0;
    }
    return ;
  }
#endif // BIHAR
 private:
#ifdef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double am_ = CppHmixDel4::am;
  const ViewInt3D    v_kmu_       = *p_v_kmu;
  const ViewDouble3D v_d2uk_      = *p_v_d2uk;
  const ViewDouble3D v_d2vk_      = *p_v_d2vk;
  const ViewDouble4D v_dlu_       = *p_v_dlu;
  const ViewDouble4D v_dlv_       = *p_v_dlv;
  const ViewDouble4D v_dm_cnsew_  = *p_v_dm_cnsew;
  const ViewDouble4D v_du_cnsewm_ = *p_v_du_cnsewm;
#endif // BIHAR
};

class FuncReadyc15Del2 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
#ifndef BIHAR
    const int iblock = 0;
    double hduk, hdvk;
    hdiffu_del2 (k, j, i, hduk, hdvk, v_up_, v_vp_);
    v_dlu_(iblock, k, j, i) += hduk;
    v_dlv_(iblock, k, j, i) += hdvk;
#endif // BIHAR
  return ;
  }
#ifndef BIHAR
  KOKKOS_INLINE_FUNCTION 
  void hdiffu_del2 (const int &k, const int &j, const int &i,
      double &hduk, double &hdvk, 
      const ViewDouble4D &v_umixk, const ViewDouble4D &v_vmixk) const {
    const int bid = 0;
    hduk = C0;
    hdvk = C0;

    if (i >= (ib_-1) && i < (ie_) && j >= (jb_-1) && j < (je_)) {
      const double cc = v_duc_(bid, j, i) + v_dum_(bid, j, i);
      hduk = am_ * ((      cc * v_umixk(bid, k, j  , i  )
          + v_dun_(bid, j, i) * v_umixk(bid, k, j-1, i  )
          + v_dus_(bid, j, i) * v_umixk(bid, k, j+1, i  )
          + v_due_(bid, j, i) * v_umixk(bid, k, j  , i+1)
          + v_duw_(bid, j, i) * v_umixk(bid, k, j  , i-1))
         + (v_dmc_(bid, j, i) * v_vmixk(bid, k, j  , i  )
          + v_dmn_(bid, j, i) * v_vmixk(bid, k, j-1, i  )
          + v_dms_(bid, j, i) * v_vmixk(bid, k, j+1, i  )
          + v_dme_(bid, j, i) * v_vmixk(bid, k, j  , i+1)
          + v_dmw_(bid, j, i) * v_vmixk(bid, k, j  , i-1))) 
              * v_viv_(bid, k, j, i);
   
      hdvk = am_ * ((      cc * v_vmixk(bid, k, j  , i  )
          + v_dun_(bid, j, i) * v_vmixk(bid, k, j-1, i  )
          + v_dus_(bid, j, i) * v_vmixk(bid, k, j+1, i  )
          + v_due_(bid, j, i) * v_vmixk(bid, k, j  , i+1)
          + v_duw_(bid, j, i) * v_vmixk(bid, k, j  , i-1))
         - (v_dmc_(bid, j, i) * v_umixk(bid, k, j  , i  )
          + v_dmn_(bid, j, i) * v_umixk(bid, k, j-1, i  )
          + v_dms_(bid, j, i) * v_umixk(bid, k, j+1, i  )
          + v_dme_(bid, j, i) * v_umixk(bid, k, j  , i+1)
          + v_dmw_(bid, j, i) * v_umixk(bid, k, j  , i-1)))
              * v_viv_(bid, k, j, i);
    }
    if (k > v_kmu_(bid, j, i) - 1) {
      hduk = C0;
      hdvk = C0;
    }
    return ;
  }
#endif // BIHAR
 private:
#ifndef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double am_ = CppHmixDel2::am;
  const ViewInt3D    v_kmu_ = *p_v_kmu;
  const ViewDouble3D v_duc_ = *p_v_duc;
  const ViewDouble3D v_dum_ = *p_v_dum;
  const ViewDouble3D v_dun_ = *p_v_dun;
  const ViewDouble3D v_dus_ = *p_v_dus;
  const ViewDouble3D v_due_ = *p_v_due;
  const ViewDouble3D v_duw_ = *p_v_duw;
  const ViewDouble3D v_dmc_ = *p_v_dmc;
  const ViewDouble3D v_dmn_ = *p_v_dmn;
  const ViewDouble3D v_dms_ = *p_v_dms;
  const ViewDouble3D v_dme_ = *p_v_dme;
  const ViewDouble3D v_dmw_ = *p_v_dmw;
  const ViewDouble4D v_up_  = *p_v_up;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_dlu_ = *p_v_dlu;
  const ViewDouble4D v_dlv_ = *p_v_dlv;
  const ViewDouble4D v_viv_ = *p_v_viv;
#endif // BIHAR
};
//-----------------------------------------------------------------
// VERTICAL INTEGRATION
//-----------------------------------------------------------------
class FunctorReadyc18 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    vinteg(j, i, v_dlu_, v_dlub_);
    vinteg(j, i, v_dlv_, v_dlvb_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void vinteg(const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble3D &v_wk2) 
          const {
    const int iblock = 0;
    double wk2 = C0;
    for (int k = 0; k < KM; ++k) {
      wk2 += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(iblock, k, j, i) *v_viv_(iblock, k, j, i);
    }
    v_wk2(iblock, j, i) = wk2;
    return ;
  }
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_ohbu_ = *p_v_ohbu;
  const ViewDouble3D v_dlub_ = *p_v_dlub;
  const ViewDouble3D v_dlvb_ = *p_v_dlvb;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
  const ViewDouble4D v_viv_  = *p_v_viv;
};

class FunctorReadyc19 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = c0f_ * std::sqrt(
        v_up_(iblock, k, j, i) * v_up_(iblock, k, j, i)
      + v_vp_(iblock, k, j, i) * v_vp_(iblock, k, j, i));
    return ;
  }
 private:
  const double c0f_ = CppPconstMod::c0f;
  const ViewDouble4D v_up_  = *p_v_up;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_wka_ = *p_v_wka;
};

class FunctorReadyc20 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    double sbcx(0.0), sbcy(0.0), bbcx(0.0), bbcy(0.0);

    const int kmb = v_kmu_(iblock, j, i) - 1;

    if (kmb >= 0) {
      sbcx = v_su_(iblock, j, i) * od0_;
      sbcy = v_sv_(iblock, j, i) * od0_;
      bbcx = v_wka_(iblock, kmb, j, i) * (
          v_up_(iblock, kmb, j, i) * cag_ + v_snlat_(iblock, j, i) 
        * v_vp_(iblock, kmb, j, i) * sag_);
         
      bbcy = v_wka_(iblock, kmb, j, i) * (- v_snlat_(iblock, j, i) 
          * v_up_(iblock, kmb, j, i) * sag_ 
          + v_vp_(iblock, kmb, j, i) * cag_);
    }
    v_dlub_(iblock, j, i) += (sbcx - bbcx) * v_ohbu_(iblock, j, i);
    v_dlvb_(iblock, j, i) += (sbcy - bbcy) * v_ohbu_(iblock, j, i);
    return ;
  }
 private:
  const double od0_ = CppPconstMod::od0;
  const double c0f_ = CppPconstMod::c0f;
  const double cag_ = CppPconstMod::cag;
  const double sag_ = CppPconstMod::sag;
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble3D v_su_    = *p_v_su;
  const ViewDouble3D v_sv_    = *p_v_sv;
  const ViewDouble3D v_snlat_ = *p_v_snlat;
  const ViewDouble4D v_up_    = *p_v_up;
  const ViewDouble4D v_vp_    = *p_v_vp;
  const ViewDouble4D v_wka_   = *p_v_wka;
  const ViewDouble3D v_ohbu_  = *p_v_ohbu;
  const ViewDouble3D v_dlub_  = *p_v_dlub;
  const ViewDouble3D v_dlvb_  = *p_v_dlvb;
};

class FunctorReadyc21 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    // d0 = 1026.0D0, od0 = 1 / d0
    const double od0 = 1.0 / 1026.0;
    // const double aidif = 0.5;
    // aidifm1 = 1.0 - aidif
    const double aidifm1 = 0.5;
    double diff_u1 = 0.0;
    double diff_u2 = 0.0;
    if (k == 0) {
      diff_u1 = v_su_(iblock, j, i) * od0 * aidifm1;
    } else {
      diff_u1 = v_akmu_(iblock, k-1, j, i) * aidifm1 
          * (v_up_(iblock, k-1, j, i) - v_up_(iblock, k, j, i)) 
          * v_odz_pt_(k, 1) * v_viv_(iblock, k, j, i) 
          + (1.0 - v_viv_(iblock, k, j, i)) 
          * v_wka_(iblock, k-1, j, i) * aidifm1 
          * (v_up_(iblock, k-1, j, i) * cag_ + v_snlat_(iblock, j, i) 
          * v_vp_(iblock, k-1, j, i) * sag_);
    }
    if (k == KM-1) {
      diff_u2 = v_wka_(iblock, k, j, i) * (v_up_(iblock, k, j, i) * cag_ 
          + v_snlat_(iblock, j, i) * v_vp_(iblock, k, j, i) * sag_) 
          * aidifm1;
    } else {
      diff_u2 = v_akmu_(iblock, k, j, i) * aidifm1 
          * (v_up_(iblock, k, j, i) - v_up_(iblock, k+1, j, i)) 
          * v_odz_pt_(k+1, 1) * v_viv_(iblock, k+1, j, i) 
          + (1.0 - v_viv_(iblock, k+1, j, i)) * v_wka_(iblock, k, j, i) 
          * aidifm1 * (v_up_(iblock, k, j, i) * cag_ 
          + v_snlat_(iblock, j, i) * v_vp_(iblock, k, j, i) * sag_);
    }
    v_dlu_(iblock, k, j, i) += v_odz_pt_(k, 0) * (diff_u1 - diff_u2);
    return ;
  }
 private:
  const double cag_ = CppPconstMod::cag;
  const double sag_ = CppPconstMod::sag;
  const ViewDouble2D v_odz_pt_ = *p_v_odz_pt;
  const ViewDouble3D v_su_     = *p_v_su;
  const ViewDouble3D v_snlat_  = *p_v_snlat;
  const ViewDouble4D v_up_     = *p_v_up;
  const ViewDouble4D v_vp_     = *p_v_vp;
  const ViewDouble4D v_dlu_    = *p_v_dlu;
  const ViewDouble4D v_viv_    = *p_v_viv;
  const ViewDouble4D v_wka_    = *p_v_wka;
  const ViewDouble4D v_akmu_   = *p_v_akmu;
};
class FunctorReadyc22 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    const double od0 = 1.0 / 1026.0;
    // const double aidif = 0.5;
    // aidifm1 = 1.0 - aidif
    const double aidifm1 = 0.5;
    double diff_v1 = 0.0;
    double diff_v2 = 0.0;
    if (k == 0) {
      diff_v1 = v_sv_(iblock, j, i) * od0 * aidifm1;
    } else {
      diff_v1 = v_akmu_(iblock, k-1, j, i) * aidifm1 
          * (v_vp_(iblock, k-1, j, i) - v_vp_(iblock, k, j, i)) 
          * v_odz_pt_(k, 1) * v_viv_(iblock, k, j, i) 
          + (1.0 - v_viv_(iblock, k, j, i)) 
          * v_wka_(iblock, k-1, j, i) * aidifm1 
          * (- v_snlat_(iblock, j, i) * v_up_(iblock, k-1, j, i) 
          * sag_ + v_vp_(iblock, k-1, j, i) * cag_);
    }
    if (k == KM-1) {
      diff_v2 = v_wka_(iblock, k, j, i) * (- v_snlat_(iblock, j, i) 
          * v_up_(iblock, k, j, i) * sag_ 
          + v_vp_(iblock, k, j, i) * cag_) * aidifm1;
    } else {
      diff_v2 = v_akmu_(iblock, k, j, i) * aidifm1 * 
          (v_vp_(iblock, k, j, i) - v_vp_(iblock, k+1, j, i)) 
          * v_odz_pt_(k+1, 1) * v_viv_(iblock, k+1, j, i) 
          + (1.0 - v_viv_(iblock, k+1, j, i)) * v_wka_(iblock, k, j, i) 
          * aidifm1 * (-v_snlat_(iblock, j, i) 
          * v_up_(iblock, k, j, i) * sag_ + v_vp_(iblock, k, j, i) * cag_);
    }
    v_dlv_(iblock, k, j, i) += v_odz_pt_(k, 0) * (diff_v1 - diff_v2);
    return ;
  }
 private:
  const double cag_ = CppPconstMod::cag;
  const double sag_ = CppPconstMod::sag;
  const ViewDouble2D v_odz_pt_ = *p_v_odz_pt;
  const ViewDouble3D v_sv_     = *p_v_sv;
  const ViewDouble3D v_snlat_  = *p_v_snlat;
  const ViewDouble4D v_up_     = *p_v_up;
  const ViewDouble4D v_vp_     = *p_v_vp;
  const ViewDouble4D v_dlv_    = *p_v_dlv;
  const ViewDouble4D v_viv_    = *p_v_viv;
  const ViewDouble4D v_wka_    = *p_v_wka;
  const ViewDouble4D v_akmu_   = *p_v_akmu;
};
//===========================

KOKKOS_REGISTER_FOR_2D(FunctorReadyc1,  FunctorReadyc1)
// KOKKOS_REGISTER_FOR_3D(FunctorReadyc2,  FunctorReadyc2)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc3,  FunctorReadyc3)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc4,  FunctorReadyc4)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc5,  FunctorReadyc5)
// KOKKOS_REGISTER_FOR_3D(FunctorReadyc51,  FunctorReadyc51)
// KOKKOS_REGISTER_FOR_2D(FunctorReadyc52,  FunctorReadyc52)
// KOKKOS_REGISTER_FOR_3D(FunctorReadyc53,  FunctorReadyc53)
// KOKKOS_REGISTER_FOR_2D(FunctorReadyc54,  FunctorReadyc54)
// KOKKOS_REGISTER_FOR_3D(FunctorReadyc55,  FunctorReadyc55)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc6,  FunctorReadyc6)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc7,  FunctorReadyc7)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc8,  FunctorReadyc8)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc9,  FunctorReadyc9)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc10, FunctorReadyc10)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc11, FunctorReadyc11)
KOKKOS_REGISTER_FOR_3D(FuncAdvMomCen1,  FuncAdvMomCen1)
KOKKOS_REGISTER_FOR_3D(FuncAdvMomFlu1,  FuncAdvMomFlu1)
KOKKOS_REGISTER_FOR_3D(FuncAdvMomCen2, FuncAdvMomCen2)
KOKKOS_REGISTER_FOR_3D(FuncAdvMomFlu2,  FuncAdvMomFlu2)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc14, FunctorReadyc14)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc15, FunctorReadyc15)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc16, FunctorReadyc16)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc17, FunctorReadyc17)
KOKKOS_REGISTER_FOR_3D(FuncReadyc15Del2, FuncReadyc15Del2)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc18, FunctorReadyc18)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc19, FunctorReadyc19)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc20, FunctorReadyc20)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc21, FunctorReadyc21)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc22, FunctorReadyc22)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYC_HPP_
