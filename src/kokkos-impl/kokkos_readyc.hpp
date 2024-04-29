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
using KokkosGrid::      p_v_fcort;
using KokkosGrid::      p_v_uarea;
using KokkosGrid::      p_v_tarea_r;
using KokkosGrid::      p_v_uarea_r;
using KokkosPconstMod:: p_v_akmu;
using KokkosPconstMod:: p_v_akt;
using KokkosPconstMod:: p_v_akmt;
using KokkosPconstMod:: p_v_ak_tide;
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

using KokkosTmpVar::p_v_curl;
using KokkosTmpVar::p_v_wp12;
using KokkosTmpVar::p_v_wp13;
using KokkosTmpVar::p_v_d2uk;
using KokkosTmpVar::p_v_d2vk;
using KokkosTmpVar::p_v_uv_ws_face;

using KokkosTmpVar::p_v_wk1;
using KokkosTmpVar::p_v_wk2;
using KokkosTmpVar::p_v_wk3;
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
class FunctorReadyc2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    const double epsln = 1.e-25;
    if (i < (NX_BLOCK-1) && j >= 1) {
      v_wka(iblock, k, j, i) = v_vit_(iblock, k, j, i) / (
          v_viv_(iblock, k, j  , i  ) * v_at0_ (iblock, j, i) +
          v_viv_(iblock, k, j-1, i  ) * v_atn_ (iblock, j, i) +
          v_viv_(iblock, k, j  , i+1) * v_ate_ (iblock, j, i) +
          v_viv_(iblock, k, j-1, i+1) * v_atne_(iblock, j, i)
              + epsln);
    }
    return ;
  }
 private:
  const ViewDouble3D v_at0_  = *p_v_at0;
  const ViewDouble3D v_atn_  = *p_v_atn;
  const ViewDouble3D v_ate_  = *p_v_ate;
  const ViewDouble3D v_atne_ = *p_v_atne;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka   = *p_v_wka;
};

class FunctorReadyc3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    ugrid_to_tgrid (iblock, k, j, i, v_wp12_, v_up_);
    ugrid_to_tgrid (iblock, k, j, i, v_wp13_, v_vp_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void ugrid_to_tgrid (const int &iblock, 
			const int &k, const int &j, const int &i,
          const ViewDouble4D &v_tgrid,
              const ViewDouble4D &v_ugrid) const {
    if (i < (NX_BLOCK-1) && j >= 1) {
      v_tgrid(iblock, k, j, i) =  v_wka_ (iblock, k, j  , i  ) * 
         (v_at0_ (iblock, j, i) * v_ugrid(iblock, k, j  , i  ) + 
          v_atn_ (iblock, j, i) * v_ugrid(iblock, k, j-1, i  ) +
          v_ate_ (iblock, j, i) * v_ugrid(iblock, k, j  , i+1) +
          v_atne_(iblock, j, i) * v_ugrid(iblock, k, j-1, i+1));
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
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_wp12_ = *p_v_wp12;
  const ViewDouble4D v_wp13_ = *p_v_wp13;
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
#ifdef CANUTO
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

#ifdef CANUTO
class FunctorReadyc51 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;
      v_wp3_(k, j, i) = 0.0;
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;
      if (k < kmt_m1 && v_vit_(iblock, k+1, j, i) > 0.0) {
        v_wp3_(k, j, i) = (v_pdensity_(iblock, k, j, i) + (v_pdensity_(iblock, k+1, j, i) 
            -  v_pdensity_(iblock, k, j, i)) * v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1))) * 1.e-3;
      }
    }
    return ;
  };
 private:
  const ViewInt3D    v_kmt_      = *p_v_kmt;
  const ViewDouble1D v_dzp_      = *p_v_dzp;
  const ViewDouble3D v_wp3_      = *p_v_wp3;
  const ViewDouble4D v_vit_      = *p_v_vit;
  const ViewDouble4D v_pdensity_ = *p_v_pdensity;
};

class FunctorReadyc52 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;
      double wp1[KM], wp4[KM];
      double wp5[KM], wp6[KM], wp7[KM], wp8[KM];
      for (int k = 0; k < KM; ++k) {
        wp1[k] = 0.0;
        wp4[k] = 0.0;
        wp5[k] = 0.0;
        wp6[k] = 0.0;
        wp7[k] = 0.0;
        wp8[k] = 0.0;
      }
      for (int k = 0; k < kmt_m1; ++k) {
        wp8[k] = - v_vit_(iblock, k+1, j, i) * v_zkp_(k+1) * 1.0e+2;
      }
      for (int k = 0; k < kmt_m1; ++k) {
        if (v_vit_(iblock, k+1, j, i) > 0.0) {
          wp1[k] = v_at_(iblock, 0, k, j, i) - (v_at_(iblock, 0, k, j, i) - v_at_(iblock, 0, k+1, j, i)) 
              * v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1));
          wp4[k] = v_rit_(iblock, k, j, i);
          wp5[k] = v_ricdttms_(iblock, k, j, i);
          wp6[k] = v_s2t_(iblock, k, j, i);
          wp7[k] = v_rict_(iblock, k, j, i);
        }
      }

      double wp10 = v_vit_(iblock, 0, j, i) * v_buoytur_(iblock, j, i) * 1.0e+4;
      double wp11 = v_vit_(iblock, 0, j, i) * v_buoysol_(iblock, j, i) * 1.0e+4;

#ifdef CANUTO2010
      // Fortran
      call canuto_2010_interface(wk1, wk2, wk3, wk4, 
          amld(i,j,iblock) ,tau_mag, wp1, wp2, wp3, wp4, wp5, 
          wp7, wp6, ulat(i,j,iblock)/degtorad , wp8, 
          kmt(i,j,iblock), iblock, j, i, isc)
#endif // CANUTO2010
      double mldtmp;

      double akm_back[KM-1];
      double akt_back[KM-1];
      double aks_back[KM-1];

      turb_2(wp8, wp1, v_wp3_, wp4, wp5, wp6, 
          dfricmx_ * 1.0e+4, dwndmix_ * 1.0e+4, 
              akm_back, akt_back, aks_back,
                  wp7, wp10, wp11, v_fcort_(iblock, j, i), 
                      mldtmp, v_wk1_, v_wk2_, v_wk3_, 
                          kmt_m1, KM-1, 1, 0, j, i);

      v_amld_(iblock, j, i) = mldtmp * 1.0e-2;
    }
    return ;
  };
 private:
  const double dfricmx_ = CppPconstMod::dfricmx;
  const double dwndmix_ = CppPconstMod::dwndmix;

  const double b1_         = CppCanutoMod::b1;
  const double dri_        = CppCanutoMod::dri;
  const double rri_        = CppCanutoMod::rri;
  const double rnd2on2_    = CppCanutoMod::rnd2on2;
  const double dand2on2_   = CppCanutoMod::dand2on2;
  const double deltheta_r_ = CppCanutoMod::deltheta_r;
  const double theta_rcrp_ = CppCanutoMod::theta_rcrp;
  const double theta_rcrn_ = CppCanutoMod::theta_rcrn;

  const ViewDouble1D v_zkp_        = *p_v_zkp;
  const ViewDouble1D v_dzp_        = *p_v_dzp;
  const ViewDouble3D v_wk1_        = *p_v_wk1;
  const ViewDouble3D v_wk2_        = *p_v_wk2;
  const ViewDouble3D v_wk3_        = *p_v_wk3;
  const ViewDouble3D v_wp3_        = *p_v_wp3;

  const ViewInt1D    v_irimax_     = *p_v_irimax;
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble1D v_rib_        = *p_v_rib;
  const ViewDouble1D v_ridb_       = *p_v_ridb;
  const ViewDouble1D v_sma1_       = *p_v_sma1;
  const ViewDouble1D v_sha1_       = *p_v_sha1;
  const ViewDouble1D v_ssa1_       = *p_v_ssa1;
  const ViewDouble1D v_and2on2a1_  = *p_v_and2on2a1;
  const ViewDouble1D v_amtaun2a1_  = *p_v_amtaun2a1;
  const ViewDouble1D v_back_ra_r_  = *p_v_back_ra_r;
  const ViewDouble1D v_sm_r1_      = *p_v_sm_r1;
  const ViewDouble1D v_sh_r1_      = *p_v_sh_r1;
  const ViewDouble1D v_ss_r1_      = *p_v_ss_r1;
  const ViewDouble1D v_slq2_r1_    = *p_v_slq2_r1;
  const ViewDouble2D v_smb_        = *p_v_smb;
  const ViewDouble2D v_shb_        = *p_v_shb;
  const ViewDouble2D v_ssb_        = *p_v_ssb;
  const ViewDouble2D v_slq2b_      = *p_v_slq2b;
  const ViewDouble3D v_amld_       = *p_v_amld;
  const ViewDouble3D v_fcort_      = *p_v_fcort;
  const ViewDouble3D v_buoytur_    = *p_v_buoytur;
  const ViewDouble3D v_buoysol_    = *p_v_buoysol;
  const ViewDouble4D v_s2t_        = *p_v_s2t;
  const ViewDouble4D v_rit_        = *p_v_rit;
  const ViewDouble4D v_rict_       = *p_v_rict;
  const ViewDouble4D v_ricdttms_   = *p_v_ricdttms;
  const ViewDouble4D v_vit_        = *p_v_vit;
  const ViewDouble5D v_at_         = *p_v_at;

  KOKKOS_INLINE_FUNCTION void turb_2(
      const double (&z)[KM],        // wp8
      const double (&t)[KM],        // wp1
      const ViewDouble3D &v_rh,     // wp3
            double (&ri)[KM],       // wp4
      const double (&rid)[KM],      // wp5
            double (&s2)[KM],       // wp6
      const double &fricmx,         // DFRICMX * 1.0e+4
      const double &wndmix,         // DWNDMIX * 1.0e+4
            double (&v_back)[KM-1], // akm_back
            double (&t_back)[KM-1], // akt_back
            double (&s_back)[KM-1], // aks_back
      const double (&an2)[KM],      // wp7
      const double &buoytur,        // wp10
      const double &buoysol,        // wp11
      const double &coriol,         // fcort[iblock][j][i]
            double &amld,           // mldtmp
      const ViewDouble3D &v_akm,    // wk1
      const ViewDouble3D &v_akh,    // wk2
      const ViewDouble3D &v_aks,    // wk3
      const int    &n,              // iwk
      const int    &nmax,           // km - 1
      const int    &isurfuse,       // 1
      const int    &ifextermld,     // 0
      const int &j, const int &i    ) const {

    double an = 0.0;
    double sh_back = 0.0;
    double sm_back = 0.0;
    double ss_back = 0.0;
    double slq2_back = 0.0;
 
    double ri1  = 0.0;
    double rid1 = 0.0;
 
    int ifbelow;
    int ifpureshear = 0;
 
    int ifnofsmall = 0;
 
    int nb = std::max(0, n);
    double visc_cbu_limit1 = fricmx;
    double diff_cbt_limit1 = fricmx;
    double buoytot = buoytur + buoysol;
 
    // if (ifextermld == 0) {
      if (n > 0) {
          // if idefmld == 0
          // idefmld = 0 in canuto_mod_2002
          formld(z, t, amld, n);
      } else if (n == 0) {
        amld = z[0];
      } else {
        amld = 0.0;
      }
    // }
 
    double al0 = 0.17 * amld;
 
    ifbelow = 0;
 
    // int icall(0);
    // int ipoint(0);
    // int iproblem(0)
    // int inegproblem(0);
 
    double sm, sh, ss;
    double slq2, epson2;

    for (int k = 0; k < n; ++k) {
  
      double amtaun2(0.0);
  
      ifnofsmall = 0;
      if (an2[k] >= 0.0) {
        an = sqrt(an2[k]);
      }
      if ((an / fabs(coriol)) < 1.0) {
        ifnofsmall = 1;
      }
  
      ri1 = ri[k];
  
      rid1 = rid[k];
      const double and2 = rid[k] / (ri[k] + 1.0e-25) * an2[k];
      double and2on2 = and2 / (an2[k] + 1.0e-25);
  
      if (an2[k] < v_rib_(0) * s2[k]) {
        ifpureshear = 1;
      } else {
        ifpureshear = 0;
      }

      if (ifpureshear == 1) {
        const int imax = MT;
  
        interp1d_expabs(and2on2, amtaun2, 
            sm, sh, ss, imax, dand2on2_, rnd2on2_);
  
        slq2 = (-amtaun2) / ((b1_ * b1_) * ri[k] + 1.0e-25);
      }

      if (ifpureshear != 1) {
        interp2d_expabs(ri1, rid1, slq2, sm, sh, ss);
      }
  
      if (ifpureshear != 1) {
        if (slq2 < 0.0) {
          return ;
        }
        if (slq2 == 0) {
          ifbelow = 1;
        }
      }
  
      // bool lifupper = (((IFEPSON2  < 2) || (ifbelow == 0) || 
      //     ((IFDEEPLAT > 0) && (ifnofsmall == 1)) || 
      //     ((ri1 < 0.0) && (k <= 2))) && (slq2 > 0.0));
      bool lifupper = (((ifbelow == 0) || (ifnofsmall == 1) || 
          ((ri1 < 0.0) && (k <= 2))) && (slq2 > 0.0));
  
      // ++ipoint;
      // if (k == 0) {
      //   ++icall;
      // }
      double epsy = 0.0;
      bool lifepsy = false;
  
      // if ((isurfuse == 1) && n > 0) {
      if (n > 0) {
        if (slq2 == 0.0) {
          epsy = 0.0;
        } else {
          epsy = - buoytot
              / ((1.0 / ((b1_ * b1_) * slq2)) - 0.5 * sm);
        }
        lifepsy = ((epsy >= 0.0) && lifupper);
        // if ((epsy[k] < 0.0) && lifupper) {
        //   // ++iproblem;
        //   // if (ri1 < 0.0) {
        //   //   ++inegproblem;
        //   // }
        // }
      }
  
      double akz = 0.4 * z[k];
      double al = akz * al0 / (al0 + akz);
  
      double al2 = al * al;
  
      // if (!(((IFEPSON2 == 2) && (ifbelow == 1)) || lifepsy[k])) {
      if (!(ifbelow == 1 || lifepsy)) {
        if (ri1 > 0.0) {
          const double anlq2 = slq2 * ri1;
          if (anlq2 > 0.281) {
            al2 = 0.281 / anlq2 * al2;
            slq2 = 0.281 / (ri1 + 1.0e-20);
          }
        }
      }
  
      double epson2_;
      if (an2[k] < 0.0) {
        epson2_ = EPSON2__;
      } else {
        double eplatidepend;
        if (ifnofsmall == 1) {
          eplatidepend = 0.0;
        } else {
	 	      eplatidepend = eplatidepend_(fabs(coriol), an);
        }
  
        eplatidepend = fmax(eplatidepend, EPLATIDEPENDMIN);
  
        epson2_ = EPSON2__ * eplatidepend;
      }
  
      epson2 = epson2_;
  
      int ifrafglt = 0;
  
      double rit(0.0), ric(0.0);
      double back_ra_r1;
      double back_ri1,  back_ric1;
      double back_rid1, back_rit1;
      double theta_r(0.0);
      double deltheta_r1;
      int itheta_r0(0), itheta_r1(0);
      int jtheta_r0(0);
      double ra_r;
  
      const double tmp_deltheta_r = 1.0 / deltheta_r_;
      if (ri[k] <= 0.0) {
        back_ra_r1 = 0.0;
        back_rit1  = 0.0;
        back_ric1  = 0.0;
        back_ri1   = 0.0;
        back_rid1  = 0.0;
      } else {
        // rit = (ri[k] + rid[k]) / 2.0;
        // ric = (ri[k] - rid[k]) / 2.0;
        rit = (ri[k] + rid[k]) * 0.5;
        ric = (ri[k] - rid[k]) * 0.5;
        ra_r = sqrt((rit * rit) + (ric * ric));
  
        if (rit == 0.0) {
          if (ric == 0.0) {
            theta_r = atan(1.0);
          } else {
            // theta_r = PI / 2.0;
            theta_r = PI * 0.5;
          }
        } else {
          theta_r = atan(ric / rit);
        }
        // if (fabs(theta_r) > (PI / 2.0)) {
        if (fabs(theta_r) > (PI * 0.5)) {
          return ;
        }
        // if (theta_r < (-PI) / 4.0) {
        if (theta_r < (-PI) * 0.25) {
          theta_r += PI;
        }
  
        jtheta_r0 = static_cast<int>(
            // (theta_r + (PI / 4.0)) / deltheta_r_);
            (theta_r + (PI * 0.25)) * tmp_deltheta_r);
  
        itheta_r0 = jtheta_r0 - N_THETA_R_OCT;
        itheta_r1 = itheta_r0 + 1;
  
        double theta_r0 = itheta_r0 * deltheta_r_;
        double theta_r1 = itheta_r1 * deltheta_r_;
  
        if ((theta_r0 <= theta_rcrp_) &&
             (theta_r > theta_rcrp_)) {
          theta_r    = theta_r1;
          theta_r0   = theta_r1;
          itheta_r0  = itheta_r1;
          itheta_r1 += 1;
          theta_r1  += deltheta_r_;
        } else if ((theta_r1 >= theta_rcrn_) &&
                   (theta_r  < theta_rcrn_)) {
          theta_r    = theta_r0;
          theta_r1   = theta_r0;
          itheta_r1  = itheta_r0;
          itheta_r0 -= 1;
          theta_r0  -= deltheta_r_;
        }
  
        if ((itheta_r1 > 3 * N_THETA_R_OCT) ||
            (itheta_r0 < - N_THETA_R_OCT)) {
          return ;
        }
  
        deltheta_r1 = theta_r - theta_r0;
  
        const double delback_ra_r = 
            v_back_ra_r_(N_THETA_R_OCT + itheta_r1) 
                - v_back_ra_r_(N_THETA_R_OCT + itheta_r0);
  
        const double dback_ra_r_o_dtheta = 
            // delback_ra_r / deltheta_r_;
            delback_ra_r * tmp_deltheta_r;
  
        back_ra_r1 = v_back_ra_r_(N_THETA_R_OCT + itheta_r0) 
            + deltheta_r1 * dback_ra_r_o_dtheta;
  
        ifrafglt = 0;
  
        if ((theta_r <= theta_rcrp_) ||
            (theta_r >= theta_rcrn_)) {
          if (back_ra_r1 > ra_r) {
            ifrafglt = 1;
            back_ra_r1 = ra_r;
          }
        }
  
        if (back_ra_r1 < 0.0) {
          return ;
        }
        back_rit1 = cos(theta_r) * back_ra_r1;
        back_ric1 = sin(theta_r) * back_ra_r1;
        back_ri1  = back_rit1 + back_ric1;
        back_rid1 = back_rit1 - back_ric1;
      }
  
      if ((IFSALBACK != 4) && (ri[k] > 0.0)) {
  
        // if ((IFBG_THETA_INTERP == 0) ||
        //     (ifrafglt == 1)) {
        if (ifrafglt == 1) {
          interp2d_expabs(back_ri1, back_rid1,
              slq2_back, sm_back, sh_back, ss_back);
        // } else if(IFBG_THETA_INTERP == 1) {
        } else if(true) {
          deltheta_r1 = theta_r - itheta_r0 * deltheta_r_;
  
          const double delsm_back = v_sm_r1_(N_THETA_R_OCT + itheta_r1)
              - v_sm_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dsm_back_o_dtheta = delsm_back / deltheta_r_;
          const double dsm_back_o_dtheta = delsm_back * tmp_deltheta_r;
  
          sm_back = v_sm_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dsm_back_o_dtheta;
  
          const double delsh_back = v_sh_r1_(N_THETA_R_OCT + itheta_r1)
              - v_sh_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dsh_back_o_dtheta = delsh_back / deltheta_r_;
          const double dsh_back_o_dtheta = delsh_back * tmp_deltheta_r;
  
          sh_back = v_sh_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dsh_back_o_dtheta;
  
          const double delss_back = v_ss_r1_(N_THETA_R_OCT + itheta_r1)
              - v_ss_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dss_back_o_dtheta = delss_back / deltheta_r_;
          const double dss_back_o_dtheta = delss_back * tmp_deltheta_r;
  
          ss_back = v_ss_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dss_back_o_dtheta;
  
          const double delslq2_back = v_slq2_r1_(N_THETA_R_OCT + itheta_r1)
              - v_slq2_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dslq2_back_o_dtheta = delslq2_back / deltheta_r_;
          const double dslq2_back_o_dtheta = delslq2_back * tmp_deltheta_r;
  
          slq2_back = v_slq2_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dslq2_back_o_dtheta;
        } else {
          return ;
        }
        if (slq2_back < 0.0) {
          return ;
        }
      } // not go to 19

      // 19 continue
      if (ri1 < 0.0) {
        sm_back = 0.0;
        sh_back = 0.0;
        ss_back = 0.0;
      }
  
      if ((sm_back < 0.0) || (sh_back < 0.0) ||
          (ss_back < 0.0)) {
        return ;
      }
      
      //double tmp_back;
  
      const double tmp_back = 0.5 * b1_ * b1_
          * back_ri1 * slq2_back * epson2;
  
      v_back[k] = tmp_back * sm_back;
      t_back[k] = tmp_back * sh_back;
      s_back[k] = tmp_back * ss_back;
  
      if ((v_back[k] < 0.0) || (t_back[k] < 0.0)
          || (s_back[k] < 0.0)) {
        return ;
      }
      if((ri[k] > 0.0) && ((v_back[k] == 0.0)
          || (t_back[k] == 0.0) || (s_back[k] == 0.0))) {
        return ;
      }
  
      double tmp(0.0);
  
      // if ((IFEPSON2 == 2) && (ifbelow == 1)) {
      //   if ((ri1 >= 0.0) || ((IFDEEPLAT == 2) && (ifnofsmall == 1))) {
      if (ifbelow == 1) {
        if (ri1 >= 0.0) {
          tmp = 0.5 * b1_ * b1_ * ri1 * slq2 * epson2;
        } else if (k > 1) {
          double delz, delrh, del2rh;
          if (k == n - 1) {
            delz = z[k] - z[k-1];
            delrh = v_rh(k, j, i) - v_rh(k-1, j, i);
            del2rh = v_rh(k, j, i) - 2.0 * v_rh(k-1, j, i) + v_rh(k-2, j, i);
          } else {
            delz = z[k+1] - z[k-1];
            delrh = v_rh(k+1, j, i) - v_rh(k-1, j, i);
            del2rh = v_rh(k+1, j, i) - 2.0 * v_rh(k, j, i) + v_rh(k-1, j, i);
          }
  
          const double dzrh = delrh / delz;
          const double d2zrh = 4.0 * del2rh / (delz * delz);
          const double rdzlndzrh = dzrh / d2zrh;
          const double al0deep = 0.17 * fabs(rdzlndzrh);
          akz = 0.4 * z[k];
          const double aldeep = akz * al0deep / (al0deep + akz);
          al2 = aldeep * aldeep;
  
          if (ifpureshear == 1) {
            // GO TO 21
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
       
            v_akm(k, j, i) = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            v_akh(k, j, i) = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            v_aks(k, j, i) = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
          } else if (IFSHEARMIN) {
            // s2[k] = fmax(s2[k], S2MIN);
            s2[k] = fmax(s2[k], S2MIN);
          }
          tmp = 0.5 * b1_ * al2
              * sqrt(s2[k] / (slq2 + 1.0e-40));
        } else {
          if (ifpureshear == 1) {
            // GO TO 21
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
       
            v_akm(k, j, i) = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            v_akh(k, j, i) = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            v_aks(k, j, i) = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
          } else if (IFSHEARMIN) {
            s2[k] = fmax(s2[k], S2MIN);
          }
          if (lifepsy) {
            tmp = 0.5 * epsy / (s2[k] + 1.0e-40);
          } else {
            tmp = 0.5 * b1_ * al2
                * sqrt(s2[k] / (slq2 + 1.0e-40));
          }
        }
      } else {
        if (ifpureshear == 1) {
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
            v_akm(k, j, i) = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            v_akh(k, j, i) = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            v_aks(k, j, i) = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
        } else if (IFSHEARMIN) {
          s2[k] = fmax(s2[k], S2MIN);
        }
        if (lifepsy) {
          tmp = 0.5 * epsy / (s2[k] + 1.0e-40);
        } else {
          tmp = 0.5 * b1_ * al2
              * sqrt(s2[k] / (slq2 + 1.0e-40));
        }
      }
  
      // 21
      if (ifpureshear == 1) {
        tmp = 0.5 * (b1_ * b1_) * al2
            * sqrt(-an2[k] / amtaun2);
      }
  
      v_akm(k, j, i) = fmin(tmp * sm + v_back[k],
          visc_cbu_limit1);
      v_akh(k, j, i) = fmin(tmp * sh + t_back[k],
          diff_cbt_limit1);
      v_aks(k, j, i) = fmin(tmp * ss + s_back[k],
          diff_cbt_limit1);
    } // End for k in 0:n

    for (int k = nb+1; k < nmax; ++k) {
      v_akm(k, j, i) = 0.0;
      v_akh(k, j, i) = 0.0;
      v_aks(k, j, i) = 0.0;
    }
 
    if (n > 0) {
      if (v_akm(0, j, i) < wndmix) {
        v_akm(0, j, i) = wndmix;
      }
      if (v_akh(0, j, i) < wndmix) {
        v_akh(0, j, i) = wndmix;
      }
      if (v_aks(0, j, i) < wndmix) {
        v_aks(0, j, i) = wndmix;
      }
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void formld(
      const double (&z)[KM], const double (&t)[KM],
          double &amld, const int &n) const {
    for (int k = 0; k < n; ++k) {
      if (fabs(t[k] - t[0]) > 0.1) {
#ifdef D_PRECISION
        const double tm = t[0] - sign_double(0.1, t[0] - t[k]);
#else
        const double tm = t[0] - sign_float(0.1,   t[0] - t[k]);
#endif // D_PRECISION
        amld = z[k] + (z[k-1] - z[k]) *
            (tm - t[k]) / (t[k-1] - t[k] + 1.e-20);
        return ;
      }
    }
    amld = z[n-1];
    return ;
  }

  KOKKOS_INLINE_FUNCTION int sign_int(const double &x, const double &y) const {
    return y >= 0 ? std::abs(static_cast<int>(x)) 
        : - std::abs(static_cast<int>(x));
  }

  KOKKOS_INLINE_FUNCTION int sign_float(const double &x, const double &y) const {
    return y >= 0 ? std::abs(static_cast<float>(x)) 
       : - std::abs(static_cast<float>(x));
  }

  KOKKOS_INLINE_FUNCTION double sign_double(const double &x, const double &y) const {
    return y >= 0.e0 ? fabs(x) : - fabs(x);
  }

  KOKKOS_INLINE_FUNCTION void interp1d_expabs(
      double &x,            // and2on2
      double &slq2,         // amtaun2
      double &sm,           // sm
      double &sh,           // sh
      double &ss,           // ss
      const int &ixmax,     // imax
      const double &delta,  // dand2on2
      const double &rat     /* rnd2on2 */ ) const {
    // printf ("interp1d_expabs\n");
  
    int lx0(0), lx1(0);
    
    if (x > v_and2on2a1_(MT + ixmax)) {
      x = v_and2on2a1_(MT + MT);
    } else if (x < v_and2on2a1_(0)) {
      x = v_and2on2a1_(0);
    }
  
    if (fabs(x) < v_and2on2a1_(MT + MT0)) {
#ifdef D_PRECISION
      lx1 = static_cast<int>(x / delta) + round(sign_double(static_cast<double>(1.0), x)); 
#else  // D_PRECISION
      lx1 = static_cast<int>(x / delta) + round(sign_float(static_cast<float>(1.0), x)); 
#endif // D_PRECISION
    } else if (fabs(x) >= v_and2on2a1_(MT + MT)) {
#ifdef D_PRECISION
      lx0 = round(sign_double(static_cast<double>(MT), x));
#else  // D_PRECISION
      lx0 = round(sign_float(static_cast<float>(MT), x));
#endif // D_PRECISION
      lx1 = lx0;
    } else {
#ifdef D_PRECISION
      const double tabindx = sign_double(static_cast<double>(MT0)
         + ((log(fabs(x)) - log(v_and2on2a1_(MT + MT0))) / log(rat)), x);
#else  // D_PRECISION
      const float tabindx = sign_double(static_cast<float>(MT0)
         + ((log(fabs(x)) - log(v_and2on2a1_(MT + MT0))) / log(rat)), x);
#endif // D_PRECISION

#ifdef D_PRECISION
      lx1 = static_cast<int>(tabindx) + round(sign_double(
          static_cast<double>(1.0), x));
#else  // D_PRECISION
      lx1 = static_cast<int>(tabindx) + round(sign_float(
          static_cast<float>(1.0), x));
#endif // D_PRECISION
    }

    if (!(fabs(x) >= v_and2on2a1_(MT + MT))) {
      if (fabs(v_and2on2a1_(MT + lx1)) < fabs(x)) {
        lx1 += sign_int(1, lx1);
      } else if (fabs(v_and2on2a1_(MT + lx1 - sign_int(1, lx1))) > fabs(x)) {
        lx1 -= sign_int(1, lx1);
      }

#ifdef D_PRECISION
      lx0 = lx1 - round(sign_double(
          static_cast<double>(1.0), x));
#else  // D_PRECISION
      lx0 = lx1 - round(sign_float(
          static_cast<float>(1.0), x));
#endif // D_PRECISION
      if (x == 0.0) {
        lx1 = 1;
      }
    }
    // 252
    if ((x > 0.0 && (x < v_and2on2a1_(MT + lx0) || x > v_and2on2a1_(MT + lx1))) ||
        (x < 0.0 && (x > v_and2on2a1_(MT + lx0) || x < v_and2on2a1_(MT + lx1)))) {
      return;
    }

    const double deltaxta = 1.0 / (v_and2on2a1_(MT + lx1) - v_and2on2a1_(MT + lx0));

    const double deltax = x - v_and2on2a1_(MT + lx0);

    double dslq2_x;
    if (lx1 == lx0) {
      dslq2_x = 0.0;
    } else {
      dslq2_x = (v_amtaun2a1_(MT + lx1) - v_amtaun2a1_(MT + lx0)) * deltaxta;
    }

    slq2 = v_amtaun2a1_(MT + lx0) + dslq2_x * deltax;

    double dsm_x;
    if (lx1 == lx0) {
      dsm_x = 0.0;
    } else {
      dsm_x = (v_sma1_(MT + lx1) - v_sma1_(MT + lx0)) * deltaxta;
    }
    sm = v_sma1_(MT + lx0) + dsm_x * deltax;

    double dsh_x;
    if (lx1 == lx0) {
      dsh_x = 0.0;
    } else {
      dsh_x = (v_sha1_(MT + lx1) - v_sha1_(MT + lx0)) * deltaxta;
    }
    sh = v_sha1_(MT + lx0) + dsh_x * deltax;

    double dss_x;
    if (lx1 == lx0) {
      dss_x = 0.0;
    } else {
      dss_x = (v_ssa1_(MT + lx1) - v_ssa1_(MT + lx0)) * deltaxta;
    }
    ss = v_ssa1_(MT + lx0) + dss_x * deltax;
  
    return ;  
  }

  KOKKOS_INLINE_FUNCTION void interp2d_expabs(double &ri, double &rid,
      double &slq2, double &sm, double &sh, double &ss) const {

    const double tmp_ri  = 1.0 /ri;
    const double tmp_rid = 1.0 /rid;

    if (ri > v_rib_(MT + MT)) {
      if (fabs(rid) <= ri) {
        // rid = v_rib_(MT + MT) * (rid / ri);
        rid = v_rib_(MT + MT) * (rid * tmp_ri);
        ri  = v_rib_(MT + MT);
      } else if (rid > ri) {
        // ri  = v_ridb_(MT + MT) * (ri / rid);
        ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
        rid = v_ridb_(MT + MT);
      } else if (rid < -ri) {
        // ri  = v_ridb_(0) * (ri / rid);
        ri  = v_ridb_(0) * (ri * tmp_rid);
        rid = v_ridb_(0);
      }
    } else if (ri < v_rib_(0)) {
      if (fabs(rid) <= -ri) {
        // rid = v_rib_(0) * (rid / ri);
        rid = v_rib_(0) * (rid * tmp_ri);
        ri  = v_rib_(0);
      } else if (rid > -ri) {
        // ri  = v_ridb_(MT + MT) * (ri / rid);
        ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
        rid = v_ridb_(MT + MT);
      } else if (rid < ri) {
        // ri  = v_ridb_(0) * (ri / rid);
        ri  = v_ridb_(0) * (ri * tmp_rid);
        rid = v_ridb_(0);
      }
    } else if (rid > v_ridb_(MT + MT)) {
      // ri  = v_ridb_(MT + MT) * (ri / rid);
      ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
      rid = v_ridb_(MT + MT);
    } else if (rid < v_ridb_(0)) {
      // ri  = v_ridb_(0) * (ri / rid);
      ri  = v_ridb_(0) * (ri * tmp_rid);
      rid = v_ridb_(0);
    }

    int lrid0(0), lrid1(0);

    if (fabs(rid) < v_ridb_(MT + MT0)) {
#ifdef D_PRECISION
      lrid1 = static_cast<int>(rid / dri_)
          + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid1 = static_cast<int>(rid / dri_)
          + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION
    } else if (fabs(rid) >= v_ridb_(MT + MT)) {
#ifdef D_PRECISION
      lrid0 = round(sign_double(static_cast<double>(MT), rid));
#else // D_PRECISION
      lrid0 = round(sign_float(static_cast<float>(MT), rid));
#endif // D_PRECISION
      lrid1 = lrid0;
    } else {
#ifdef D_PRECISION
      const double tabindrid = sign_double(static_cast<double>(MT0)
          + ((log(fabs(rid)) - log(v_ridb_(MT + MT0)))
              / log(rri_)), rid);
#else // D_PRECISION
      const double tabindrid = sign_float(static_cast<float>(MT0)
          + ((log(fabs(rid)) - log(v_ridb_(MT + MT0)))
              / log(rri)), rid);
#endif // D_PRECISION

#ifdef D_PRECISION
      lrid1 = static_cast<int>(tabindrid)
          + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid1 = static_cast<int>(tabindrid)
          + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION

    }

    if (!(fabs(rid) >= v_ridb_(MT + MT))) {

      if (fabs(v_ridb_(MT + lrid1)) < fabs(rid)) {
        lrid1 += sign_int(1, lrid1);
      } else if (fabs(v_ridb_(MT + lrid1 - sign_int(1, lrid1))) > 
          fabs(rid)) {
        lrid1 -= sign_int(1, lrid1);
      }

#ifdef D_PRECISION
      lrid0 = lrid1 - round(sign_double(
          static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid0 = lrid1 - round(sign_float(
          static_cast<float>(1.0), rid));
#endif // D_PRECISION
      if (rid == 0.0) {
        lrid1 = 1;
      }
    }

    // 252
    if ((rid > 0.0 && (rid < v_ridb_(MT + lrid0) || rid > v_ridb_(MT + lrid1))) ||
        (rid < 0.0 && (rid > v_ridb_(MT + lrid0) || rid < v_ridb_(MT + lrid1)))) {
      return ;
    }
 
    if (ri > fmin(v_rib_(MT + v_irimax_(MT + lrid0)),
                  v_rib_(MT + v_irimax_(MT + lrid1)))) {
      slq2 = 0.0;
      sm   = 0.0;
      sh   = 0.0;
      ss   = 0.0;
      return ;
    }

    int lri0(0), lri1(0);

    if (fabs(ri) < v_rib_(MT + MT0)) {
#ifdef D_PRECISION
      lri1 = static_cast<int>(ri / dri_)
          + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri1 = static_cast<int>(ri / dri_)
          + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
    } else if (fabs(ri) >= v_rib_(MT + MT)) {
#ifdef D_PRECISION
      lri0 = round(sign_double(static_cast<double>(MT), ri));
#else // D_PRECISION
      lri0 = round(sign_float(static_cast<float>(MT), ri));
#endif // D_PRECISION
      lri1 = lri0;
    } else {
#ifdef D_PRECISION
      const double tabindri = sign_double(static_cast<double>(MT0)
          + ((log(fabs(ri)) - log(v_rib_(MT + MT0)))
              / log(rri_)), ri);
#else // D_PRECISION
      const double tabindri = sign_float(static_cast<float>(MT0)
          + ((log(fabs(ri)) - log(v_rib_(MT + MT0)))
              / log(rri)), ri);
#endif // D_PRECISION

#ifdef D_PRECISION
      lri1 = static_cast<int>(tabindri)
          + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri1 = static_cast<int>(tabindri)
          + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
    }
 
    if (!(fabs(ri) >= v_rib_(MT + MT))) {
      if (fabs(v_rib_(MT + lri1)) < fabs(ri)) {
        lri1 += sign_int(1, lri1);
      } else if (fabs(v_rib_(MT + lri1 - sign_int(1, lri1))) >
          fabs(ri)) {
        lri1 -= sign_int(1, lri1);
      }

#ifdef D_PRECISION
      lri0 = lri1 - round(sign_double(
          static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri0 = lri1 - round(sign_float(
          static_cast<float>(1.0), ri));
#endif // D_PRECISION

      if (ri == 0.0) {
        lri1 = 1;
      }
    }

    if ((ri > 0.0 && (ri < v_rib_(MT + lri0) || ri > v_rib_(MT + lri1))) ||
        (ri < 0.0 && (ri > v_rib_(MT + lri0) || ri < v_rib_(MT + lri1)))) {
      return ;
    }

    const double deltaridta = 1.0 / (v_ridb_(MT + lrid1) - v_ridb_(MT + lrid0));
    const double deltarita  = 1.0 / (v_rib_(MT + lri1) - v_rib_(MT + lri0));
    const double deltarid = rid - v_ridb_(MT + lrid0);
    const double deltari  = ri  - v_rib_(MT + lri0);

    double dslq2_rid;

    if (lrid1 == lrid0) {
      dslq2_rid = 0.0;
    } else {
      dslq2_rid = (v_slq2b_(MT+lrid1, MT+lri0) - v_slq2b_(MT+lrid0, MT+lri0))
          * deltaridta;
    }

    double dslq2_ri;
    if (lri1 == lri0) {
      dslq2_ri = 0.0;
    } else {
      dslq2_ri = (v_slq2b_(MT+lrid0, MT+lri1) - v_slq2b_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    slq2 = v_slq2b_(MT+lrid0, MT+lri0)
        + dslq2_ri * deltari + dslq2_rid * deltarid;

    double dsm_rid;
    if (lrid1 == lrid0) {
      dsm_rid = 0.0;
    } else {
      dsm_rid = (v_smb_(MT+lrid1, MT+lri0) - v_smb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }

    double dsm_ri;

    if (lri1 == lri0) {
      dsm_ri = 0.0;
    } else {
      dsm_ri = (v_smb_(MT+lrid0, MT+lri1) - v_smb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    sm = v_smb_(MT+lrid0, MT+lri0) 
        + dsm_ri * deltari + dsm_rid * deltarid;

    double dsh_rid;
    if (lrid1 == lrid0) {
      dsh_rid = 0.0;
    } else {
      dsh_rid = (v_shb_(MT+lrid1, MT+lri0) - v_shb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }
    double dsh_ri;
    if (lri1 == lri0) {
      dsh_ri = 0.0;
    } else {
      dsh_ri = (v_shb_(MT+lrid0, MT+lri1) - v_shb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    sh = v_shb_(MT+lrid0, MT+lri0) 
        + dsh_ri * deltari + dsh_rid * deltarid;
    
    double dss_rid;
    if (lrid1 == lrid0) {
      dss_rid = 0.0;
    } else {
      dss_rid = (v_ssb_(MT+lrid1, MT+lri0) - v_ssb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }
    double dss_ri;
    if (lri1 == lri0) {
      dss_ri = 0.0;
    } else {
      dss_ri = (v_ssb_(MT+lrid0, MT+lri1) - v_ssb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    ss = v_ssb_(MT+lrid0, MT+lri0) 
        + dss_ri * deltari + dss_rid * deltarid;

    return ;
  }

  KOKKOS_INLINE_FUNCTION double acosh1(const double &x) const {
    return log(x + sqrt((x * x) - 1.0));
  }
  KOKKOS_INLINE_FUNCTION double wavelat(const double &xf, const double &yn) 
      const {
    return xf * acosh1(yn / xf);
  }

  KOKKOS_INLINE_FUNCTION double eplatidepend_(const double &f, const double &an) 
      const {
    double f_30;
    double anum, den;
    double pi1, omega;
    const double an0 = 5.24e-3;
    pi1   = 4.0 * atan(1.0);
    omega = pi1 / 43082.0e0;
    f_30  = omega;
    den   = wavelat(f_30, an0);
    anum  = wavelat(f, an);
    return anum / den;
  }
};
class FunctorReadyc53 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;

#ifdef TIDEMIX
      v_ak_tide_(iblock, k, j, i) = 0.0;
      if (k < kmt_m1) {
        v_ak_tide_(iblock, k, j, i) = BACK_TIDALMIXING + MIXING_EF * 
            LOCAL_MIXING_FRACTION * v_wave_dis_(iblock, j, i) * 
            v_fz_tide_(iblock, k, j, i) /
            (fmax(v_rict_(iblock, k, j, i), 1.0e-8) * v_wp3_(k, j, i) * 1000.0);

        v_ak_tide_(iblock, k, j, i) = fmin(v_ak_tide_(iblock, k, j, i),
            MAX_TIDALMIXING);

        v_richardson_(iblock, k, j, i) = v_rict_(iblock, k, j, i);
        v_fztidal_   (iblock, k, j, i) = v_fz_tide_(iblock, k, j, i);
        v_wp3_tidal_ (iblock, k, j, i) = v_wp3_(k, j, i);
      }
#endif // TIDEMIX
    }
    return ;
  };
 private:
  const ViewDouble3D v_wp3_        = *p_v_wp3;
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble3D v_wave_dis_   = *p_v_wave_dis;
  const ViewDouble4D v_vit_        = *p_v_vit;
  const ViewDouble4D v_rict_       = *p_v_rict;
  const ViewDouble4D v_fztidal_    = *p_v_fztidal;
  const ViewDouble4D v_ak_tide_    = *p_v_ak_tide;
  const ViewDouble4D v_fz_tide_    = *p_v_fz_tide;
  const ViewDouble4D v_wp3_tidal_  = *p_v_wp3_tidal;
  const ViewDouble4D v_richardson_ = *p_v_richardson;
};

class FunctorReadyc54 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
#ifdef TIDEMIX
#ifdef CANUTOMIXOUT
#endif // CANUTOMIXOUT
      const int iblock = 0;
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;
      double pre_ak_tide = v_ak_tide_(iblock, kmt_m1-1, j, i);
      for (int k = kmt_m1 - 2; k >= 0; --k) {
        v_ak_tide_(iblock, k, j, i) = fmin(
            v_ak_tide_(iblock, k, j, i), pre_ak_tide);
        pre_ak_tide = v_ak_tide_(iblock, k, j, i);
      }
#endif // TIDEMIX
    }
    return ;
  };
 private:
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble4D v_ak_tide_ = *p_v_ak_tide;
};

class FunctorReadyc55 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;
      const float tmp_ncc = 1.0 / static_cast<float>(ncc_);

#ifdef TIDEMIX
      v_akmt_(iblock, k, j, i)    =  v_wk1_(k, j, i) * 1.0e-4 +
          v_ak_tide_(iblock, k, j, i) * 5.0;
      v_akt_(iblock, 0, k, j, i) += (v_wk2_(k, j, i) * 1.0e-4 + 
          v_ak_tide_(iblock, k, j, i)) * tmp_ncc;
      v_akt_(iblock, 1, k, j, i) += (v_wk3_(k, j, i) * 1.0e-4 + 
          v_ak_tide_(iblock, k, j, i)) * tmp_ncc;
#ifdef BCKMEX
      v_akmt_(iblock, k, j, i)   += diff_back[iblock][j][i] * 
          10.0 * 1.0e-4;
      v_akt_(iblock, 0, k, j, i) += diff_back[iblock][j][i] /
          static_cast<float>(ncc) * 1.0e-4;
      v_akt_(iblock, 1, k, j, i) += diff_back[iblock][j][i] /
          static_cast<float>(ncc) * 1.0e-4;
#endif // BCKMEX
#else // TIDEMIX
      v_akmt_(iblock, k, j, i)    = v_wk1_(k, j, i) * 1.0e-4;
      v_akt_(iblock, 0, k, j, i) += v_wk2_(k, j, i) * 1.0e-4 / 
          static_cast<float>(ncc_);
      v_akt_(iblock, 1, k, j, i) += v_wk3_(k, j, i) * 1.0e-4 / 
          static_cast<float>(ncc_);
#ifdef BCKMEX
      v_akmt_(iblock, k, j, i)   += diff_back[iblock][j][i] * 
          10.0 * 1.0e-4;
      v_akt_(iblock, 0, k, j, i) += diff_back[iblock][j][i] / 
          static_cast<float>(ncc_) * 1.0e-4;
      v_akt_(iblock, 1, k, j, i) += diff_back[iblock][j][i] / 
          static_cast<float>(ncc_) * 1.0e-4;
#endif // BCKMEX
#endif // TIDEMIX
    }
    return ;
  };
 private:
  const int ncc_ = CppPconstMod::ncc;
  const ViewDouble3D v_wk1_     = *p_v_wk1;
  const ViewDouble3D v_wk2_     = *p_v_wk2;
  const ViewDouble3D v_wk3_     = *p_v_wk3;
  const ViewDouble4D v_akmt_    = *p_v_akmt;
  const ViewDouble4D v_ak_tide_ = *p_v_ak_tide;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble5D v_akt_     = *p_v_akt;
};

class FunctorReadyc5 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    if (v_vit_(0, 0, j, i) > 0.5) {
      const int iblock = 0;
// #ifdef __sw_slave__
//       double* wk1 = (double*) ldm_malloc (sizeof(double) * (KM-1));
//       double* wk2 = (double*) ldm_malloc (sizeof(double) * (KM-1));
//       double* wk3 = (double*) ldm_malloc (sizeof(double) * (KM-1));

//       double* wp1 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp2 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp3 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp4 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp5 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp6 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp7 = (double*) ldm_malloc (sizeof(double) * KM);
//       double* wp8 = (double*) ldm_malloc (sizeof(double) * KM);

//       double* akm_back = (double*) ldm_malloc (sizeof(double) * (KM-1));
//       double* akt_back = (double*) ldm_malloc (sizeof(double) * (KM-1));
//       double* aks_back = (double*) ldm_malloc (sizeof(double) * (KM-1));
// #else
      double wk1[KM-1], wk2[KM-1], wk3[KM-1];
      double wp1[KM], wp3[KM], wp4[KM];
      double wp5[KM], wp6[KM], wp7[KM], wp8[KM];
// #endif // __sw_slave__

      for (int k = 0; k < KM; ++k) {
        wp1[k] = 0.0;
        wp3[k] = 0.0;
        wp4[k] = 0.0;
        wp5[k] = 0.0;
        wp6[k] = 0.0;
        wp7[k] = 0.0;
        wp8[k] = 0.0;
      }
      const int kmt_m1 = v_kmt_(iblock, j, i) - 1;
      for (int k = 0; k < kmt_m1; ++k) {
        wp8[k] = - v_vit_(iblock, k+1, j, i) * v_zkp_(k+1) * 1.0e+2;
      }

      for (int k = 0; k < kmt_m1; ++k) {
        if (v_vit_(iblock, k+1, j, i) > 0.0) {
          const double tmp = v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1));
          wp1[k] = v_at_(iblock, 0, k, j, i) - (v_at_(iblock, 0, k, j, i) - v_at_(iblock, 0, k+1, j, i)) * tmp;
          wp4[k] = v_rit_(iblock, k, j, i);
          wp5[k] = v_ricdttms_(iblock, k, j, i);
          wp6[k] = v_s2t_(iblock, k, j, i);
          wp7[k] = v_rict_(iblock, k, j, i);

          wp3[k] = (v_pdensity_(iblock, k, j, i) + (v_pdensity_(iblock, k+1, j, i) -  v_pdensity_(iblock, k, j, i)) * tmp) * 1.e-3;
        }

        // wp1[k] = v_vit_(iblock, k+1, j, i) * 
        //     (v_at_(iblock, 0, k, j, i) - (v_at_(iblock, 0, k, j, i) - v_at_(i, j, k+1, 0, iblock)) * 
        //         v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1)));
        // wp4[k] =  v_vit_(iblock, k+1, j, i) * v_rit_(iblock, k, j, i);
        // wp5[k] =  v_vit_(iblock, k+1, j, i) * v_ricdttms_(iblock, k, j, i);
        // wp6[k] =  v_vit_(iblock, k+1, j, i) * v_s2t_(iblock, k, j, i);
        // wp7[k] =  v_vit_(iblock, k+1, j, i) * v_rict_(iblock, k, j, i);

        // wp3[k] =  v_vit_(iblock, k+1, j, i) * (v_pdensity_(iblock, k, j, i) + 
        //     (v_pdensity_(iblock, k+1, j, i) -  v_pdensity_(iblock, k, j, i)) * 
        //         v_dzp_(k) / (v_dzp_(k) + v_dzp_(k+1))) * 1.e-3;
      }
      // TODO vit can reduce
      double wp10 = v_vit_(iblock, 0, j, i) * v_buoytur_(iblock, j, i) * 1.0e+4;
      double wp11 = v_vit_(iblock, 0, j, i) * v_buoysol_(iblock, j, i) * 1.0e+4;

#ifdef CANUTOMIXOUT
      // To do something
#endif // CANUTOMIXOUT
#ifdef CANUTO2010
      // Fortran
      call canuto_2010_interface(wk1, wk2, wk3, wk4, 
          amld(i,j,iblock) ,tau_mag, wp1, wp2, wp3, wp4, wp5, 
          wp7, wp6, ulat(i,j,iblock)/degtorad , wp8, 
          kmt(i,j,iblock), iblock, j, i, isc)
#endif // CANUTO2010
      double mldtmp;
// #ifdef __sw_slave__
// #else
      double akm_back[KM-1];
      double akt_back[KM-1];
      double aks_back[KM-1];
// #endif

      turb_2(wp8, wp1, wp3, wp4, wp5, wp6, 
          dfricmx_ * 1.0e+4, dwndmix_ * 1.0e+4, 
              akm_back, akt_back, aks_back,
                  wp7, wp10, wp11, v_fcort_(iblock, j, i), 
                      mldtmp, wk1, wk2, wk3, kmt_m1, KM-1, 1, 0);

      v_amld_(iblock, j, i) = mldtmp * 1.0e-2;

#ifdef TIDEMIX
      for (int k = 0; k < KM; ++k) {
        v_ak_tide_(iblock, k, j, i) = 0.0;
      }
      for (int k = 0; k < kmt_m1; ++k) {
        v_ak_tide_(iblock, k, j, i) = BACK_TIDALMIXING + MIXING_EF * 
            LOCAL_MIXING_FRACTION * v_wave_dis_(iblock, j, i) * 
            v_fz_tide_(iblock, k, j, i) /
            (fmax(v_rict_(iblock, k, j, i), 1.0e-8) * wp3[k] * 1000.0);

        v_ak_tide_(iblock, k, j, i) = fmin(v_ak_tide_(iblock, k, j, i),
            MAX_TIDALMIXING);

        v_richardson_(iblock, k, j, i) = v_rict_(iblock, k, j, i);
        v_fztidal_(iblock, k, j, i)    = v_fz_tide_(iblock, k, j, i);
        v_wp3_tidal_(iblock, k, j, i)  = wp3[k];
      }
#ifdef CANUTOMIXOUT
#endif // CANUTOMIXOUT
      for (int k = kmt_m1 - 2; k >= 0; --k) {
        v_ak_tide_(iblock, k, j, i) = fmin(v_ak_tide_(iblock, k, j, i),
            v_ak_tide_(iblock, k+1, j, i));
      }
#endif // TIDEMIX
      const float tmp_ncc = 1.0 / static_cast<float>(ncc_);

      for (int k = 0; k < KM - 1; ++k) {
#ifdef TIDEMIX
        v_akmt_(iblock, k, j, i) = wk1[k] * 1.0e-4 +
            v_ak_tide_(iblock, k, j, i) * 5.0;
        // v_akt_(iblock, 0, k, j, i) += (wk2[k] * 1.0e-4 + 
        //     v_ak_tide_(iblock, k, j, i)) / 
        //         static_cast<float>(ncc_);
        // v_akt_(iblock, 1, k, j, i) += (wk3[k] * 1.0e-4 + 
        //     v_ak_tide_(iblock, k, j, i)) / 
        //         static_cast<float>(ncc_);
        v_akt_(iblock, 0, k, j, i) += (wk2[k] * 1.0e-4 + 
            v_ak_tide_(iblock, k, j, i)) * tmp_ncc;
        v_akt_(iblock, 1, k, j, i) += (wk3[k] * 1.0e-4 + 
            v_ak_tide_(iblock, k, j, i)) * tmp_ncc;
#ifdef BCKMEX
        v_akmt_(iblock, k, j, i) += diff_back[iblock][j][i] * 
            10.0 * 1.0e-4;
        v_akt_(iblock, 0, k, j, i) += diff_back[iblock][j][i] /
                static_cast<float>(ncc) * 1.0e-4;
        v_akt_(iblock, 1, k, j, i) += diff_back[iblock][j][i] /
                static_cast<float>(ncc) * 1.0e-4;
#endif // BCKMEX
#else // TIDEMIX
        v_akmt_(iblock, k, j, i) = wk1[k] * 1.0e-4;
        v_akt_(iblock, 0, k, j, i) += wk2[k] * 1.0e-4 / 
            static_cast<float>(ncc_);
        v_akt_(iblock, 1, k, j, i) += wk3[k] * 1.0e-4 / 
            static_cast<float>(ncc_);
#ifdef BCKMEX
        v_akmt_(iblock, k, j, i)   += diff_back[iblock][j][i] * 
            10.0 * 1.0e-4;
        v_akt_(iblock, 0, k, j, i) += diff_back[iblock][j][i] / 
            static_cast<float>(ncc_) * 1.0e-4;
        v_akt_(iblock, 1, k, j, i) += diff_back[iblock][j][i] / 
            static_cast<float>(ncc_) * 1.0e-4;
#endif // BCKMEX
#endif // TIDEMIX
      }
// #ifdef __sw_slave__
//       ldm_free (wk1, sizeof(double) * (KM-1));
//       ldm_free (wk2, sizeof(double) * (KM-1));
//       ldm_free (wk3, sizeof(double) * (KM-1));

//       ldm_free (wp1, sizeof(double) * KM);
//       ldm_free (wp2, sizeof(double) * KM);
//       ldm_free (wp3, sizeof(double) * KM);
//       ldm_free (wp4, sizeof(double) * KM);
//       ldm_free (wp5, sizeof(double) * KM);
//       ldm_free (wp6, sizeof(double) * KM);
//       ldm_free (wp7, sizeof(double) * KM);
//       ldm_free (wp8, sizeof(double) * KM);

//       ldm_free (akm_back, sizeof(double) * (KM-1));
//       ldm_free (akt_back, sizeof(double) * (KM-1));
//       ldm_free (aks_back, sizeof(double) * (KM-1));
// #endif // __sw_slave__
    }
    return ;
  };
 private:
  const int ncc_     = CppPconstMod::ncc;
  const double dfricmx_ = CppPconstMod::dfricmx;
  const double dwndmix_ = CppPconstMod::dwndmix;

  const double b1_         = CppCanutoMod::b1;
  const double dri_        = CppCanutoMod::dri;
  const double rri_        = CppCanutoMod::rri;
  const double rnd2on2_    = CppCanutoMod::rnd2on2;
  const double dand2on2_   = CppCanutoMod::dand2on2;
  const double deltheta_r_ = CppCanutoMod::deltheta_r;
  const double theta_rcrp_ = CppCanutoMod::theta_rcrp;
  const double theta_rcrn_ = CppCanutoMod::theta_rcrn;

  const ViewInt1D    v_irimax_     = *p_v_irimax;
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble1D v_dzp_        = *p_v_dzp;
  const ViewDouble1D v_zkp_        = *p_v_zkp;
  const ViewDouble1D v_rib_        = *p_v_rib;
  const ViewDouble1D v_ridb_       = *p_v_ridb;
  const ViewDouble1D v_sma1_       = *p_v_sma1;
  const ViewDouble1D v_sha1_       = *p_v_sha1;
  const ViewDouble1D v_ssa1_       = *p_v_ssa1;
  const ViewDouble1D v_and2on2a1_  = *p_v_and2on2a1;
  const ViewDouble1D v_amtaun2a1_  = *p_v_amtaun2a1;
  const ViewDouble1D v_back_ra_r_  = *p_v_back_ra_r;
  const ViewDouble1D v_sm_r1_      = *p_v_sm_r1;
  const ViewDouble1D v_sh_r1_      = *p_v_sh_r1;
  const ViewDouble1D v_ss_r1_      = *p_v_ss_r1;
  const ViewDouble1D v_slq2_r1_    = *p_v_slq2_r1;
  const ViewDouble2D v_smb_        = *p_v_smb;
  const ViewDouble2D v_shb_        = *p_v_shb;
  const ViewDouble2D v_ssb_        = *p_v_ssb;
  const ViewDouble2D v_slq2b_      = *p_v_slq2b;
  const ViewDouble3D v_amld_       = *p_v_amld;
  const ViewDouble3D v_fcort_      = *p_v_fcort;
  const ViewDouble3D v_buoytur_    = *p_v_buoytur;
  const ViewDouble3D v_buoysol_    = *p_v_buoysol;
  const ViewDouble3D v_wave_dis_   = *p_v_wave_dis;
  const ViewDouble4D v_rit_        = *p_v_rit;
  const ViewDouble4D v_rict_       = *p_v_rict;
  const ViewDouble4D v_akmt_       = *p_v_akmt;
  const ViewDouble4D v_fztidal_    = *p_v_fztidal;
  const ViewDouble4D v_wp3_tidal_  = *p_v_wp3_tidal;
  const ViewDouble4D v_pdensity_   = *p_v_pdensity;
  const ViewDouble4D v_ricdttms_   = *p_v_ricdttms;
  const ViewDouble4D v_richardson_ = *p_v_richardson;
  const ViewDouble4D v_ak_tide_    = *p_v_ak_tide;
  const ViewDouble4D v_fz_tide_    = *p_v_fz_tide;
  const ViewDouble4D v_s2t_        = *p_v_s2t;
  const ViewDouble4D v_vit_        = *p_v_vit;
  const ViewDouble5D v_at_         = *p_v_at;
  const ViewDouble5D v_akt_        = *p_v_akt;

  KOKKOS_INLINE_FUNCTION void turb_2(
      const double (&z)[KM],        // wp8
      const double (&t)[KM],        // wp1
      const double (&rh)[KM],       // wp3
            double (&ri)[KM],       // wp4
      const double (&rid)[KM],      // wp5
            double (&s2)[KM],       // wp6
      const double &fricmx,         // DFRICMX * 1.0e+4
      const double &wndmix,         // DWNDMIX * 1.0e+4
            double (&v_back)[KM-1], // akm_back
            double (&t_back)[KM-1], // akt_back
            double (&s_back)[KM-1], // aks_back
      const double (&an2)[KM],      // wp7
      const double &buoytur,        // wp10
      const double &buoysol,        // wp11
      const double &coriol,         // fcort[iblock][j][i]
            double &amld,           // mldtmp
            double (&akm)[KM-1],    // wk1
            double (&akh)[KM-1],    // wk2
            double (&aks)[KM-1],    // wk3
      const int    &n,              // iwk
      const int    &nmax,           // km - 1
      const int    &isurfuse,       // 1
      const int    &ifextermld      // 0
  ) const {

    double an = 0.0;
    double sh_back = 0.0;
    double sm_back = 0.0;
    double ss_back = 0.0;
    double slq2_back = 0.0;
 
    double ri1  = 0.0;
    double rid1 = 0.0;
 
    int ifbelow;
    int ifpureshear = 0;
 
    int ifnofsmall = 0;
 
    int nb = std::max(0, n);
    double visc_cbu_limit1 = fricmx;
    double diff_cbt_limit1 = fricmx;
    double buoytot = buoytur + buoysol;
 
    // if (ifextermld == 0) {
    //   if (n > 0) {
    //       formld(z, t, amld, n);
    //   } else if (n == 0) {
    //     amld = z[0];
    //   } else {
    //     amld = 0.0;
    //   }
    // }
    // if (ifextermld == 0) {
      if (n > 0) {
          formld(z, t, amld, n);
      } else if (n == 0) {
        amld = z[0];
      } else {
        amld = 0.0;
      }
    // }
 
    double al0 = 0.17 * amld;
 
    ifbelow = 0;
 
    // int icall(0);
    // int ipoint(0);
    // int iproblem(0)
    // int inegproblem(0);
 
    double sm, sh, ss;
    double slq2, epson2;

    double epsy[KM];
    bool lifepsy[KM];

    // !020219D REFERENCE NOTE: START OF PRIMARY LOOP OVER DEPTH LEVELS.
    for (int k = 0; k < n; ++k) {
  
      double amtaun2(0.0);
  
      ifnofsmall = 0;
      if (an2[k] >= 0.0) {
        an = sqrt(an2[k]);
      }
      if ((an / fabs(coriol)) < 1.0) {
        ifnofsmall = 1;
      }
  
      ri1 = ri[k];
  
      // if (ifsali == 1), ifsali = 1 in canuto_mod_2002
      // 981104 Use Ri_d = Ri_C - Ri_T in salinity-temperature turbulence model.
      rid1 = rid[k];
      const double and2 = rid[k] / (ri[k] + 1.0e-25) * an2[k];
      double and2on2 = and2 / (an2[k] + 1.0e-25);
  
      if (an2[k] < v_rib_(0) * s2[k]) {
        ifpureshear = 1;
      } else {
        ifpureshear = 0;
      }

      // ifastexpabs = 1 in canuto_mod_2002
      if (ifpureshear == 1) {
        const int imax = MT;
  
        interp1d_expabs(and2on2, amtaun2, 
            sm, sh, ss, imax, dand2on2_, rnd2on2_);
  
        slq2 = (-amtaun2) / ((b1_ * b1_) * ri[k] + 1.0e-25);
      }

      if (ifpureshear != 1) {
        interp2d_expabs(ri1, rid1, slq2, sm, sh, ss);
      }
  
      if (ifpureshear != 1) {
        if (slq2 < 0.0) {
          return ;
        }
        if (slq2 == 0) {
          ifbelow = 1;
        }
      }
  
      // bool lifupper = (((IFEPSON2  < 2) || (ifbelow == 0) || 
      //     ((IFDEEPLAT > 0) && (ifnofsmall == 1)) || 
      //     ((ri1 < 0.0) && (k <= 2))) && (slq2 > 0.0));
      bool lifupper = (((ifbelow == 0) || (ifnofsmall == 1) || 
          ((ri1 < 0.0) && (k <= 2))) && (slq2 > 0.0));
  
      // ++ipoint;
      // if (k == 0) {
      //   ++icall;
      // }
      epsy[k] = 0.0;
      lifepsy[k] = false;
  
      // if ((isurfuse == 1) && n > 0) {
      if (n > 0) {
        if (slq2 == 0.0) {
          epsy[k] = 0.0;
        } else {
          epsy[k] = - buoytot
              / ((1.0 / ((b1_ * b1_) * slq2)) - 0.5 * sm);
        }
        lifepsy[k] = ((epsy[k] >= 0.0) && lifupper);
        // if ((epsy[k] < 0.0) && lifupper) {
        //   ++iproblem;
        //   if (ri1 < 0.0) {
        //     ++inegproblem;
        //   }
        // }
      }

      // if ilomega == 1, ilomega = 0 in canuto_mod_2002
  
      double akz = 0.4 * z[k];
      double al = akz * al0 / (al0 + akz);
  
      double al2 = al * al;
  
      // if (!(((IFEPSON2 == 2) && (ifbelow == 1)) || lifepsy[k])) {
      if (!(ifbelow == 1 || lifepsy[k])) {
        if (ri1 > 0.0) {
          const double anlq2 = slq2 * ri1;
          if (anlq2 > 0.281) {
            al2  = 0.281 / anlq2 * al2;
            // if icondear == 0, icondear = 0 in canuto_mod_2002
            slq2 = 0.281 / (ri1 + 1.0e-20);
          }
        }
      }
  
      double epson2_;
      if (an2[k] < 0.0) {
        // Just to "hold the place". Value irrelevent here.
        epson2_ = EPSON2__;
      } else {
        double eplatidepend;
        // if ifdeeplat, ifdeeplat = 1 in canuto_mod_2002
        if (ifnofsmall == 1) {
          eplatidepend = 0.0;
        } else {
	 	      eplatidepend = eplatidepend_(fabs(coriol), an);
        }
  
        eplatidepend = fmax(eplatidepend, EPLATIDEPENDMIN);
  
        epson2_ = EPSON2__ * eplatidepend;
      }
  
      // if ifepson2 >= 1, ifepson2 = 2 in canuto_mod_2002
      // if ifbotenhance == 0, ifbotenhance = 0 in canuto_mod_2002
      epson2 = epson2_;
  
      // if ifback >= 4 and ifsali == 0, ifback = 5 and ifsali = 1 in canuto_mod_2002
      // if ifsali > 0, ifsali = 1 in canuto_mod_2002
      // if ifsalback == 1, 2, 3, ifsaliback = 5 in canuto_mod_2002
      // if ifsalback >= 4, ifsaliback = 5 in canuto_mod_2002
      // if ifsalback == 4, ifsaliback = 5 in canuto_mod_2002
      // else

      int ifrafglt = 0;
  
      double rit(0.0), ric(0.0);
      double back_ra_r1;
      double back_ri1,  back_ric1;
      double back_rid1, back_rit1;
      double theta_r(0.0);
      double deltheta_r1;
      int itheta_r0(0), itheta_r1(0);
      int jtheta_r0(0);
      double ra_r;
  
      const double tmp_deltheta_r = 1.0 / deltheta_r_;
      if (ri[k] <= 0.0) {
        back_ra_r1 = 0.0;
        back_rit1  = 0.0;
        back_ric1  = 0.0;
        back_ri1   = 0.0;
        back_rid1  = 0.0;
      } else {
        // rit = (ri[k] + rid[k]) / 2.0;
        // ric = (ri[k] - rid[k]) / 2.0;
        rit = (ri[k] + rid[k]) * 0.5;
        ric = (ri[k] - rid[k]) * 0.5;
        ra_r = sqrt((rit * rit) + (ric * ric));
  
        if (rit == 0.0) {
          if (ric == 0.0) {
            theta_r = atan(1.0);
          } else {
            // theta_r = PI / 2.0;
            theta_r = PI * 0.5;
          }
        } else {
          theta_r = atan(ric / rit);
        }
        // if (fabs(theta_r) > (PI / 2.0)) {
        if (fabs(theta_r) > (PI * 0.5)) {
          return ;
        }
        // if (theta_r < (-PI) / 4.0) {
        if (theta_r < (-PI) * 0.25) {
          theta_r += PI;
        }
  
        jtheta_r0 = static_cast<int>(
            // (theta_r + (PI / 4.0)) / deltheta_r_);
            (theta_r + (PI * 0.25)) * tmp_deltheta_r);
  
        itheta_r0 = jtheta_r0 - N_THETA_R_OCT;
        itheta_r1 = itheta_r0 + 1;
  
        double theta_r0 = itheta_r0 * deltheta_r_;
        double theta_r1 = itheta_r1 * deltheta_r_;
  
        if ((theta_r0 <= theta_rcrp_) &&
             (theta_r > theta_rcrp_)) {
          theta_r    = theta_r1;
          theta_r0   = theta_r1;
          itheta_r0  = itheta_r1;
          itheta_r1 += 1;
          theta_r1  += deltheta_r_;
        } else if ((theta_r1 >= theta_rcrn_) &&
                   (theta_r  < theta_rcrn_)) {
          theta_r    = theta_r0;
          theta_r1   = theta_r0;
          itheta_r1  = itheta_r0;
          itheta_r0 -= 1;
          theta_r0  -= deltheta_r_;
        }
  
        if ((itheta_r1 > 3 * N_THETA_R_OCT) ||
            (itheta_r0 < - N_THETA_R_OCT)) {
          return ;
        }
  
        deltheta_r1 = theta_r - theta_r0;
  
        const double delback_ra_r = 
            v_back_ra_r_(N_THETA_R_OCT + itheta_r1) 
                - v_back_ra_r_(N_THETA_R_OCT + itheta_r0);
  
        const double dback_ra_r_o_dtheta = 
            // delback_ra_r / deltheta_r_;
            delback_ra_r * tmp_deltheta_r;
  
        back_ra_r1 = v_back_ra_r_(N_THETA_R_OCT + itheta_r0) 
            + deltheta_r1 * dback_ra_r_o_dtheta;
  
        ifrafglt = 0;

        // if ifrafgmax == 1, ifrafgmax = 1 in canuto_mod_2002
  
        if ((theta_r <= theta_rcrp_) ||
            (theta_r >= theta_rcrn_)) {
          if (back_ra_r1 > ra_r) {
            ifrafglt = 1;
            back_ra_r1 = ra_r;
          }
        }
  
        if (back_ra_r1 < 0.0) {
          return ;
        }
        back_rit1 = cos(theta_r) * back_ra_r1;
        back_ric1 = sin(theta_r) * back_ra_r1;
        back_ri1  = back_rit1 + back_ric1;
        back_rid1 = back_rit1 - back_ric1;
      }
      // End if ifsalback > 4
  
      if (ri[k] > 0.0) {
  
        // ifbg_theta_interp = 1 in canuto_mod_2002
        // if ((IFBG_THETA_INTERP == 0) ||
        //     (ifrafglt == 1)) {
        if (ifrafglt == 1) {
          // if ifastexpabs == 1, ifastexpabs = 1 in canuto_mod_2002
          interp2d_expabs(back_ri1, back_rid1,
              slq2_back, sm_back, sh_back, ss_back);

        // ifbg_theta_interp = 1 in canuto_mod_2002
        // } else if(IFBG_THETA_INTERP == 1) {
        } else if(true) {
          // !000315 Interpolate 1D table of background vs. theta_r instead.
          deltheta_r1 = theta_r - itheta_r0 * deltheta_r_;
  
          const double delsm_back = v_sm_r1_(N_THETA_R_OCT + itheta_r1)
              - v_sm_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dsm_back_o_dtheta = delsm_back / deltheta_r_;
          const double dsm_back_o_dtheta = delsm_back * tmp_deltheta_r;
  
          sm_back = v_sm_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dsm_back_o_dtheta;
  
          const double delsh_back = v_sh_r1_(N_THETA_R_OCT + itheta_r1)
              - v_sh_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dsh_back_o_dtheta = delsh_back / deltheta_r_;
          const double dsh_back_o_dtheta = delsh_back * tmp_deltheta_r;
  
          sh_back = v_sh_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dsh_back_o_dtheta;
  
          const double delss_back = v_ss_r1_(N_THETA_R_OCT + itheta_r1)
              - v_ss_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dss_back_o_dtheta = delss_back / deltheta_r_;
          const double dss_back_o_dtheta = delss_back * tmp_deltheta_r;
  
          ss_back = v_ss_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dss_back_o_dtheta;
  
          const double delslq2_back = v_slq2_r1_(N_THETA_R_OCT + itheta_r1)
              - v_slq2_r1_(N_THETA_R_OCT + itheta_r0);
  
          // const double dslq2_back_o_dtheta = delslq2_back / deltheta_r_;
          const double dslq2_back_o_dtheta = delslq2_back * tmp_deltheta_r;
  
          slq2_back = v_slq2_r1_(N_THETA_R_OCT + itheta_r0)
              + deltheta_r1 * dslq2_back_o_dtheta;
        } else {
          return ;
        }
        if (slq2_back < 0.0) {
          return ;
        }
      } // not go to 19

      // 19 continue
      if (ri1 < 0.0) {
        sm_back = 0.0;
        sh_back = 0.0;
        ss_back = 0.0;
      }
  
      if ((sm_back < 0.0) || (sh_back < 0.0) ||
          (ss_back < 0.0)) {
        return ;
      }
      
      //double tmp_back;

      // if ifepson2 == 0, ifepson2 = 2 in canuto_mod_2002
      // if ifepson2 >  0, ifepson2 = 2 in canuto_mod_2002
      const double tmp_back = 0.5 * b1_ * b1_
          * back_ri1 * slq2_back * epson2;
  
      v_back[k] = tmp_back * sm_back;
      t_back[k] = tmp_back * sh_back;
      s_back[k] = tmp_back * ss_back;
  
      if ((v_back[k] < 0.0) || (t_back[k] < 0.0)
          || (s_back[k] < 0.0)) {
        return ;
      }
      if((ri[k] > 0.0) && ((v_back[k] == 0.0)
          || (t_back[k] == 0.0) || (s_back[k] == 0.0))) {
        return ;
      }
  
      double aldeep[NBIG];
  
      aldeep[k] = 0.0;
  
      double tmp(0.0);
  
      // if ifepson2 == 2, ifepson2 = 2 in canuto_mod_2002
      // if ((IFEPSON2 == 2) && (ifbelow == 1)) {
      //   if ((ri1 >= 0.0) || ((IFDEEPLAT == 2) && (ifnofsmall == 1))) {
      if (ifbelow == 1) {
        if (ri1 >= 0.0) {
          tmp = 0.5 * b1_ * b1_ * ri1 * slq2 * epson2;
        } else if (k > 1) {
          double delz, delrh, del2rh;
          if (k == n - 1) {
            delz   = z [k] - z [k-1];
            delrh  = rh[k] - rh[k-1];
            del2rh = rh[k] - 2.0 * rh[k-1] + rh[k-2];
          } else {
            delz   = z [k+1] - z [k-1];
            delrh  = rh[k+1] - rh[k-1];
            del2rh = rh[k+1] - 2.0 * rh[k] + rh[k-1];
          }
  
          const double dzrh = delrh / delz;
          const double d2zrh = 4.0 * del2rh / (delz * delz);
          const double rdzlndzrh = dzrh / d2zrh;
          const double al0deep = 0.17 * fabs(rdzlndzrh);
          akz = 0.4 * z[k];
          aldeep[k] = akz * al0deep / (al0deep + akz);
          al2 = aldeep[k] * aldeep[k];
  
          if (ifpureshear == 1) {
            // GO TO 21
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
       
            akm[k] = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            akh[k] = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            aks[k] = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
          } else if (IFSHEARMIN) {
            s2[k] = fmax(s2[k], S2MIN);
          }
          tmp = 0.5 * b1_ * al2
              * sqrt(s2[k] / (slq2 + 1.0e-40));
        } else {
          if (ifpureshear == 1) {
            // GO TO 21
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
       
            akm[k] = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            akh[k] = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            aks[k] = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
          } else if (IFSHEARMIN) {
            s2[k] = fmax(s2[k], S2MIN);
          }
          if (lifepsy[k]) {
            tmp = 0.5 * epsy[k] / (s2[k] + 1.0e-40);
          } else {
            tmp = 0.5 * b1_ * al2
                * sqrt(s2[k] / (slq2 + 1.0e-40));
          }
        }
      } else {
        if (ifpureshear == 1) {
            // GO TO 21
            tmp = 0.5 * (b1_ * b1_) * al2
                * sqrt(-an2[k] / amtaun2);
            akm[k] = fmin(tmp * sm + v_back[k],
                visc_cbu_limit1);
            akh[k] = fmin(tmp * sh + t_back[k],
                diff_cbt_limit1);
            aks[k] = fmin(tmp * ss + s_back[k],
                diff_cbt_limit1);
            continue ;
        } else if (IFSHEARMIN) {
          s2[k] = fmax(s2[k], S2MIN);
        }
        if (lifepsy[k]) {
          tmp = 0.5 * epsy[k] / (s2[k] + 1.0e-40);
        } else {
          tmp = 0.5 * b1_ * al2
              * sqrt(s2[k] / (slq2 + 1.0e-40));
        }
      }
      // 21
      if (ifpureshear == 1) {
        tmp = 0.5 * (b1_ * b1_) * al2
            * sqrt(-an2[k] / amtaun2);
      }
  
      akm[k] = fmin(tmp * sm + v_back[k],
          visc_cbu_limit1);
      akh[k] = fmin(tmp * sh + t_back[k],
          diff_cbt_limit1);
      aks[k] = fmin(tmp * ss + s_back[k],
          diff_cbt_limit1);
    } // End for k in 0:n

    for (int k = nb+1; k < nmax; ++k) {
      akm[k] = 0.0;
      akh[k] = 0.0;
      aks[k] = 0.0;
    }
 
    if (n > 0) {
      if (akm[0] < wndmix) {
        akm[0] = wndmix;
      }
      if (akh[0] < wndmix) {
        akh[0] = wndmix;
      }
      if (aks[0] < wndmix) {
        aks[0] = wndmix;
      }
    }
    return ;
  }

  KOKKOS_INLINE_FUNCTION void formld(
      const double (&z)[KM], const double (&t)[KM],
          double &amld, const int &n) const {
    for (int k = 0; k < n; ++k) {
      if (fabs(t[k] - t[0]) > 0.1) {
#ifdef D_PRECISION
        const double tm = t[0] - sign_double(0.1, t[0] - t[k]);
#else
        const double tm = t[0] - sign_float(0.1,   t[0] - t[k]);
#endif // D_PRECISION
        amld = z[k] + (z[k-1] - z[k]) *
            (tm - t[k]) / (t[k-1] - t[k] + 1.e-20);
        return ;
      }
    }
    amld = z[n-1];
    return ;
  }

  KOKKOS_INLINE_FUNCTION int sign_int(const double &x, const double &y) const {
    return y >= 0 ? std::abs(static_cast<int>(x)) 
        : - std::abs(static_cast<int>(x));
  }

  KOKKOS_INLINE_FUNCTION int sign_float(const double &x, const double &y) const {
    return y >= 0.0 ? std::abs(static_cast<float>(x)) 
       : - std::abs(static_cast<float>(x));
  }

  KOKKOS_INLINE_FUNCTION double sign_double(const double &x, const double &y) const {
    return y >= 0.0 ? fabs(x) : - fabs(x);
  }

  KOKKOS_INLINE_FUNCTION void interp1d_expabs(
      double &x,            // and2on2
      double &slq2,         // amtaun2
      double &sm,           // sm
      double &sh,           // sh
      double &ss,           // ss
      const int &ixmax,     // imax
      const double &delta,  // dand2on2
      const double &rat     /* rnd2on2 */ ) const {
    // printf ("interp1d_expabs\n");
  
    int lx0(0), lx1(0);
    
    if (x > v_and2on2a1_(MT + ixmax)) {
      x = v_and2on2a1_(MT + MT);
    } else if (x < v_and2on2a1_(0)) {
      x = v_and2on2a1_(0);
    }
  
    if (fabs(x) < v_and2on2a1_(MT + MT0)) {
#ifdef D_PRECISION
      lx1 = static_cast<int>(x / delta) + round(sign_double(static_cast<double>(1.0), x)); 
#else  // D_PRECISION
      lx1 = static_cast<int>(x / delta) + round(sign_float(static_cast<float>(1.0), x)); 
#endif // D_PRECISION
    } else if (fabs(x) >= v_and2on2a1_(MT + MT)) {
#ifdef D_PRECISION
      lx0 = round(sign_double(static_cast<double>(MT), x));
#else  // D_PRECISION
      lx0 = round(sign_float(static_cast<float>(MT), x));
#endif // D_PRECISION
      lx1 = lx0;
    } else {
#ifdef D_PRECISION
      const double tabindx = sign_double(static_cast<double>(MT0)
         + ((log(fabs(x)) - log(v_and2on2a1_(MT + MT0))) / log(rat)), x);
#else  // D_PRECISION
      const float tabindx = sign_double(static_cast<float>(MT0)
         + ((log(fabs(x)) - log(v_and2on2a1_(MT + MT0))) / log(rat)), x);
#endif // D_PRECISION

#ifdef D_PRECISION
      lx1 = static_cast<int>(tabindx) + round(sign_double(
          static_cast<double>(1.0), x));
#else  // D_PRECISION
      lx1 = static_cast<int>(tabindx) + round(sign_float(
          static_cast<float>(1.0), x));
#endif // D_PRECISION
    }

    if (!(fabs(x) >= v_and2on2a1_(MT + MT))) {
      if (fabs(v_and2on2a1_(MT + lx1)) < fabs(x)) {
        lx1 += sign_int(1, lx1);
      } else if (fabs(v_and2on2a1_(MT + lx1 - sign_int(1, lx1))) > fabs(x)) {
        lx1 -= sign_int(1, lx1);
      }

#ifdef D_PRECISION
      lx0 = lx1 - round(sign_double(
          static_cast<double>(1.0), x));
#else  // D_PRECISION
      lx0 = lx1 - round(sign_float(
          static_cast<float>(1.0), x));
#endif // D_PRECISION
      if (x == 0.0) {
        lx1 = 1;
      }
    }
    if ((x > 0.0 && (x < v_and2on2a1_(MT + lx0) || x > v_and2on2a1_(MT + lx1))) ||
        (x < 0.0 && (x > v_and2on2a1_(MT + lx0) || x < v_and2on2a1_(MT + lx1)))) {
      return;
    }

    const double deltaxta = 1.0 / (v_and2on2a1_(MT + lx1) - v_and2on2a1_(MT + lx0));

    const double deltax = x - v_and2on2a1_(MT + lx0);

    double dslq2_x;
    if (lx1 == lx0) {
      dslq2_x = 0.0;
    } else {
      dslq2_x = (v_amtaun2a1_(MT + lx1) - v_amtaun2a1_(MT + lx0)) * deltaxta;
    }

    slq2 = v_amtaun2a1_(MT + lx0) + dslq2_x * deltax;

    double dsm_x;
    if (lx1 == lx0) {
      dsm_x = 0.0;
    } else {
      dsm_x = (v_sma1_(MT + lx1) - v_sma1_(MT + lx0)) * deltaxta;
    }
    sm = v_sma1_(MT + lx0) + dsm_x * deltax;

    double dsh_x;
    if (lx1 == lx0) {
      dsh_x = 0.0;
    } else {
      dsh_x = (v_sha1_(MT + lx1) - v_sha1_(MT + lx0)) * deltaxta;
    }
    sh = v_sha1_(MT + lx0) + dsh_x * deltax;

    double dss_x;
    if (lx1 == lx0) {
      dss_x = 0.0;
    } else {
      dss_x = (v_ssa1_(MT + lx1) - v_ssa1_(MT + lx0)) * deltaxta;
    }
    ss = v_ssa1_(MT + lx0) + dss_x * deltax;
  
    return ;  
  }

  KOKKOS_INLINE_FUNCTION void interp2d_expabs(double &ri, double &rid,
      double &slq2, double &sm, double &sh, double &ss) const {

    const double tmp_ri  = 1.0 /ri;
    const double tmp_rid = 1.0 /rid;

    if (ri > v_rib_(MT + MT)) {
      if (fabs(rid) <= ri) {
        // rid = v_rib_(MT + MT) * (rid / ri);
        rid = v_rib_(MT + MT) * (rid * tmp_ri);
        ri  = v_rib_(MT + MT);
      } else if (rid > ri) {
        // ri  = v_ridb_(MT + MT) * (ri / rid);
        ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
        rid = v_ridb_(MT + MT);
      } else if (rid < -ri) {
        // ri  = v_ridb_(0) * (ri / rid);
        ri  = v_ridb_(0) * (ri * tmp_rid);
        rid = v_ridb_(0);
      }
    } else if (ri < v_rib_(0)) {
      if (fabs(rid) <= -ri) {
        // rid = v_rib_(0) * (rid / ri);
        rid = v_rib_(0) * (rid * tmp_ri);
        ri  = v_rib_(0);
      } else if (rid > -ri) {
        // ri  = v_ridb_(MT + MT) * (ri / rid);
        ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
        rid = v_ridb_(MT + MT);
      } else if (rid < ri) {
        // ri  = v_ridb_(0) * (ri / rid);
        ri  = v_ridb_(0) * (ri * tmp_rid);
        rid = v_ridb_(0);
      }
    } else if (rid > v_ridb_(MT + MT)) {
      // ri  = v_ridb_(MT + MT) * (ri / rid);
      ri  = v_ridb_(MT + MT) * (ri * tmp_rid);
      rid = v_ridb_(MT + MT);
    } else if (rid < v_ridb_(0)) {
      // ri  = v_ridb_(0) * (ri / rid);
      ri  = v_ridb_(0) * (ri * tmp_rid);
      rid = v_ridb_(0);
    }

    int lrid0(0), lrid1(0);

    if (fabs(rid) < v_ridb_(MT + MT0)) {
#ifdef D_PRECISION
      lrid1 = static_cast<int>(rid / dri_)
          + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid1 = static_cast<int>(rid / dri_)
          + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION
    } else if (fabs(rid) >= v_ridb_(MT + MT)) {
#ifdef D_PRECISION
      lrid0 = round(sign_double(static_cast<double>(MT), rid));
#else // D_PRECISION
      lrid0 = round(sign_float(static_cast<float>(MT), rid));
#endif // D_PRECISION
      lrid1 = lrid0;
    } else {
#ifdef D_PRECISION
      const double tabindrid = sign_double(static_cast<double>(MT0)
          + ((log(fabs(rid)) - log(v_ridb_(MT + MT0)))
              / log(rri_)), rid);
#else // D_PRECISION
      const double tabindrid = sign_float(static_cast<float>(MT0)
          + ((log(fabs(rid)) - log(v_ridb_(MT + MT0)))
              / log(rri)), rid);
#endif // D_PRECISION

#ifdef D_PRECISION
      lrid1 = static_cast<int>(tabindrid)
          + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid1 = static_cast<int>(tabindrid)
          + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION

    }

    if (!(fabs(rid) >= v_ridb_(MT + MT))) {

      if (fabs(v_ridb_(MT + lrid1)) < fabs(rid)) {
        lrid1 += sign_int(1, lrid1);
      } else if (fabs(v_ridb_(MT + lrid1 - sign_int(1, lrid1))) > 
          fabs(rid)) {
        lrid1 -= sign_int(1, lrid1);
      }

#ifdef D_PRECISION
      lrid0 = lrid1 - round(sign_double(
          static_cast<double>(1.0), rid));
#else // D_PRECISION
      lrid0 = lrid1 - round(sign_float(
          static_cast<float>(1.0), rid));
#endif // D_PRECISION
      if (rid == 0.0) {
        lrid1 = 1;
      }
    }

    // 252
    if ((rid > 0.0 && (rid < v_ridb_(MT + lrid0) || rid > v_ridb_(MT + lrid1))) ||
        (rid < 0.0 && (rid > v_ridb_(MT + lrid0) || rid < v_ridb_(MT + lrid1)))) {
      return ;
    }
 
    if (ri > fmin(v_rib_(MT + v_irimax_(MT + lrid0)),
                  v_rib_(MT + v_irimax_(MT + lrid1)))) {
      slq2 = 0.0;
      sm   = 0.0;
      sh   = 0.0;
      ss   = 0.0;
      return ;
    }

    int lri0(0), lri1(0);

    if (fabs(ri) < v_rib_(MT + MT0)) {
#ifdef D_PRECISION
      lri1 = static_cast<int>(ri / dri_)
          + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri1 = static_cast<int>(ri / dri_)
          + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
    } else if (fabs(ri) >= v_rib_(MT + MT)) {
#ifdef D_PRECISION
      lri0 = round(sign_double(static_cast<double>(MT), ri));
#else // D_PRECISION
      lri0 = round(sign_float(static_cast<float>(MT), ri));
#endif // D_PRECISION
      lri1 = lri0;
    } else {
#ifdef D_PRECISION
      const double tabindri = sign_double(static_cast<double>(MT0)
          + ((log(fabs(ri)) - log(v_rib_(MT + MT0)))
              / log(rri_)), ri);
#else // D_PRECISION
      const double tabindri = sign_float(static_cast<float>(MT0)
          + ((log(fabs(ri)) - log(v_rib_(MT + MT0)))
              / log(rri)), ri);
#endif // D_PRECISION

#ifdef D_PRECISION
      lri1 = static_cast<int>(tabindri)
          + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri1 = static_cast<int>(tabindri)
          + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
    }
 
    if (!(fabs(ri) >= v_rib_(MT + MT))) {
      if (fabs(v_rib_(MT + lri1)) < fabs(ri)) {
        lri1 += sign_int(1, lri1);
      } else if (fabs(v_rib_(MT + lri1 - sign_int(1, lri1))) >
          fabs(ri)) {
        lri1 -= sign_int(1, lri1);
      }

#ifdef D_PRECISION
      lri0 = lri1 - round(sign_double(
          static_cast<double>(1.0), ri));
#else // D_PRECISION
      lri0 = lri1 - round(sign_float(
          static_cast<float>(1.0), ri));
#endif // D_PRECISION

      if (ri == 0.0) {
        lri1 = 1;
      }
    }

    if ((ri > 0.0 && (ri < v_rib_(MT + lri0) || ri > v_rib_(MT + lri1))) ||
        (ri < 0.0 && (ri > v_rib_(MT + lri0) || ri < v_rib_(MT + lri1)))) {
      return ;
    }

    const double deltaridta = 1.0 / (v_ridb_(MT + lrid1) - v_ridb_(MT + lrid0));
    const double deltarita  = 1.0 / (v_rib_(MT + lri1) - v_rib_(MT + lri0));
    const double deltarid = rid - v_ridb_(MT + lrid0);
    const double deltari  = ri  - v_rib_(MT + lri0);

    double dslq2_rid;

    if (lrid1 == lrid0) {
      dslq2_rid = 0.0;
    } else {
      dslq2_rid = (v_slq2b_(MT+lrid1, MT+lri0) - v_slq2b_(MT+lrid0, MT+lri0))
          * deltaridta;
    }

    double dslq2_ri;
    if (lri1 == lri0) {
      dslq2_ri = 0.0;
    } else {
      dslq2_ri = (v_slq2b_(MT+lrid0, MT+lri1) - v_slq2b_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    slq2 = v_slq2b_(MT+lrid0, MT+lri0)
        + dslq2_ri * deltari + dslq2_rid * deltarid;

    double dsm_rid;
    if (lrid1 == lrid0) {
      dsm_rid = 0.0;
    } else {
      dsm_rid = (v_smb_(MT+lrid1, MT+lri0) - v_smb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }

    double dsm_ri;

    if (lri1 == lri0) {
      dsm_ri = 0.0;
    } else {
      dsm_ri = (v_smb_(MT+lrid0, MT+lri1) - v_smb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    sm = v_smb_(MT+lrid0, MT+lri0) 
        + dsm_ri * deltari + dsm_rid * deltarid;

    double dsh_rid;
    if (lrid1 == lrid0) {
      dsh_rid = 0.0;
    } else {
      dsh_rid = (v_shb_(MT+lrid1, MT+lri0) - v_shb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }
    double dsh_ri;
    if (lri1 == lri0) {
      dsh_ri = 0.0;
    } else {
      dsh_ri = (v_shb_(MT+lrid0, MT+lri1) - v_shb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    sh = v_shb_(MT+lrid0, MT+lri0) 
        + dsh_ri * deltari + dsh_rid * deltarid;
    
    double dss_rid;
    if (lrid1 == lrid0) {
      dss_rid = 0.0;
    } else {
      dss_rid = (v_ssb_(MT+lrid1, MT+lri0) - v_ssb_(MT+lrid0, MT+lri0))
          * deltaridta;
    }
    double dss_ri;
    if (lri1 == lri0) {
      dss_ri = 0.0;
    } else {
      dss_ri = (v_ssb_(MT+lrid0, MT+lri1) - v_ssb_(MT+lrid0, MT+lri0))
          * deltarita;
    }

    ss = v_ssb_(MT+lrid0, MT+lri0) 
        + dss_ri * deltari + dss_rid * deltarid;

    return ;
  }

  KOKKOS_INLINE_FUNCTION double acosh1(const double &x) const {
    return log(x + sqrt((x * x) - 1.0));
  }
  KOKKOS_INLINE_FUNCTION double wavelat(const double &xf, const double &yn) 
      const {
    return xf * acosh1(yn / xf);
  }

  KOKKOS_INLINE_FUNCTION double eplatidepend_(const double &f, const double &an) 
      const {
    double f_30;
    double anum, den;
    double pi1, omega;
    const double an0 = 5.24e-3;
 
    pi1   = 4.0 * atan(1.0);
    omega = pi1 / 43082.0e0;
 
    f_30  = omega;
 
    den   = wavelat(f_30, an0);
 
    anum  = wavelat(f, an);
    return anum / den;
  }

};

#endif // CANUTO

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
    // if (i >= 1 && j < (JMT-1)) {
      v_uk_(iblock, k, j, i) = (1.0 + v_work_(iblock, j, i) 
          * v_ohbu_(iblock, j, i)) * v_uwk(iblock, k, j, i);
      v_vk_(iblock, k, j, i) = (1.0 + v_work_(iblock, j, i) 
          * v_ohbu_(iblock, j, i)) * v_vwk(iblock, k, j, i);
    // }
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

    v_work_(iblock, j, i) = C0;

    if (i >= 1 && i < (IMT-1) && j >= 1 && j < (JMT-1)) {
      for (int k = 0; k < KM; ++k) {
        v_work_(iblock, j, i) -= v_dzp_(k)
            * v_wka_(iblock, k, j, i)
                * v_vit_(iblock, k, j, i);
      }
      double pre_ws = v_ws_(iblock, 0, j, i);
      for (int k = 1; k < KM; ++k) {
        v_ws_(iblock, k, j, i) = v_vit_(iblock, k, j, i) 
            * (pre_ws + v_dzp_(k-1) 
                * (v_work_(iblock, j, i) * v_ohbt_(iblock, j, i)
                    + v_wka_(iblock, k-1, j, i  )));
        pre_ws = v_ws_(iblock, k, j, i);
      }

      v_work_(iblock, j, i) = 1.0 / (1.0 + v_h0wk(iblock, j, i)
          * v_ohbt_(iblock, j, i));

      for (int k = 1; k < KM; ++k) {
        v_ws_(iblock, k, j, i) *= v_work_(iblock, j, i);
      }
    }
    return ;
  }
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_work_ = *p_v_work;
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
         + v_uuu(iblock, k, j, i  ) * v_dyu_(iblock, j, i  )) 
              * P25;
      v_uv_ws_face_(k, j, i, 1) = 
          (v_vvv(iblock, k, j  , i) * v_dxu_(iblock, j  , i)  
         + v_vvv(iblock, k, j+1, i) * v_dxu_(iblock, j+1, i))
              * P25;
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
    hdiffu_del4_1(k, j, i, v_up_, v_vp_);
    return ;
  }
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
};

class FunctorReadyc15 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {

    hdiffu_del4_2 (k, j, i, v_up_, v_vp_);

    return ;
  }
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
          + (v_dm_cnsew_ (bid, j, i, 0) * v_umixk(bid, k, j  , i  )
           + v_dm_cnsew_ (bid, j, i, 1) * v_umixk(bid, k, j-1, i  )
           + v_dm_cnsew_ (bid, j, i, 2) * v_umixk(bid, k, j+1, i  )
           + v_dm_cnsew_ (bid, j, i, 3) * v_umixk(bid, k, j  , i+1)
           + v_dm_cnsew_ (bid, j, i, 4) * v_umixk(bid, k, j  , i-1));
    }
    return ;
  }
 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ib;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const ViewDouble3D v_d2uk_      = *p_v_d2uk;
  const ViewDouble3D v_d2vk_      = *p_v_d2vk;
  const ViewDouble4D v_up_        = *p_v_up;
  const ViewDouble4D v_vp_        = *p_v_vp;
  const ViewDouble4D v_dm_cnsew_  = *p_v_dm_cnsew;
  const ViewDouble4D v_du_cnsewm_ = *p_v_du_cnsewm;
};
class FunctorReadyc16 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    hdiffu_del4_3 (iblock, k, j, i);

    return ;
  }
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
 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ib;
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
};
class FunctorReadyc17 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      double hduk, hdvk;

      hdiffu_del4_4 (iblock, k, j, i, hduk, hdvk);

      v_dlu_(iblock, k, j, i) += hduk;
      v_dlv_(iblock, k, j, i) += hdvk;
    }
    return ;
  }
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
         + (v_dm_cnsew_( bid, j, i, 0) * v_d2uk_(k, j  , i  )
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
 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ib;
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
    v_wk2(iblock, j, i) = C0;
    for (int k = 0; k < KM; ++k) {
      v_wk2(iblock, j, i) += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(iblock, k, j, i) *v_viv_(iblock, k, j, i);
    }
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
    v_wka_(iblock, k, j, i) = c0f_ * sqrt(
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
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    // d0 = 1026.0D0, od0 = 1 / d0
    const double od0 = 1.0 / 1026.0;
    // const double aidif = 0.5;
    // aidifm1 = 1.0 - aidif
    const double aidifm1 = 0.5;
    double diff_u1, diff_u2;
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
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    const double od0 = 1.0 / 1026.0;
    // const double aidif = 0.5;
    // aidifm1 = 1.0 - aidif
    const double aidifm1 = 0.5;
    double diff_v1, diff_v2;
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
KOKKOS_REGISTER_FOR_3D(FunctorReadyc2,  FunctorReadyc2)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc3,  FunctorReadyc3)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc4,  FunctorReadyc4)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc5,  FunctorReadyc5)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc51,  FunctorReadyc51)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc52,  FunctorReadyc52)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc53,  FunctorReadyc53)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc54,  FunctorReadyc54)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc55,  FunctorReadyc55)
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
KOKKOS_REGISTER_FOR_2D(FunctorReadyc18, FunctorReadyc18)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc19, FunctorReadyc19)
KOKKOS_REGISTER_FOR_2D(FunctorReadyc20, FunctorReadyc20)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc21, FunctorReadyc21)
KOKKOS_REGISTER_FOR_3D(FunctorReadyc22, FunctorReadyc22)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYC_HPP_
