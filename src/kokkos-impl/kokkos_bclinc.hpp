#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BCLINC_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BCLINC_HPP_
#include "../head/def-undef.h"

#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_work_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_grid.h"
#include "../head/kokkos_tmp_var.h" 
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

using CppParamMod::KM;  
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JET;
using CppParamMod::JST;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppPconstMod::isc;
using CppConstantMod::C0;
using CppConstantMod::G;
using CppConstantMod::P5;
using CppConstantMod::P25;

using KokkosPconstMod::p_v_snlat;
using KokkosPconstMod::p_v_vit;
using KokkosPconstMod::p_v_viv;
using KokkosPconstMod::p_v_dzp;
using KokkosPconstMod::p_v_epea;
using KokkosPconstMod::p_v_epeb;
using KokkosPconstMod::p_v_epla;
using KokkosPconstMod::p_v_eplb;
using KokkosPconstMod::p_v_odzp;
using KokkosPconstMod::p_v_odzt;
using KokkosPconstMod::p_v_ohbu;
using KokkosPconstMod::p_v_ohbt;
using KokkosPconstMod::p_v_akmu;
using KokkosPconstMod::p_v_zkt;
using KokkosGrid::p_v_kmu;
using KokkosGrid::p_v_fcor;
using KokkosGrid::p_v_dxyur;
using KokkosDynMod::p_v_sbcx;
using KokkosDynMod::p_v_sbcy;
using KokkosDynMod::p_v_bbcx;
using KokkosDynMod::p_v_bbcy;
using KokkosDynMod::p_v_up;
using KokkosDynMod::p_v_vp;
using KokkosDynMod::p_v_dlu;
using KokkosDynMod::p_v_dlv;
using KokkosDynMod::p_v_h0bf;
using KokkosDynMod::p_v_h0bl;
using KokkosDynMod::p_v_gg;
using KokkosDynMod::p_v_u;
using KokkosDynMod::p_v_v;
using KokkosDynMod::p_v_ub;
using KokkosDynMod::p_v_vb;
using KokkosDynMod::p_v_utf;
using KokkosDynMod::p_v_vtf;
using KokkosWorkMod::p_v_work;
using KokkosWorkMod::p_v_work_1;
using KokkosWorkMod::p_v_wka;
using KokkosForcMod::p_v_psa;
using KokkosForcMod::p_v_su;
using KokkosForcMod::p_v_sv;

using KokkosTmpVar::p_v_wp3;
using KokkosTmpVar::p_v_wp4;
using KokkosTmpVar::p_v_wp5;
using KokkosTmpVar::p_v_wp6;
using KokkosTmpVar::p_v_wk1;
using KokkosTmpVar::p_v_wk2;

using KokkosTmpVar::p_v_work1;
using KokkosTmpVar::p_v_work2;
using KokkosTmpVar::p_v_work_merged;

class FunctorBclinc1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {   
    const int iblock = 0;
#if (defined CANUTO) || (defined CANUTO2010)
    const int kmb = v_kmu_(iblock, j, i);
    if (kmb >= 1) {
      v_sbcx_(iblock, j, i) = v_su_(iblock, j, i) * od0_;
      v_sbcy_(iblock, j, i) = v_sv_(iblock, j, i) * od0_;
    }
#endif // CANUTO
    return ;
  }
 private:
  const double od0_ = CppPconstMod::od0;
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble3D v_su_   = *p_v_su;
  const ViewDouble3D v_sv_   = *p_v_sv;
  const ViewDouble3D v_sbcx_ = *p_v_sbcx;
  const ViewDouble3D v_sbcy_ = *p_v_sbcy;
};
class FunctorBclinc2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {   

    const int iblock = 0;
#if (defined CANUTO) || (defined CANUTO2010)
    const int kmb = v_kmu_(iblock, j, i);
    if (kmb >= 1) {
      const double tmp = c0f_ * sqrt(
          v_up_(iblock, kmb-1, j, i) * v_up_(iblock, kmb-1, j, i)
        + v_vp_(iblock, kmb-1, j, i) * v_vp_(iblock, kmb-1, j, i));

      v_bbcx_(iblock, j, i) = tmp * (v_up_(iblock, kmb-1, j, i) * cag_ 
          + v_snlat_(iblock, j, i) * v_vp_(iblock, kmb-1, j, i) * sag_);

      v_bbcy_(iblock, j, i) = tmp * (v_vp_(iblock, kmb-1, j, i) * cag_ 
          - v_snlat_(iblock, j, i) * v_up_(iblock, kmb-1, j, i) * sag_);
    }
#endif // CANUTO
    return ;
  }
 private:
  const double c0f_ = CppPconstMod::c0f;
  const double cag_ = CppPconstMod::cag;
  const double sag_ = CppPconstMod::sag;
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble3D v_bbcx_  = *p_v_bbcx;
  const ViewDouble3D v_bbcy_  = *p_v_bbcy;
  const ViewDouble3D v_snlat_ = *p_v_snlat;
  const ViewDouble4D v_up_    = *p_v_up;
  const ViewDouble4D v_vp_    = *p_v_vp;
};

// !---------------------------------------------------------------------
// !     ADVECTION + DIFFUSION + CORIOLIS
// !---------------------------------------------------------------------
class FunctorBclinc3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {   
    const int iblock = 0;
    v_dlu_(iblock, k, j, i) -= v_fcor_(iblock, j, i) * v_vp_(iblock, k, j, i);
    v_dlv_(iblock, k, j, i) += v_fcor_(iblock, j, i) * v_up_(iblock, k, j, i);
    return ;
  }
 private:
  const ViewDouble3D v_fcor_ = *p_v_fcor;
  const ViewDouble4D v_up_   = *p_v_up;
  const ViewDouble4D v_vp_   = *p_v_vp;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
};

// !---------------------------------------------------------------------
// !     PRESSURE GRADIENT FORCES
// !---------------------------------------------------------------------
class FunctorBclinc4 {
 public:
  FunctorBclinc4 (const double &aa) : aa_(aa) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {   
    const int iblock = 0;

    v_h0bf_(iblock, j, i) = v_h0bf_(iblock, j, i) * onbb_;

    v_work_(iblock, j, i) = aa_ * v_h0bf_(iblock, j, i) 
                  + (1.0 - aa_) * v_h0bl_(iblock, j, i);

    double wkk1, wkk2;
    wkk1 = (v_psa_(iblock, j, i) * od0_ + v_work_(iblock, j, i) * G) 
        * v_vit_(iblock, 0, j, i);
    for (int k = 0; k < KM; ++k) {
      wkk2 = wkk1 - v_gg_(iblock, k, j, i) 
          * v_dzp_(k) * v_vit_(iblock, k, j, i);
      v_wka_(iblock, k, j, i) = P5 * (wkk1 + wkk2);
      wkk1 = wkk2;
    }
    return ;
  }
 private:
  const double aa_;
  const double od0_  = CppPconstMod::od0;
  const double onbb_ = CppPconstMod::onbb;

  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_psa_  = *p_v_psa;
  const ViewDouble3D v_h0bf_ = *p_v_h0bf;
  const ViewDouble3D v_h0bl_ = *p_v_h0bl;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_gg_   = *p_v_gg;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka_  = *p_v_wka;
};
class FunctorBclinc5 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const  {   
    const int iblock = 0;
    double gradx, grady;
    grad(iblock, k, j, i, gradx, grady, v_wka_);
    v_dlu_(iblock, k, j, i) -= gradx;
    v_dlv_(iblock, k, j, i) -= grady;
    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &iblock, const int &k, 
      const int &j, const int &i, double &gradx, double &grady,
          const ViewDouble4D &v_f) const {
    const int bid = 0;
    gradx = 0.0;
    grady = 0.0;
    if (i >= 1 && j < (NY_BLOCK - 1)) {
      if (k <= (v_kmu_(bid, j, i) - 1)) {
        gradx = v_dxyur_(bid, j, i, 0) * P5 
            * (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j  , i-1) 
             - v_f(iblock, k, j+1, i-1) + v_f(iblock, k, j  , i  ));
      
        grady = v_dxyur_(bid, j, i, 1) * P5 
            * (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j  , i-1) 
             + v_f(iblock, k, j+1, i-1) - v_f(iblock, k, j  , i  ));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
};

class FunctorBclinc6 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = (1.0 + v_ohbt_(iblock, j, i) * v_zkt_(k)) 
        * v_work_(iblock, j, i) * v_vit_(iblock, k, j, i);
    return ;
  }
 private:
  const ViewDouble1D v_zkt_  = *p_v_zkt;
  const ViewDouble3D v_ohbt_ = *p_v_ohbt;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

class FunctorBclinc7 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    double gradx(0.0), grady(0.0);
    grad(iblock, k, j, i, gradx, grady, v_wka_);
    const double ggu = P25 * (
        v_gg_(iblock, k, j  , i) + v_gg_(iblock, k, j  , i-1) 
      + v_gg_(iblock, k, j+1, i) + v_gg_(iblock, k, j+1, i-1));
    
    v_dlu_(iblock, k, j, i) += ggu * gradx;
    v_dlv_(iblock, k, j, i) += ggu * grady;

    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &iblock, const int &k, 
      const int &j, const int &i, double &gradx, double &grady,
          const ViewDouble4D &v_f) const {
    const int bid = 0;
    gradx = 0.0;
    grady = 0.0;
    if (i >=1 && j < (NY_BLOCK - 1)) {
      if (k <= (v_kmu_(bid, j, i) - 1)) {
        gradx = v_dxyur_(bid, j, i, 0) * P5 
            * (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j  , i-1) 
             - v_f(iblock, k, j+1, i-1) + v_f(iblock, k, j  , i  ));
      
        grady = v_dxyur_(bid, j, i, 1) * P5 
            * (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j  , i-1) 
             + v_f(iblock, k, j+1, i-1) - v_f(iblock, k, j  , i  ));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble4D v_wka_   = *p_v_wka;
  const ViewDouble4D v_gg_    = *p_v_gg;
  const ViewDouble4D v_dlu_   = *p_v_dlu;
  const ViewDouble4D v_dlv_   = *p_v_dlv;
  const ViewDouble4D v_viv_   = *p_v_viv;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
};

// !---------------------------------------------------------------------
// !     CORIOLIS ADJUSTMENT
// !---------------------------------------------------------------------
class FunctorBclinc8 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    if (isc_ == 0) {
      const double wk1 = v_epea_(iblock, j, i) * v_dlv_(iblock, k, j, i) 
                       + v_epeb_(iblock, j, i) * v_dlu_(iblock, k, j, i);
      const double wk2 = v_epea_(iblock, j, i) * v_dlu_(iblock, k, j, i) 
                       - v_epeb_(iblock, j, i) * v_dlv_(iblock, k, j, i); 
    
      v_dlv_(iblock, k, j, i) = wk1 * v_viv_(iblock, k, j, i);
      v_dlu_(iblock, k, j, i) = wk2 * v_viv_(iblock, k, j, i); 
    } else {
      const double wk1 = v_epla_(iblock, j, i) * v_dlv_(iblock, k, j, i) 
                       + v_eplb_(iblock, j, i) * v_dlu_(iblock, k, j, i);
      const double wk2 = v_epla_(iblock, j, i) * v_dlu_(iblock, k, j, i) 
                       - v_eplb_(iblock, j, i) * v_dlv_(iblock, k, j, i); 
    
      v_dlv_(iblock, k, j, i) = wk1 * v_viv_(iblock, k, j, i);
      v_dlu_(iblock, k, j, i) = wk2 * v_viv_(iblock, k, j, i); 
    }
    return ;
  }
 private:
  const int isc_ = CppPconstMod::isc;
  const ViewDouble3D v_epea_ = *p_v_epea;
  const ViewDouble3D v_epeb_ = *p_v_epeb;
  const ViewDouble3D v_epla_ = *p_v_epla;
  const ViewDouble3D v_eplb_ = *p_v_eplb;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_vit_  = *p_v_vit;
};

// !---------------------------------------------------------------------
// !     PREDICTING VC & UC
// !---------------------------------------------------------------------
class FunctorBclinc9 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_u_(iblock, k, j, i) = v_up_(iblock, k, j, i) 
        + v_dlu_(iblock, k, j, i) * dtc_;
    v_v_(iblock, k, j, i) = v_vp_(iblock, k, j, i) 
        + v_dlv_(iblock, k, j, i) * dtc_;
    return ;
  }
 private:
  const double dtc_ = CppPconstMod::dtc;

  const ViewDouble4D v_u_   = *p_v_u;
  const ViewDouble4D v_v_   = *p_v_v;
  const ViewDouble4D v_up_  = *p_v_up;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_dlu_ = *p_v_dlu;
  const ViewDouble4D v_dlv_ = *p_v_dlv;
};

class FunctorBclinc10 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
#if (defined CANUTO) || (defined CANUTO2010)
    const double aidif = 0.5;
#else  // CANUTO
    const double aidif = 0.0; 
#endif // CANUTO
    invtriu(j, i, v_u_, v_sbcx_, v_bbcx_, v_akmu_, aidif, dtc_);
    invtriu(j, i, v_v_, v_sbcy_, v_bbcy_, v_akmu_, aidif, dtc_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void invtriu (const int &j, const int &i,
      const ViewDouble4D &v_wk,
      const ViewDouble3D &v_topbc,
      const ViewDouble3D &v_bombc,
      const ViewDouble4D &v_dcb,
      const double &aidif, const double &c2dtc) const {
    int k;
    const int iblock = 0;

    const double c2dtc_times_aidif = c2dtc * aidif;

    if (v_kmu_(iblock, j, i) > 0) {
      const int kz = v_kmu_(iblock, j, i);
      for (k = 1; k < kz; ++k) {
        v_a8_(k, j, i) = v_dcb(iblock, k-1, j, i) * v_odzt_(k) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      }
      for (k = 1; k < kz - 1 ; ++k) { 
        v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k + 1) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_b8_(k, j, i) = 1.0 +  v_a8_(k, j, i) + v_c8_(k, j, i);
        v_e8_(k, j, i) = 0.0;
        v_f8_(k, j, i) = 0.0;
      }
      //B.C. AT TOP
      k = 0;
      v_a8_(k, j, i) = v_odzp_(k) * c2dtc_times_aidif;
      v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k+1) * v_odzp_(k) 
          * c2dtc_times_aidif;
      v_b8_(k, j, i) = 1.0 + v_c8_(k, j, i);
      v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      v_e8_(k, j, i) = 0.0;
      v_f8_(k, j, i) = 0.0;
      //B.C. AT BOTTOM
      v_b8_(kz-1, j, i) = 1.0 + v_a8_(kz-1, j, i);
      v_c8_(kz-1, j, i) = v_odzp_(kz-1) * c2dtc_times_aidif;
      v_e8_(kz, j, i) = 0.0;
      v_f8_(kz, j, i) = 0.0;
      v_d8_(kz-1, j, i) = v_wk(iblock, kz-1, j, i) 
          - v_bombc(iblock, j, i) * v_odzp_(kz-1) * c2dtc_times_aidif;
      //NOW INVERT
      for (k = kz - 1; k >= 0; --k) {
        const double g0 = 1.0 / (v_b8_(k, j, i) - v_c8_(k, j, i) * v_e8_(k+1, j, i));
        v_e8_(k, j, i) = v_a8_(k, j, i) * g0;
        v_f8_(k, j, i) = (v_d8_(k, j, i) + v_c8_(k, j, i) * v_f8_(k+1, j, i)) * g0;
      }
     //B.C. AT SURFACE 
      k = 0;
      double wk = (v_e8_(k, j, i) * v_topbc(iblock, j, i) + v_f8_(k, j, i)) 
          * v_viv_(iblock, k, j, i); 
      v_wk(iblock, k, j, i) = wk;
      for (k = 1; k < kz ; ++k) {
        wk = (v_e8_(k, j, i) * wk + v_f8_(k, j, i)) * v_viv_(iblock, k, j, i);
        v_wk(iblock, k, j, i) = wk;
      } 
    }
    return ;
  }
 private:
  const double dtc_ = CppPconstMod::dtc;

  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_a8_   = *p_v_wp3;
  const ViewDouble3D v_b8_   = *p_v_wp4;
  const ViewDouble3D v_c8_   = *p_v_wp5;
  const ViewDouble3D v_d8_   = *p_v_wp6;
  const ViewDouble3D v_e8_   = *p_v_wk1;
  const ViewDouble3D v_f8_   = *p_v_wk2;
  const ViewDouble3D v_sbcx_ = *p_v_sbcx;
  const ViewDouble3D v_bbcx_ = *p_v_bbcx;
  const ViewDouble3D v_sbcy_ = *p_v_sbcy;
  const ViewDouble3D v_bbcy_ = *p_v_bbcy;
  const ViewDouble4D v_u_    = *p_v_u;
  const ViewDouble4D v_v_    = *p_v_v;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
};

// !---------------------------------------------------------------------
// !@@@  INTERACTION BETWEEN BAROTROPIC AND BAROCLINIC MODES
// !---------------------------------------------------------------------
class FunctorBclinc12 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    // work_1[0][j][i]
    vinteg_1 (j, i, v_u_, v_work_1_);
    // work_1[1][j][i]
    vinteg_2 (j, i, v_v_, v_work_1_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void vinteg_1 (const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble4D &v_wk2) 
          const {
    const int iblock = 0;
    double wk2 = C0;
    for (int k = 0; k < KM; ++k) {
      wk2 += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(iblock, k, j, i) *v_viv_(iblock, k, j, i);
    }
    v_wk2(iblock, 0, j, i) = wk2;
    return ;
  }
  KOKKOS_INLINE_FUNCTION void vinteg_2 (const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble4D &v_wk2) 
          const {
    const int iblock = 0;
    double wk2 = C0;
    for (int k = 0; k < KM; ++k) {
      wk2 += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(iblock, k, j, i) *v_viv_(iblock, k, j, i);
    }
    v_wk2(iblock, 1, j, i) = wk2;
    return ;
  }
 private:
  const ViewDouble1D v_dzp_    = *p_v_dzp;
  const ViewDouble3D v_ohbu_   = *p_v_ohbu;
  const ViewDouble4D v_u_      = *p_v_u;
  const ViewDouble4D v_v_      = *p_v_v;
  const ViewDouble4D v_viv_    = *p_v_viv;
  const ViewDouble4D v_work_1_ = *p_v_work_1;
};
class FunctorBclinc13 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_u_(iblock, k, j, i) = (v_u_(iblock, k, j, i) 
        - v_work_1_(iblock, 0, j, i) + v_ub_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);

    v_v_(iblock, k, j, i) = (v_v_(iblock, k, j, i) 
        - v_work_1_(iblock, 1, j, i) + v_vb_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);
    // if (j >= (JST-1) && j < JET) {
    v_utf_(iblock, k, j, i) += v_u_(iblock, k, j, i);
    v_vtf_(iblock, k, j, i) += v_v_(iblock, k, j, i);
    // }
    return ;
  }
 private:
  const ViewDouble3D v_ub_     = *p_v_ub;
  const ViewDouble3D v_vb_     = *p_v_vb;
  const ViewDouble4D v_u_      = *p_v_u;
  const ViewDouble4D v_v_      = *p_v_v;
  const ViewDouble4D v_utf_    = *p_v_utf;
  const ViewDouble4D v_vtf_    = *p_v_vtf;
  const ViewDouble4D v_viv_    = *p_v_viv;
  const ViewDouble4D v_work_1_ = *p_v_work_1;
};

class FunctorBclinc14 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = v_vp_(iblock, k, j, i)
        + v_dlv_(iblock, k, j, i) * dtc2_;
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_dlv_ = *p_v_dlv;
  const ViewDouble4D v_wka_ = *p_v_wka;
};

class FunctorBclinc15 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
#if (defined CANUTO) || (defined CANUTO2010)
    const double aidif = 0.5;
#else  // CANUTO
    const double aidif = 0.0; 
#endif // CANUTO
    invtriu(j, i, v_wka_, v_sbcy_, v_bbcy_, v_akmu_, aidif, dtc2_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void invtriu (const int &j, const int &i,
      const ViewDouble4D &v_wk,
      const ViewDouble3D &v_topbc,
      const ViewDouble3D &v_bombc,
      const ViewDouble4D &v_dcb,
      const double &aidif, const double &c2dtc) const {
    int k;
    const int iblock = 0;

    const double c2dtc_times_aidif = c2dtc * aidif;

    if (v_kmu_(iblock, j, i) > 0) {
      const int kz = v_kmu_(iblock, j, i);
      for (k = 1; k < kz; ++k) {
        v_a8_(k, j, i) = v_dcb(iblock, k-1, j, i) * v_odzt_(k) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      }
      for (k = 1; k < kz - 1 ; ++k) { 
        v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k + 1) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_b8_(k, j, i) = 1.0 +  v_a8_(k, j, i) + v_c8_(k, j, i);
        v_e8_(k, j, i) = 0.0;
        v_f8_(k, j, i) = 0.0;
      }
      //B.C. AT TOP
      k = 0;
      v_a8_(k, j, i) = v_odzp_(k) * c2dtc_times_aidif;
      v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k+1) * v_odzp_(k) 
          * c2dtc_times_aidif;
      v_b8_(k, j, i) = 1.0 + v_c8_(k, j, i);
      v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      v_e8_(k, j, i) = 0.0;
      v_f8_(k, j, i) = 0.0;
      //B.C. AT BOTTOM
      v_b8_(kz-1, j, i) = 1.0 + v_a8_(kz-1, j, i);
      v_c8_(kz-1, j, i) = v_odzp_(kz-1) * c2dtc_times_aidif;
      v_e8_(kz, j, i) = 0.0;
      v_f8_(kz, j, i) = 0.0;
      v_d8_(kz-1, j, i) = v_wk(iblock, kz-1, j, i) 
          - v_bombc(iblock, j, i) * v_odzp_(kz-1) * c2dtc_times_aidif;
      //NOW INVERT
      for (k = kz - 1; k >= 0; --k) {
        const double g0 = 1.0 / (v_b8_(k, j, i) - v_c8_(k, j, i) * v_e8_(k+1, j, i));
        v_e8_(k, j, i) = v_a8_(k, j, i) * g0;
        v_f8_(k, j, i) = (v_d8_(k, j, i) + v_c8_(k, j, i) * v_f8_(k+1, j, i)) * g0;
      }
     //B.C. AT SURFACE 
      k = 0;
      double wk = (v_e8_(k, j, i) * v_topbc(iblock, j, i) + v_f8_(k, j, i)) 
          * v_viv_(iblock, k, j, i); 
      v_wk(iblock, k, j, i) = wk;
      for (k = 1; k < kz ; ++k) {
        wk = (v_e8_(k, j, i) * wk + v_f8_(k, j, i)) * v_viv_(iblock, k, j, i);
        v_wk(iblock, k, j, i) = wk;
      } 
    }
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_a8_   = *p_v_wp3;
  const ViewDouble3D v_b8_   = *p_v_wp4;
  const ViewDouble3D v_c8_   = *p_v_wp5;
  const ViewDouble3D v_d8_   = *p_v_wp6;
  const ViewDouble3D v_e8_   = *p_v_wk1;
  const ViewDouble3D v_f8_   = *p_v_wk2;
  const ViewDouble3D v_sbcy_ = *p_v_sbcy;
  const ViewDouble3D v_bbcy_ = *p_v_bbcy;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
};

class FunctorBclinc16 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    vinteg(j, i, v_wka_, v_work_);
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
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_wka_  = *p_v_wka;
};
// !---------------------------------------------------------------------
// !     FILTER FORCING AT HIGT LATITUDES
// !---------------------------------------------------------------------
class FunctorBclinc17 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = (v_wka_(iblock, k, j, i) 
        - v_work_(iblock, j, i) + v_vb_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);
    v_vp_(iblock, k, j, i) = afc2_ * v_v_(iblock, k, j, i) 
        + afc1_ * (v_vp_(iblock, k, j, i) 
            + v_wka_(iblock, k, j, i));
    
    v_v_(iblock, k, j, i) = v_wka_(iblock, k, j, i);

    if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      v_wka_(iblock, k, j, i) = v_up_(iblock, k, j, i) 
          + v_dlu_(iblock, k, j, i) * dtc2_;
    }

    return ;
  }
 private:
  const double afc1_ = CppPconstMod::afc1;
  const double afc2_ = CppPconstMod::afc2;
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewDouble3D v_vb_   = *p_v_vb;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_v_    = *p_v_v;
  const ViewDouble4D v_up_   = *p_v_up;
  const ViewDouble4D v_vp_   = *p_v_vp;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_wka_  = *p_v_wka;
};
class FunctorBclinc18 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
#if (defined CANUTO) || (defined CANUTO2010)
    const double aidif = 0.5;
#else  // CANUTO
    const double aidif = 0.0; 
#endif // CANUTO
    invtriu(j, i, v_wka_, v_sbcx_, v_bbcx_, v_akmu_, aidif, dtc2_);
    return ;
  }
KOKKOS_INLINE_FUNCTION void invtriu (const int &j, const int &i,
      const ViewDouble4D &v_wk,
      const ViewDouble3D &v_topbc,
      const ViewDouble3D &v_bombc,
      const ViewDouble4D &v_dcb,
      const double &aidif, const double &c2dtc) const {
    int k;
    const int iblock = 0;

    const double c2dtc_times_aidif = c2dtc * aidif;

    if (v_kmu_(iblock, j, i) > 0) {
      const int kz = v_kmu_(iblock, j, i);
      for (k = 1; k < kz; ++k) {
        v_a8_(k, j, i) = v_dcb(iblock, k-1, j, i) * v_odzt_(k) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      }
      for (k = 1; k < kz - 1 ; ++k) { 
        v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k + 1) 
            * v_odzp_(k) * c2dtc_times_aidif;
        v_b8_(k, j, i) = 1.0 +  v_a8_(k, j, i) + v_c8_(k, j, i);
        v_e8_(k, j, i) = 0.0;
        v_f8_(k, j, i) = 0.0;
      }
      //B.C. AT TOP
      k = 0;
      v_a8_(k, j, i) = v_odzp_(k) * c2dtc_times_aidif;
      v_c8_(k, j, i) = v_dcb(iblock, k, j, i) * v_odzt_(k+1) * v_odzp_(k) 
          * c2dtc_times_aidif;
      v_b8_(k, j, i) = 1.0 + v_c8_(k, j, i);
      v_d8_(k, j, i) = v_wk(iblock, k, j, i);
      v_e8_(k, j, i) = 0.0;
      v_f8_(k, j, i) = 0.0;
      //B.C. AT BOTTOM
      v_b8_(kz-1, j, i) = 1.0 + v_a8_(kz-1, j, i);
      v_c8_(kz-1, j, i) = v_odzp_(kz-1) * c2dtc_times_aidif;
      v_e8_(kz, j, i) = 0.0;
      v_f8_(kz, j, i) = 0.0;
      v_d8_(kz-1, j, i) = v_wk(iblock, kz-1, j, i) 
          - v_bombc(iblock, j, i) * v_odzp_(kz-1) * c2dtc_times_aidif;
      //NOW INVERT
      for (k = kz - 1; k >= 0; --k) {
        const double g0 = 1.0 / (v_b8_(k, j, i) - v_c8_(k, j, i) * v_e8_(k+1, j, i));
        v_e8_(k, j, i) = v_a8_(k, j, i) * g0;
        v_f8_(k, j, i) = (v_d8_(k, j, i) + v_c8_(k, j, i) * v_f8_(k+1, j, i)) * g0;
      }
     //B.C. AT SURFACE 
      k = 0;
      double wk = (v_e8_(k, j, i) * v_topbc(iblock, j, i) + v_f8_(k, j, i)) 
          * v_viv_(iblock, k, j, i); 
      v_wk(iblock, k, j, i) = wk;
      for (k = 1; k < kz ; ++k) {
        wk = (v_e8_(k, j, i) * wk + v_f8_(k, j, i)) * v_viv_(iblock, k, j, i);
        v_wk(iblock, k, j, i) = wk;
      } 
    }
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble3D v_a8_   = *p_v_wp3;
  const ViewDouble3D v_b8_   = *p_v_wp4;
  const ViewDouble3D v_c8_   = *p_v_wp5;
  const ViewDouble3D v_d8_   = *p_v_wp6;
  const ViewDouble3D v_e8_   = *p_v_wk1;
  const ViewDouble3D v_f8_   = *p_v_wk2;
  const ViewDouble3D v_sbcx_ = *p_v_sbcx;
  const ViewDouble3D v_bbcx_ = *p_v_bbcx;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
};
class FunctorBclinc19 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i)  const {
    vinteg(j, i, v_wka_, v_work_);
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
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble3D v_ohbu_ = *p_v_ohbu;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

// !---------------------------------------------------------------------
// !     FILTER FORCING AT HIGT LATITUDES
// !---------------------------------------------------------------------
class FunctorBclinc20 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i)  const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = (v_wka_(iblock, k, j, i) 
        - v_work_(iblock, j, i) + v_ub_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);

    v_up_(iblock, k, j, i) = afc2_ * v_u_(iblock, k, j, i) 
      + afc1_ * (v_up_(iblock, k, j, i) + v_wka_(iblock, k, j, i));

    v_u_(iblock, k, j, i) = v_wka_(iblock, k, j, i);

    // if (j >= (JST-1) && j < JET) {
    v_utf_(iblock, k, j, i) += v_u_(iblock, k, j, i);
    v_vtf_(iblock, k, j, i) += v_v_(iblock, k, j, i);
    // }
    return ;
  }
 private:
  const double afc1_ = CppPconstMod::afc1;
  const double afc2_ = CppPconstMod::afc2;
  const ViewDouble3D v_ub_   = *p_v_ub;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_u_    = *p_v_u;
  const ViewDouble4D v_v_    = *p_v_v;
  const ViewDouble4D v_up_   = *p_v_up;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_utf_  = *p_v_utf;
  const ViewDouble4D v_vtf_  = *p_v_vtf;
};
// merged version

class FunctorBclincMer14 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_work_merged_(0, k, j, i) = v_up_(iblock, k, j, i)
        + v_dlu_(iblock, k, j, i) * dtc2_;
    v_work_merged_(1, k, j, i) = v_vp_(iblock, k, j, i)
        + v_dlv_(iblock, k, j, i) * dtc2_;
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewDouble4D v_up_  = *p_v_up;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_dlu_ = *p_v_dlu;
  const ViewDouble4D v_dlv_ = *p_v_dlv;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;
};

class FunctorBclincMer15 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &j, const int &i) const {
#if (defined CANUTO) || (defined CANUTO2010)
    const double aidif = 0.5;
#else  // CANUTO
    const double aidif = 0.0; 
#endif // CANUTO
    invtriu_u (j, i, v_work_merged_, v_sbcx_, v_bbcx_, v_akmu_, aidif, dtc2_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void invtriu_u(const int &j, const int &i,
      const ViewDouble4D &v_wk,
      const ViewDouble3D &v_topbc,
      const ViewDouble3D &v_bombc,
      const ViewDouble4D &v_dcb,
      const double &aidif, const double &c2dtc) const {
    int k;
    const int idx_u  = 0;
    const int iblock = 0;
    double a8[KM], b8[KM], c8[KM], d8[KM];
    double e8[KM+1], f8[KM+1];
    const double c2dtc_times_aidif = c2dtc * aidif;

    if (v_kmu_(iblock, j, i) > 0) {
      const int kz = v_kmu_(iblock, j, i);
      for (k = 1; k < kz; ++k) {
        a8[k] = v_dcb(iblock, k-1, j, i) * v_odzt_(k) 
            * v_odzp_(k) * c2dtc_times_aidif;
        d8[k] = v_wk(idx_u, k, j, i);
      }
      for (k = 1; k < kz - 1 ; ++k) { 
        c8[k] = v_dcb(iblock, k, j, i) * v_odzt_(k + 1) 
            * v_odzp_(k) * c2dtc_times_aidif;
        b8[k] = 1.0 +  a8[k] + c8[k];
        e8[k] = 0.0;
        f8[k] = 0.0;
      }
      //B.C. AT TOP
      k = 0;
      a8[k] = v_odzp_(k) * c2dtc_times_aidif;
      c8[k] = v_dcb(iblock, k, j, i) * v_odzt_(k+1) * v_odzp_(k) 
          * c2dtc_times_aidif;
      b8[k] = 1.0 + c8[k];
      d8[k] = v_wk(idx_u, k, j, i);
      e8[k] = 0.0;
      f8[k] = 0.0;
      //B.C. AT BOTTOM
      b8[kz-1] = 1.0 + a8[kz-1];
      c8[kz-1] = v_odzp_[kz-1] * c2dtc_times_aidif;
      e8[kz] = 0.0;
      f8[kz] = 0.0;
      d8[kz-1] = v_wk(idx_u, kz-1, j, i) 
          - v_bombc(iblock, j, i) * v_odzp_(kz-1) * c2dtc_times_aidif;
      //NOW INVERT
      for (k = kz - 1; k >= 0; --k) {
        const double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
        e8[k] = a8[k] * g0;
        f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
      }
     //B.C. AT SURFACE 
      k = 0;
      double wk = (e8[k] * v_topbc(iblock, j, i) + f8[k]) 
          * v_viv_(iblock, k, j, i); 
      v_wk(idx_u, k, j, i) = wk;
      for (k = 1; k < kz ; ++k) {
        wk = (e8[k] * wk + f8[k]) * v_viv_(iblock, k, j, i);
        v_wk(idx_u, k, j, i) = wk;
      } 
    }
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_sbcx_ = *p_v_sbcx;
  const ViewDouble3D v_bbcx_ = *p_v_bbcx;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;
};

class FunctorBclincMer16 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &j, const int &i) const {
#if (defined CANUTO) || (defined CANUTO2010)
    const double aidif = 0.5;
#else  // CANUTO
    const double aidif = 0.0; 
#endif // CANUTO
    invtriu_v (j, i, v_work_merged_, v_sbcy_, v_bbcy_, v_akmu_, aidif, dtc2_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void invtriu_v (const int &j, const int &i,
      const ViewDouble4D &v_wk,
      const ViewDouble3D &v_topbc,
      const ViewDouble3D &v_bombc,
      const ViewDouble4D &v_dcb,
      const double &aidif, const double &c2dtc) const {
    int k;
    const int idx_v  = 1;
    const int iblock = 0;
    double a8[KM], b8[KM], c8[KM], d8[KM];
    double e8[KM+1], f8[KM+1];
    const double c2dtc_times_aidif = c2dtc * aidif;

    if (v_kmu_(iblock, j, i) > 0) {
      const int kz = v_kmu_(iblock, j, i);
      for (k = 1; k < kz; ++k) {
        a8[k] = v_dcb(iblock, k-1, j, i) * v_odzt_(k) 
            * v_odzp_(k) * c2dtc_times_aidif;
        d8[k] = v_wk(idx_v, k, j, i);
      }
      for (k = 1; k < kz - 1 ; ++k) { 
        c8[k] = v_dcb(iblock, k, j, i) * v_odzt_(k + 1) 
            * v_odzp_(k) * c2dtc_times_aidif;
        b8[k] = 1.0 +  a8[k] + c8[k];
        e8[k] = 0.0;
        f8[k] = 0.0;
      }
      //B.C. AT TOP
      k = 0;
      a8[k] = v_odzp_(k) * c2dtc_times_aidif;
      c8[k] = v_dcb(iblock, k, j, i) * v_odzt_(k+1) * v_odzp_(k) 
          * c2dtc_times_aidif;
      b8[k] = 1.0 + c8[k];
      d8[k] = v_wk(idx_v, k, j, i);
      e8[k] = 0.0;
      f8[k] = 0.0;
      //B.C. AT BOTTOM
      b8[kz-1] = 1.0 + a8[kz-1];
      c8[kz-1] = v_odzp_[kz-1] * c2dtc_times_aidif;
      e8[kz] = 0.0;
      f8[kz] = 0.0;
      d8[kz-1] = v_wk(idx_v, kz-1, j, i) 
          - v_bombc(iblock, j, i) * v_odzp_(kz-1) * c2dtc_times_aidif;
      //NOW INVERT
      for (k = kz - 1; k >= 0; --k) {
        const double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
        e8[k] = a8[k] * g0;
        f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
      }
     //B.C. AT SURFACE 
      k = 0;
      double wk = (e8[k] * v_topbc(iblock, j, i) + f8[k]) 
          * v_viv_(iblock, k, j, i); 
      v_wk(idx_v, k, j, i) = wk;
      for (k = 1; k < kz ; ++k) {
        wk = (e8[k] * wk + f8[k]) * v_viv_(iblock, k, j, i);
        v_wk(idx_v, k, j, i) = wk;
      } 
    }
    return ;
  }
 private:
  const double dtc2_ = CppPconstMod::dtc2;
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_sbcy_ = *p_v_sbcy;
  const ViewDouble3D v_bbcy_ = *p_v_bbcy;
  const ViewDouble4D v_viv_  = *p_v_viv;
  const ViewDouble4D v_akmu_ = *p_v_akmu;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;
};

class FunctorBclincMer17 {
 public:
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &j, const int &i) const {
    vinteg_u (j, i, v_work_merged_, v_work_u_);
    vinteg_v (j, i, v_work_merged_, v_work_v_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void vinteg_u (const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble3D &v_wk2) 
          const {
    const int idx_u  = 0;
    const int iblock = 0;
    double wk2 = C0;
    for (int k = 0; k < KM; ++k) {
      wk2 += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(idx_u, k, j, i) * v_viv_(iblock, k, j, i);
    }
    v_wk2(iblock, j, i) = wk2;
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void vinteg_v (const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble3D &v_wk2) 
          const {
    const int idx_v  = 1;
    const int iblock = 0;
    double wk2 = C0;
    for (int k = 0; k < KM; ++k) {
      wk2 += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(idx_v, k, j, i) * v_viv_(iblock, k, j, i);
    }
    v_wk2(iblock, j, i) = wk2;
    return ;
  }
 private:
  const ViewDouble1D v_dzp_    = *p_v_dzp;
  const ViewDouble3D v_ohbu_   = *p_v_ohbu;
  const ViewDouble3D v_work_u_ = *p_v_work1;
  const ViewDouble3D v_work_v_ = *p_v_work2;
  const ViewDouble4D v_viv_    = *p_v_viv;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;
};

class FunctorBclincMer18 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i)  const {
    const int idx_u  = 0;
    const int iblock = 0;
    v_work_merged_(idx_u, k, j, i) = (v_work_merged_(idx_u, k, j, i) 
        - v_work_u_(iblock, j, i) + v_ub_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);

    // v_work_merged_(idx_u, k, j, i) = (v_work_merged_(idx_u, k, j, i) 
    //     - v_work_u_(j, i) + v_ub_(iblock, j, i)) 
    //         * v_viv_(iblock, k, j, i);
    // v_work_merged_(idx_u, k, j, i) = (v_wka_(iblock, k, j, i) 
    //     - v_work_(iblock, j, i) + v_ub_(iblock, j, i)) 
    //         * v_viv_(iblock, k, j, i);

    // v_up_(iblock, k, j, i) = afc2_ * v_u_(iblock, k, j, i) 
    //   + afc1_ * (v_up_(iblock, k, j, i) + v_wka_(iblock, k, j, i));

    // v_u_(iblock, k, j, i) = v_wka_(iblock, k, j, i);

    v_up_(iblock, k, j, i) = afc2_ * v_u_(iblock, k, j, i) 
      + afc1_ * (v_up_(iblock, k, j, i) + v_work_merged_(idx_u, k, j, i));

    v_u_(iblock, k, j, i) = v_work_merged_(idx_u, k, j, i);

    // if (j >= (JST-1) && j < JET) {
    v_utf_(iblock, k, j, i) += v_u_(iblock, k, j, i);
    // }
    return ;
  }
 private:
  const double afc1_ = CppPconstMod::afc1;
  const double afc2_ = CppPconstMod::afc2;
  const ViewDouble3D v_ub_      = *p_v_ub;
  const ViewDouble3D v_work_u_  = *p_v_work1;
  const ViewDouble4D v_u_       = *p_v_u;
  const ViewDouble4D v_up_      = *p_v_up;
  const ViewDouble4D v_viv_     = *p_v_viv;
  const ViewDouble4D v_utf_     = *p_v_utf;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;

  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_wka_ = *p_v_wka;
};

class FunctorBclincMer19 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i)  const {
    const int idx_v  = 1;
    const int iblock = 0;
    v_work_merged_(idx_v, k, j, i) = (v_work_merged_(idx_v, k, j, i) 
        - v_work_v_(iblock, j, i) + v_vb_(iblock, j, i)) 
            * v_viv_(iblock, k, j, i);

    v_vp_(iblock, k, j, i) = afc2_ * v_v_(iblock, k, j, i) 
      + afc1_ * (v_vp_(iblock, k, j, i) + v_work_merged_(idx_v, k, j, i));

    v_v_(iblock, k, j, i) = v_work_merged_(idx_v, k, j, i);

    // if (j >= (JST-1) && j < JET) {
    v_vtf_(iblock, k, j, i) += v_v_(iblock, k, j, i);
    // }
    return ;
  }
 private:
  const double afc1_ = CppPconstMod::afc1;
  const double afc2_ = CppPconstMod::afc2;
  const ViewDouble3D v_vb_      = *p_v_vb;
  const ViewDouble3D v_work_v_  = *p_v_work2;
  const ViewDouble4D v_v_       = *p_v_v;
  const ViewDouble4D v_vp_      = *p_v_vp;
  const ViewDouble4D v_viv_     = *p_v_viv;
  const ViewDouble4D v_vtf_     = *p_v_vtf;
  const ViewDouble4D v_work_merged_ = *p_v_work_merged;
};



KOKKOS_REGISTER_FOR_2D(FunctorBclinc1,  FunctorBclinc1)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc2,  FunctorBclinc2)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc3,  FunctorBclinc3)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc4,  FunctorBclinc4)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc5,  FunctorBclinc5)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc6,  FunctorBclinc6)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc7,  FunctorBclinc7)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc8,  FunctorBclinc8)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc9,  FunctorBclinc9)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc10, FunctorBclinc10)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc12, FunctorBclinc12)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc13, FunctorBclinc13)

// Original
KOKKOS_REGISTER_FOR_3D(FunctorBclinc14, FunctorBclinc14)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc15, FunctorBclinc15)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc16, FunctorBclinc16)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc17, FunctorBclinc17)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc18, FunctorBclinc18)
KOKKOS_REGISTER_FOR_2D(FunctorBclinc19, FunctorBclinc19)
KOKKOS_REGISTER_FOR_3D(FunctorBclinc20, FunctorBclinc20)

// Merged
KOKKOS_REGISTER_FOR_3D(FunctorBclincMer14, FunctorBclincMer14)
KOKKOS_REGISTER_FOR_2D(FunctorBclincMer15, FunctorBclincMer15)
KOKKOS_REGISTER_FOR_2D(FunctorBclincMer16, FunctorBclincMer16)
KOKKOS_REGISTER_FOR_2D(FunctorBclincMer17, FunctorBclincMer17)
KOKKOS_REGISTER_FOR_3D(FunctorBclincMer18, FunctorBclincMer18)
KOKKOS_REGISTER_FOR_3D(FunctorBclincMer19, FunctorBclincMer19)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BCLINC_HPP_
