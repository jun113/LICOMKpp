#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYT_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYT_HPP_

#include "../head/def-undef.h"

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_grid.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_pmix_mod.h"
#include "../head/kokkos_tracer_mod.h"
#include "../head/kokkos_tmp_var.h"
#include "../head/kokkos_work_mod.h"

#include "../head/kokkos_config.hpp"

#include "../head/fortran_extern_functions.h"

#include "Kokkos_Core.hpp"

#include "math.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KMM1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::NTRA;
using CppParamMod::JST;
using CppParamMod::JET;

using CppConstantMod::P5;
using CppConstantMod::P25;
using CppConstantMod::C0;
using CppConstantMod::G;


using KokkosDynMod::p_v_dlu;
using KokkosDynMod::p_v_dlv;
using KokkosDynMod::p_v_gg;
using KokkosDynMod::p_v_h0;
using KokkosDynMod::p_v_h0l;
using KokkosDynMod::p_v_h0f;
using KokkosDynMod::p_v_u;
using KokkosDynMod::p_v_v;
using KokkosDynMod::p_v_utl;
using KokkosDynMod::p_v_vtl;
using KokkosDynMod::p_v_utf;
using KokkosDynMod::p_v_vtf;

using KokkosForcMod::p_v_psa;
using KokkosForcMod::p_v_swv;
using KokkosForcMod::p_v_nswv;
using KokkosForcMod::p_v_buoytur;
using KokkosForcMod::p_v_buoysol;

// using KokkosGrid::p_v_au0;
// using KokkosGrid::p_v_aus;
// using KokkosGrid::p_v_auw;
// using KokkosGrid::p_v_ausw;
using KokkosGrid::p_v_kmu;
using KokkosGrid::p_v_dxur;
using KokkosGrid::p_v_dyur;
using KokkosGrid::p_v_dxyur;

using KokkosPconstMod::p_v_c;
using KokkosPconstMod::p_v_dzp;
using KokkosPconstMod::p_v_hbx;
using KokkosPconstMod::p_v_hby;
using KokkosPconstMod::p_v_to;
using KokkosPconstMod::p_v_po;
using KokkosPconstMod::p_v_so;
using KokkosPconstMod::p_v_akt;
using KokkosPconstMod::p_v_viv;
using KokkosPconstMod::p_v_vit;
using KokkosPconstMod::p_v_odzt;
using KokkosPconstMod::p_v_ohbu;
using KokkosPconstMod::p_v_ohbt;
using KokkosPconstMod::p_v_zkt;
#ifdef CANUTOMIXOUT
using KokkosPconstMod::p_v_alpha_canuto;
using KokkosPconstMod::p_v_beta_canuto;
#endif // CANUTOMIXOUT

// using KokkosPmixMod::p_v_ric;
using KokkosPmixMod::p_v_rit;
using KokkosPmixMod::p_v_rict;
using KokkosPmixMod::p_v_ricdt;
using KokkosPmixMod::p_v_ricdttms;
// using KokkosPmixMod::p_v_rict_replace;

using KokkosWorkMod::p_v_pax;
using KokkosWorkMod::p_v_pay;
using KokkosWorkMod::p_v_pxb;
using KokkosWorkMod::p_v_pyb;
using KokkosWorkMod::p_v_whx;
using KokkosWorkMod::p_v_why;
using KokkosWorkMod::p_v_wgp;
using KokkosWorkMod::p_v_work;

using KokkosTracerMod::p_v_ax;
using KokkosTracerMod::p_v_ay;
using KokkosTracerMod::p_v_az;
using KokkosTracerMod::p_v_at;
using KokkosTracerMod::p_v_atb;
using KokkosTracerMod::p_v_pdensity;

using KokkosTmpVar::p_v_pp;
using KokkosTmpVar::p_v_ppa;
using KokkosTmpVar::p_v_ppb;
using KokkosTmpVar::p_v_ppc;
using KokkosTmpVar::p_v_work1;
using KokkosTmpVar::p_v_work2;
using KokkosTmpVar::p_v_alpha;
using KokkosTmpVar::p_v_beta;


class FunctorReadyt1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    for (int n = 0; n < NTRA; ++n) {
      v_akt_(iblock, n, k, j, i) = 0.0;
    }

    v_utl_(iblock, k, j, i) = v_utf_(iblock, k, j, i);
    v_vtl_(iblock, k, j, i) = v_vtf_(iblock, k, j, i);
    v_utf_(iblock, k, j, i) = v_u_(iblock, k, j, i);
    v_vtf_(iblock, k, j, i) = v_v_(iblock, k, j, i);

    return ;
  };

 private:
  const ViewDouble4D v_u_   = *p_v_u;
  const ViewDouble4D v_v_   = *p_v_v;
  const ViewDouble4D v_utl_ = *p_v_utl;
  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_utf_ = *p_v_utf;
  const ViewDouble4D v_vtf_ = *p_v_vtf;
  const ViewDouble5D v_akt_ = *p_v_akt;
};

class FunctorReadyt2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    if (k < KMM1) {
      const double tup = v_atb_(iblock, 0, k+1, j, i) - v_to_(k+1);
      const double sup = v_atb_(iblock, 1, k+1, j, i) - v_so_(k+1);
      const double tlo = v_atb_(iblock, 0, k+2, j, i) - v_to_(k+1);
      const double slo = v_atb_(iblock, 1, k+2, j, i) - v_so_(k+1);
    
      const double rhoup = dens(tup, sup, k+1);
      const double rholo = dens(tlo, slo, k+1);
      v_rict_(iblock, k, j, i) = v_vit_(iblock, k+1, j, i) 
          * od0_ * G * (rholo - rhoup) * v_odzt_(k+1);
    }
    // jst = 1, jet = ny_block
    // if (j >= (JST-1) && j < JET) {
      density(k, j, i);
    // }
    return ;
  };
  KOKKOS_INLINE_FUNCTION void density (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    if (v_vit_(iblock, k, j, i) > 0.0) {
      const double tq = v_at_(iblock, 0, k, j, i) - v_to_(k); 
      const double sq = v_at_(iblock, 1, k, j, i) - v_so_(k); 

      v_pdensity_(iblock, k, j, i) = 1.0e+3 + v_po_(k) +
          (v_c_(0, k) + (v_c_(3, k) + v_c_(6, k) * sq) *sq +
          (v_c_(2, k) +  v_c_(7, k) * sq + v_c_(5, k) * tq) * tq) * tq +
          (v_c_(1, k) + (v_c_(4, k) + v_c_(8, k) * sq) * sq) * sq;
    } else {
      v_pdensity_(iblock, k, j, i) = 0.0;
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION double dens (
      const double &tq, const double &sq, const int &kk) const {
    double dens;
    dens = (v_c_(0, kk) + (v_c_(3, kk) + v_c_(6, kk) * sq) * sq +
           (v_c_(2, kk) +  v_c_(7, kk) * sq + v_c_(5, kk) * tq) * tq) * tq +
           (v_c_(1, kk) + (v_c_(4, kk) + v_c_(8, kk) * sq) * sq) * sq;
    return dens;
  }

 private:
  const double od0_ = CppPconstMod::od0;
  const ViewDouble1D v_to_       = *p_v_to;
  const ViewDouble1D v_so_       = *p_v_so;
  const ViewDouble1D v_po_       = *p_v_po;
  const ViewDouble1D v_odzt_     = *p_v_odzt;
  const ViewDouble2D v_c_        = *p_v_c;
  const ViewDouble4D v_vit_      = *p_v_vit;
  const ViewDouble4D v_rict_     = *p_v_rict;
  const ViewDouble4D v_pdensity_ = *p_v_pdensity;
  const ViewDouble5D v_at_       = *p_v_at;
  const ViewDouble5D v_atb_      = *p_v_atb;
};

class FunctorReadyt3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_gg_(iblock, k, j, i) = - od0_ * G * v_pdensity_(iblock, k, j, i) 
      * v_vit_(iblock, k, j, i);

    if (k != 0) {
      // calculate pressure at T-grid
      v_ppb_(iblock, k, j, i) = v_vit_(iblock, k, j, i) *
          (v_at_(iblock, 0, k-1, j, i) -
          (v_at_(iblock, 0, k-1, j, i) - v_at_(iblock, 0, k, j, i)) *
           v_dzp_(k-1) / (v_dzp_(k-1) + v_dzp_(k)));

      v_ppc_(iblock, k, j, i) = v_vit_(iblock, k, j, i) *
          (v_at_(iblock, 1, k-1, j, i) -
          (v_at_(iblock, 1, k-1, j, i) - v_at_(iblock, 1, k, j, i)) *
           v_dzp_(k-1) / (v_dzp_(k-1) + v_dzp_(k)));
    } else {
      // calculate pressure at level 0 (Level 1 in Fortran)
      v_ppb_(iblock, 0, j, i) = v_at_(iblock, 0, 0, j, i) * 
          v_vit_(iblock, 0, j, i);
      v_ppc_(iblock, 0, j, i) = v_at_(iblock, 1, 0, j, i) * 
          v_vit_(iblock, 0, j, i);
    }
    return ;
  };

 private:
  const double od0_ = CppPconstMod::od0;
  const ViewDouble1D v_dzp_      = *p_v_dzp;
  const ViewDouble4D v_gg_       = *p_v_gg;
  const ViewDouble4D v_vit_      = *p_v_vit;
  const ViewDouble4D v_ppb_      = *p_v_ppb;
  const ViewDouble4D v_ppc_      = *p_v_ppc;
  const ViewDouble4D v_pdensity_ = *p_v_pdensity;
  const ViewDouble5D v_at_       = *p_v_at;
};

class FunctorReadyt4 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_h0l_(iblock, j, i) = v_h0f_(iblock, j, i);
    v_h0f_(iblock, j, i) = v_h0_(iblock, j, i);

    v_pp_(iblock, 0, j, i)  = v_gg_(iblock, 0, j, i) * 0.5 
        * v_dzp_(0) * v_vit_(iblock, 0, j, i);
    v_ppa_(iblock, 0, j, i) = v_psa_(iblock, j, i) 
        * v_vit_(iblock, 0, j, i);

    for (int k = 1; k < KM; ++k) {
      // calculate pressure at T-grid
      v_pp_(iblock, k, j, i) = v_vit_(iblock, k, j, i) *
          (v_pp_(iblock, k-1, j, i) + 0.5 * 
          (v_gg_(iblock, k  , j, i) * v_dzp_(k) +
           v_gg_(iblock, k-1, j, i) * v_dzp_(k-1)));

      v_ppa_(iblock, k, j, i) = v_vit_(iblock, k, j, i) *
          (v_ppa_(iblock, k-1, j, i) + v_gg_(iblock, k-1, j, i) *
           v_dzp_(k-1));
    }

    return ;
  };
 private:
  const ViewDouble1D v_dzp_ = *p_v_dzp;
  const ViewDouble3D v_h0_  = *p_v_h0;
  const ViewDouble3D v_h0l_ = *p_v_h0l;
  const ViewDouble3D v_h0f_ = *p_v_h0f;
  const ViewDouble3D v_psa_ = *p_v_psa;
  const ViewDouble4D v_gg_  = *p_v_gg;
  const ViewDouble4D v_pp_  = *p_v_pp;
  const ViewDouble4D v_ppa_ = *p_v_ppa;
  const ViewDouble4D v_vit_ = *p_v_vit;
};
class FunctorReadyt5 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    thermal(k, j, i, v_ppb_, v_ppc_, v_ppa_, 
        v_alpha_, v_beta_, v_vit_);
    v_gg_(iblock, k, j, i) = - od0_ * G *
        (v_pdensity_(iblock, k, j, i) - v_po_(k) - 1000.0) *
         v_vit_(iblock, k, j, i);

    return ;
  };
  KOKKOS_INLINE_FUNCTION void thermal (const int &k, const int &j, const int &i,
      const ViewDouble4D &v_tt,
      const ViewDouble4D &v_ss,
      const ViewDouble4D &v_pp,
      const ViewDouble4D &v_aa,
      const ViewDouble4D &v_bb,
      const ViewDouble4D &v_mask) const {

    const int iblock = 0;

    const double tmp1 = - v_pp(iblock, k, j, i) 
        / od0_ / 10000.0 * v_mask(iblock, k, j, i);
    const double tmp2 = tmp1 * tmp1;
    const double tmp3 = tmp1 * tmp2;

    const double tt1 = v_tt(iblock, k, j, i);
    const double tt2 = tt1 * tt1;
    const double tt3 = tt1 * tt2;
    const double tt4 = tt1 * tt3;

    const double ss1 = v_ss(iblock, k, j, i) * 1000.0;
    const double ss2 = ss1 * ss1;
    
    v_bb(iblock, k, j, i) = (0.785567e-3 - 0.301985e-5 * tt1 +
        0.555579e-7 * tt2 - 0.415613e-9 * tt3 + ss1 * 
        (-0.356603e-6 + 0.788212e-8 * tt1 + 0.408195e-10 * 
        tmp1 - 0.602281e-15 * tmp2) + ss2 * (0.515032e-8) + 
        tmp1 * (-0.121555e-7 + 0.192867e-9 * tt1 - 0.2131127e-11 *
        tt2) + tmp2 * (0.176621e-12 - 0.175379e-14 *
        tt1) + tmp3 * (0.12155e-17)) * v_mask(iblock, k, j, i);
    
    v_aa(iblock, k, j, i) = (0.665157e-1 + 0.170907e-1 * tt1 -
        0.203814e-3 * tt2 + 0.298357e-5 * tt3 - 0.255019e-7 * tt4 +
        ss1 * (0.378110e-2 - 0.846960e-4 * tt1 - 0.164759e-6 * tmp1 - 
        0.251520e-11 * tmp2) + ss2 * (-0.678662e-5) + tmp1 * 
        (0.380374e-4 - 0.933746e-6 * tt1 + 0.791325e-8 *
        tt2) + 0.512857e-12 * tmp2 * tt2 - 0.302285e-13 * tmp3) *
        v_bb(iblock, k, j, i) * v_mask(iblock, k, j, i);
    return ;
  }
 private:
  const double od0_ = CppPconstMod::od0;
  const ViewDouble1D v_po_       = *p_v_po;
  const ViewDouble4D v_gg_       = *p_v_gg;
  const ViewDouble4D v_ppa_      = *p_v_ppa;
  const ViewDouble4D v_ppb_      = *p_v_ppb;
  const ViewDouble4D v_ppc_      = *p_v_ppc;
  const ViewDouble4D v_vit_      = *p_v_vit;
  const ViewDouble4D v_alpha_    = *p_v_alpha;
  const ViewDouble4D v_beta_     = *p_v_beta;
  const ViewDouble4D v_pdensity_ = *p_v_pdensity;
};

class FunctorReadyt6 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    // if (k < KMM1) {
      if (i == 0 || i == IMT-1) {
        v_ricdt_(iblock, k, j, i)    = 0.0;
        v_ricdttms_(iblock, k, j, i) = 0.0;
      }
      if (i >= 1 && i < (IMT-1)) {
        const double epsln = 1.0e-25;
        v_ricdttms_(iblock, k, j, i) = v_vit_(iblock, k+1, j, i) * G *
            ((v_at_(iblock, 0, k, j, i) - v_at_(iblock, 0, k+1, j, i)) * 
                 v_alpha_(iblock, k+1, j, i) + 1000.0 * 
             (v_at_(iblock, 1, k, j, i) - v_at_(iblock, 1, k+1, j, i)) *
                 v_beta_(iblock, k+1, j, i)) * v_odzt_(k+1);
    
        v_ricdt_(iblock, k, j, i) = v_vit_(iblock, k+1, j, i) / 
            ((v_at_(iblock, 0, k, j, i) - v_at_(iblock, 0, k+1, j, i) + epsln) * 
              v_alpha_(iblock, k+1, j, i)) * 1000.0 *
            ((v_at_(iblock, 1, k, j, i) - v_at_(iblock, 1, k+1, j, i)) * 
              v_beta_(iblock, k+1, j, i));
#ifdef CANUTOMIXOUT
        v_alpha_canuto_(iblock, k, j, i) = v_alpha_(iblock, k, j, i);
        v_beta_canuto_(iblock, k, j, i)  = v_beta_(iblock, k, j, i);
#endif // CANUTOMIXOUT
      }
    // }
    return ;
  }
 private:
  const ViewDouble1D v_odzt_     = *p_v_odzt;
  const ViewDouble4D v_vit_      = *p_v_vit;
  const ViewDouble4D v_ricdt_    = *p_v_ricdt;
  const ViewDouble4D v_alpha_    = *p_v_alpha;
  const ViewDouble4D v_beta_     = *p_v_beta;
  const ViewDouble4D v_ricdttms_ = *p_v_ricdttms;
  const ViewDouble5D v_at_       = *p_v_at;
#ifdef CANUTOMIXOUT
  const ViewDouble4D v_alpha_canuto_ = *p_v_alpha_canuto;
  const ViewDouble4D v_beta_canuto_  = *p_v_beta_canuto;
#endif // CANUTOMIXOUT
};
class FunctorReadyt7 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    grad(iblock, k, j, i, v_dlu_, v_dlv_, v_pp_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad(const int &iblock, const int &k,
      const int &j, const int &i, 
      const ViewDouble4D &v_gradx,
      const ViewDouble4D &v_grady,
      const ViewDouble4D &v_f) const {
    const int bid = 0;
    v_gradx(iblock, k, j, i) = C0;
    v_grady(iblock, k, j, i) = C0;
    if (i >= 1 && j < (NY_BLOCK-1)) {
      if (k <= v_kmu_(bid, j, i) - 1) {
        v_gradx(iblock, k, j, i) = v_dxyur_(bid, j, i, 0) * P5 * 
            (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j, i-1) - 
             v_f(iblock, k, j+1, i-1) + v_f(iblock, k, j, i  ));
      
        v_grady(iblock, k, j, i) = v_dxyur_(bid, j, i, 1) * P5 * 
            (v_f(iblock, k, j+1, i  ) - v_f(iblock, k, j, i-1) + 
             v_f(iblock, k, j+1, i-1) - v_f(iblock, k, j, i  ));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmu_  = *p_v_kmu;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
  const ViewDouble4D v_pp_   = *p_v_pp;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
};

class FunctorReadyt8 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_buoytur_(iblock, j, i) = v_vit_(iblock, 0, j, i) *
        v_nswv_(iblock, j, i) * G * v_alpha_(iblock, 0, j, i) * od0cp_;
    v_buoysol_(iblock, j, i) = v_vit_(iblock, 0, j, i) *
        v_swv_(iblock, j, i) * G * v_alpha_(iblock, 0, j, i) * od0cp_;
    return ;
  }

 private:
  const double od0cp_ = CppPconstMod::od0cp;
  const ViewDouble3D v_swv_     = *p_v_swv;
  const ViewDouble3D v_nswv_    = *p_v_nswv;
  const ViewDouble3D v_buoytur_ = *p_v_buoytur;
  const ViewDouble3D v_buoysol_ = *p_v_buoysol;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble4D v_alpha_   = *p_v_alpha;
};
class FunctorReadyt9 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    vinteg(j, i, v_dlu_, v_pxb_);
    vinteg(j, i, v_dlv_, v_pyb_);
    return ;
  }

  KOKKOS_INLINE_FUNCTION void vinteg (const int &j, const int &i,
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
  const ViewDouble3D v_pxb_  = *p_v_pxb;
  const ViewDouble3D v_pyb_  = *p_v_pyb;
  const ViewDouble3D v_ohbu_ = *p_v_ohbu;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
  const ViewDouble4D v_viv_  = *p_v_viv;
};

class FunctorReadyt10 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {
    const int iblock = 0;
    v_dlu_(iblock, 0, j, i) = 0.0;
    v_dlu_(iblock, 1, j, i) = 0.0;
    for (int k = 0; k < KM; ++k) {
      const double abcd = v_gg_(iblock, k, j, i) 
          * v_ohbt_(iblock, j, i) * v_dzp_(k);
      v_dlu_(iblock, 0, j, i) += abcd;
      v_dlu_(iblock, 1, j, i) += abcd * v_zkt_(k);
    }
    for (int k = 0; k < KM; ++k) {
      v_dlv_(iblock, 0, j, i) = (v_dlu_(iblock, 0, j, i) 
          + v_dlu_(iblock, 1, j, i)
              * v_ohbt_(iblock, j, i)) / G;
      v_dlv_(iblock, 1, j, i) = v_dlu_(iblock, 1, j, i) 
          * v_ohbt_(iblock, j, i) * v_ohbt_(iblock, j, i);
    }
    return ;
  }
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble1D v_zkt_  = *p_v_zkt;
  const ViewDouble3D v_ohbt_ = *p_v_ohbt;
  const ViewDouble4D v_gg_   = *p_v_gg;
  const ViewDouble4D v_dlu_  = *p_v_dlu;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
};

class FunctorReadyt11 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {

    tgrid_to_ugrid(j, i, 0, v_wgp_,  v_dlv_);
    tgrid_to_ugrid(j, i, 1, v_work_, v_dlv_);

    return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
      const int &j, const int &i, const int &n,
          const ViewDouble3D &v_ugrid, const ViewDouble4D &v_tgrid) 
              const {
    const int iblock = 0;
    if (i >= 1 && j <(NY_BLOCK-1)) {
      // v_ugrid(iblock, j, i) = 
      //     v_au0_ (iblock, j, i) * v_tgrid(i  , j  , n, iblock) +
      //     v_aus_ (iblock, j, i) * v_tgrid(i  , j+1, n, iblock) +
      //     v_auw_ (iblock, j, i) * v_tgrid(i-1, j  , n, iblock) +
      //     v_ausw_(iblock, j, i) * v_tgrid(i-1, j+1, n, iblock);
      // v_ugrid(iblock, j, i) = P25 * v_tgrid(i  , j  , n, iblock) +
      //                         P25 * v_tgrid(i  , j+1, n, iblock) +
      //                         P25 * v_tgrid(i-1, j  , n, iblock) +
      //                         P25 * v_tgrid(i-1, j+1, n, iblock);
      v_ugrid(iblock, j, i) = P25 * (
             v_tgrid(iblock, n, j  , i  )
           + v_tgrid(iblock, n, j+1, i  )
           + v_tgrid(iblock, n, j  , i-1)
           + v_tgrid(iblock, n, j+1, i-1));
    // } else {
    //   v_ugrid(iblock, j, i) = 0;
    }
    if (i == 0 || j == (NY_BLOCK-1)) {
      v_ugrid(iblock, j, i) = 0;
    }
    return ;
  }
 private:
  // au0 aus auw ausw = P25
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble3D v_wgp_  = *p_v_wgp;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_dlv_  = *p_v_dlv;
};

// using team_policy = Kokkos::TeamPolicy<>;
// using TeamHandle = typename team_policy::member_type;
// class FunctorReadyt111 {
//  public:
//   KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {
//     const int tid = team_member.league_rank() * team_member.team_size() 
//         + team_member.team_rank();
//     const int i = tid % IMT;
//     const int j = tid / IMT;
//     tgrid_to_ugrid(i, j, 0, v_wgp_,  v_dlv_);
//     tgrid_to_ugrid(i, j, 1, v_work_, v_dlv_);
//   }
//   KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
//       const int &j, const int &i, const int &n,
//           const ViewDouble3D &v_ugrid, const ViewDouble4D &v_tgrid) 
//               const {
//     const int iblock = 0;
//     if (i < IMT && j < JMT) {
//       if (i >= 1 && j <(NY_BLOCK-1)) {
//         v_ugrid(iblock, j, i) = P25 * v_tgrid(i  , j  , n, iblock) +
//                                 P25 * v_tgrid(i  , j+1, n, iblock) +
//                                 P25 * v_tgrid(i-1, j  , n, iblock) +
//                                 P25 * v_tgrid(i-1, j+1, n, iblock);
//       }
//       if (i == 0 || j == (NY_BLOCK-1)) {
//         v_ugrid(iblock, j, i) = 0.0;
//       }
//     }
//     return ;
//   }

//  private:
//   const ViewDouble3D v_wgp_  = *p_v_wgp;
//   const ViewDouble3D v_work_ = *p_v_work;
//   const ViewDouble4D v_dlv_  = *p_v_dlv;
// };
// class FunctorReadyt112 {
//  public:
//   KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {

//     using ScratchView = Kokkos::View<double*, 
//         Kokkos::DefaultExecutionSpace::scratch_memory_space>;

//     ScratchView v_tile(team_member.team_scratch(0), team_member.team_size());

//     const int x = team_member.team_rank() % 16;
//     const int y = team_member.team_rank() / 16;

//     const int i_tmp = (IMT + 14) / 15;

//     const int ii = team_member.league_rank() % i_tmp;
//     const int jj = team_member.league_rank() / i_tmp;

//     const int i = 15 * ii + x;
//     const int j = 7  * jj + y;

//     const int iblock = 0;

//     if (i < IMT && j < JMT) {
//       int k = 0;
//       v_tile(y * 16 + x) = v_dlv_(iblock, k, j, i);

//       team_member.team_barrier();

//       if (x != 0 && y < 7) {
//         double tmp           =  v_tile( y   * 16 + x-1);
//         tmp                 +=  v_tile( y   * 16 + x  );
//         tmp                 +=  v_tile((y+1)* 16 + x-1);
//         v_wgp_(iblock, j, i) = (v_tile((y+1)* 16 + x  ) + tmp) * P25;
//       }
//       if (i == 0 || j == JMT-1) {
//         v_wgp_(iblock, j, i) = 0.0;
//       }

//       k = 1;
//       v_tile(y * 16 + x) = v_dlv_(iblock, k, j, i);

//       team_member.team_barrier();

//       if (x != 0 && y < 7) {
//         double tmp            =  v_tile( y   * 16 + x-1);
//         tmp                  +=  v_tile( y   * 16 + x  );
//         tmp                  +=  v_tile((y+1)* 16 + x-1);
//         v_work_(iblock, j, i) = (v_tile((y+1)* 16 + x  ) + tmp) * P25;
//       }
//       if (i == 0 || j == JMT-1) {
//         v_work_(iblock, j, i) = 0.0;
//       }
//     }
//     return ;
//   }
//  private:
//   const ViewDouble3D v_wgp_  = *p_v_wgp;
//   const ViewDouble3D v_work_ = *p_v_work;
//   const ViewDouble4D v_dlv_  = *p_v_dlv;
// };

class FunctorReadyt12 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_whx_(iblock, j, i) = v_hbx_(iblock, j, i) * 
        v_work_(iblock, j, i) * v_viv_(iblock, 0, j, i);
    v_why_(iblock, j, i) = v_hby_(iblock, j, i) * 
        v_work_(iblock, j, i) * v_viv_(iblock, 0, j, i);

    return ;
  }
 private:
  const ViewDouble3D v_hbx_  = *p_v_hbx;
  const ViewDouble3D v_hby_  = *p_v_hby;
  const ViewDouble3D v_whx_  = *p_v_whx;
  const ViewDouble3D v_why_  = *p_v_why;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_viv_  = *p_v_viv;
};

class FunctorReadyt13 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;

    grad(iblock, 0, j, i, v_work1_, v_work2_, v_psa_);

    v_pay_(iblock, j, i) = - od0_ * v_work2_(iblock, j, i);
    v_pax_(iblock, j, i) = - od0_ * v_work1_(iblock, j, i);

    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &iblock,
      const int &k, const int &j, const int &i, 
      const ViewDouble3D &v_gradx,
      const ViewDouble3D &v_grady, 
      const ViewDouble3D &v_f) const {
    const int bid = 0;
    v_gradx(iblock, j, i) = C0;
    v_grady(iblock, j, i) = C0;
    if (i >= 1 && j < (NY_BLOCK-1)) {
      if (k <= v_kmu_(bid, j, i) - 1) {
        v_gradx(iblock, j, i) = v_dxyur_(bid, j, i, 0) * P5 * 
            (v_f(iblock, j+1, i  ) - v_f(iblock, j, i-1) - 
             v_f(iblock, j+1, i-1) + v_f(iblock, j, i  ));
      
        v_grady(iblock, j, i) = v_dxyur_(bid, j, i, 1) * P5 * 
            (v_f(iblock, j+1, i  ) - v_f(iblock, j, i-1) + 
             v_f(iblock, j+1, i-1) - v_f(iblock, j, i  ));
      }
    }
    return ;
  }
 private:
  const double od0_ = CppPconstMod::od0;
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble3D v_pax_   = *p_v_pax;
  const ViewDouble3D v_pay_   = *p_v_pay;
  const ViewDouble3D v_psa_   = *p_v_psa;
  const ViewDouble3D v_work1_ = *p_v_work1;
  const ViewDouble3D v_work2_ = *p_v_work2;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
};

KOKKOS_REGISTER_FOR_3D(FunctorReadyt1,  FunctorReadyt1)
KOKKOS_REGISTER_FOR_3D(FunctorReadyt2,  FunctorReadyt2)
KOKKOS_REGISTER_FOR_3D(FunctorReadyt3,  FunctorReadyt3)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt4,  FunctorReadyt4)
KOKKOS_REGISTER_FOR_3D(FunctorReadyt5,  FunctorReadyt5)
KOKKOS_REGISTER_FOR_3D(FunctorReadyt6,  FunctorReadyt6)
KOKKOS_REGISTER_FOR_3D(FunctorReadyt7,  FunctorReadyt7)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt8,  FunctorReadyt8)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt9,  FunctorReadyt9)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt10, FunctorReadyt10)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt11, FunctorReadyt11)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt12, FunctorReadyt12)
KOKKOS_REGISTER_FOR_2D(FunctorReadyt13, FunctorReadyt13)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_READYT_HPP_
