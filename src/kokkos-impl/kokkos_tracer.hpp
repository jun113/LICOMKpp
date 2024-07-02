#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_TRACER_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_TRACER_HPP_
#include "../head/def-undef.h"

#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"

#ifndef BIHAR
#include "../head/cpp_hmix_del2.h"
#else  // BIHAR
#include "../head/cpp_hmix_del4.h"
#endif // BIHAR

#ifdef ISO
#include "../head/cpp_isopyc_mod.h"
#endif // ISO

#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"

#ifdef COUP
#include "../head/kokkos_buf_mod.h"
#endif // COUP                  

#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_grid.h"

#ifndef BIHAR
#include "../head/kokkos_hmix_del2.h"
#else  // BIHAR
#include "../head/kokkos_hmix_del4.h"
#endif // BIHAR

#ifdef ISO
#include "../head/kokkos_isopyc_mod.h"
#endif // ISO

#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_pmix_mod.h"
#include "../head/kokkos_work_mod.h"
#include "../head/kokkos_tracer_mod.h" 
#include "../head/kokkos_tmp_var.h" 
#include "../head/kokkos_config.hpp"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include "Kokkos_Core.hpp"

#include <string>

using CppConstantMod::C0;
using CppConstantMod::P5;
using CppConstantMod::P25;
using CppParamMod::KM;
// using CppParamMod::KMM1;
// using CppParamMod::KMP1;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JST;
using CppParamMod::JET;
using CppParamMod::NTRA;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
// using CppParamMod::MAX_BLOCKS_CLINIC;

using KokkosDynMod   ::p_v_h0f;
using KokkosDynMod   ::p_v_h0l;
using KokkosDynMod   ::p_v_utf;
using KokkosDynMod   ::p_v_vtf;
using KokkosDynMod   ::p_v_utl;
using KokkosDynMod   ::p_v_vtl;
using KokkosDynMod   ::p_v_ws;
using KokkosForcMod  ::p_v_sss;
using KokkosForcMod  ::p_v_sst;
using KokkosForcMod  ::p_v_swv;
using KokkosForcMod  ::p_v_nswv;
using KokkosForcMod  ::p_v_restore;
// using KokkosGrid     ::p_v_au0;
// using KokkosGrid     ::p_v_aus;
// using KokkosGrid     ::p_v_auw;
// using KokkosGrid     ::p_v_ausw;
using KokkosGrid     ::p_v_dxu;
using KokkosGrid     ::p_v_dyu;
using KokkosGrid     ::p_v_hue;
using KokkosGrid     ::p_v_hun;
using KokkosGrid     ::p_v_hts;
using KokkosGrid     ::p_v_htw;
using KokkosGrid     ::p_v_h_tu_swen;
using KokkosGrid     ::p_v_kmt;
using KokkosGrid     ::p_v_kmtn;
using KokkosGrid     ::p_v_kmts;
using KokkosGrid     ::p_v_kmte;
using KokkosGrid     ::p_v_kmtw;
using KokkosGrid     ::p_v_kmt_nsew;
using KokkosGrid     ::p_v_tlat;
using KokkosGrid     ::p_v_tarea;
using KokkosGrid     ::p_v_tarea_r;

using KokkosPconstMod:: p_v_c;
using KokkosPconstMod:: p_v_to;
using KokkosPconstMod:: p_v_so;
using KokkosPconstMod:: p_v_akt;
using KokkosPconstMod:: p_v_dzp;
using KokkosPconstMod:: p_v_vit;
using KokkosPconstMod:: p_v_zkp;
using KokkosPconstMod:: p_v_zkt;
using KokkosPconstMod:: p_v_rrd1;
using KokkosPconstMod:: p_v_rrd2;
using KokkosPconstMod:: p_v_odzp;
using KokkosPconstMod:: p_v_odzt;
using KokkosPconstMod:: p_v_odz_pt;
using KokkosPconstMod:: p_v_ohbu;
using KokkosPconstMod:: p_v_ohbt;

using KokkosWorkMod  :: p_v_tf;
using KokkosWorkMod  :: p_v_uk;
using KokkosWorkMod  :: p_v_vk;
using KokkosWorkMod  :: p_v_stf;
using KokkosWorkMod  :: p_v_wka;
using KokkosWorkMod  :: p_v_wkb;
using KokkosWorkMod  :: p_v_wkc;
using KokkosWorkMod  :: p_v_wkd;
using KokkosWorkMod  :: p_v_work;

// using KokkosTracerMod:: p_v_ax;
// using KokkosTracerMod:: p_v_ay;
// using KokkosTracerMod:: p_v_az;
using KokkosTracerMod:: p_v_at;
// using KokkosTracerMod:: p_v_dx;
// using KokkosTracerMod:: p_v_dz;
using KokkosTracerMod:: p_v_atb;
using KokkosTracerMod:: p_v_net;
using KokkosTracerMod:: p_v_tend;
using KokkosTracerMod:: p_v_dt_diff;
using KokkosTracerMod:: p_v_penetrate;

#ifdef BIHAR
using KokkosHmixDel4::  p_v_ahf;
using KokkosHmixDel4::  p_v_dtn;
using KokkosHmixDel4::  p_v_dts;
using KokkosHmixDel4::  p_v_dte;
using KokkosHmixDel4::  p_v_dtw;
using KokkosHmixDel4::  p_v_dt_nsew;
#else // BIHAR
using KokkosHmixDel2::  p_v_dtn;
using KokkosHmixDel2::  p_v_dts;
using KokkosHmixDel2::  p_v_dte;
using KokkosHmixDel2::  p_v_dtw;
#endif // BIHAR

#ifdef SOLAR
using KokkosPmixMod  :: p_v_pen;
#endif // SOLAR
#ifdef SOLARCHLORO
using KokkosPmixMod  :: p_v_pen_chl;
#endif // SOLARCHLORO

#ifdef COUP                  
using KokkosBufMod   ::p_v_ifrac;
using KokkosForcMod  ::p_v_ssf;
using KokkosForcMod  ::p_v_tsf;
#else  // COUP
#ifdef FRC_CORE
using KokkosForcMod  ::p_v_fresh;
using KokkosForcMod  ::p_v_seaice;
#else  // FRC_CORE
using KokkosForcMod  ::p_v_dqdt;
#endif // FRC_CORE
#endif // COUP

#ifdef ISO
using CppIsopycMod::NRPL;
using KokkosWorkMod::p_v_temp;
using KokkosWorkMod::p_v_work_1;
using KokkosWorkMod::p_v_work_2;
using KokkosWorkMod::p_v_work_3;

using KokkosIsopycMod::p_v_e;
using KokkosIsopycMod::p_v_k1;
using KokkosIsopycMod::p_v_k2;
using KokkosIsopycMod::p_v_k3;
using KokkosIsopycMod::p_v_dzr;
using KokkosIsopycMod::p_v_dzw;
using KokkosIsopycMod::p_v_dzwr;
using KokkosIsopycMod::p_v_rhoi;
using KokkosIsopycMod::p_v_tmask;
using KokkosIsopycMod::p_v_ahisop;
using KokkosIsopycMod::p_v_athkdf;
using KokkosIsopycMod::p_v_fzisop;
using KokkosIsopycMod::p_v_krplin;
using KokkosIsopycMod::p_v_kisrpl;
using KokkosIsopycMod::p_v_adv_vbtiso;
using KokkosIsopycMod::p_v_adv_vetiso;
using KokkosIsopycMod::p_v_adv_vntiso;

using KokkosTracerMod::p_v_ax_iso;
using KokkosTracerMod::p_v_ay_iso;
using KokkosTracerMod::p_v_az_iso;
using KokkosTracerMod::p_v_dx_iso;
using KokkosTracerMod::p_v_dy_iso;
using KokkosTracerMod::p_v_dz_iso;
using KokkosTracerMod::p_v_ddy_iso;

#endif // ISO

using KokkosTmpVar::p_v_nn;
using KokkosTmpVar::p_v_xs;
using KokkosTmpVar::p_v_at_00_max_min;
using KokkosTmpVar::p_v_adv_tt;
using KokkosTmpVar::p_v_uv_ws_face;
using KokkosTmpVar::p_v_vtl_ori;

#ifdef BIHAR
using KokkosTmpVar::p_v_dt2k;
using KokkosTmpVar::p_v_c_cnsew;
#endif // BIHAR

using KokkosTmpVar::p_v_wp1;
using KokkosTmpVar::p_v_wp2;
using KokkosTmpVar::p_v_wp3;
using KokkosTmpVar::p_v_wp4;
using KokkosTmpVar::p_v_wk1;
using KokkosTmpVar::p_v_wk2;

// !---------------------------------------------------------------------
// !     PREPARATION FOR VERTICAL ADVECTIVE TERM
// !---------------------------------------------------------------------
class FunctorTracer1 {
 public: 
  FunctorTracer1 (const double &aa) : aa_(aa) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_h0f_(iblock, j, i) *= onbc_;

    v_stf_(iblock, j, i) = aa_ * v_h0f_(iblock, j, i)
        + (1.0 - aa_) * v_h0l_(iblock, j, i);
    return ;
  }
 private:
  const double aa_;
  const double onbc_ = CppPconstMod::onbc;
  const ViewDouble3D v_h0f_ = *p_v_h0f;
  const ViewDouble3D v_h0l_ = *p_v_h0l;
  const ViewDouble3D v_stf_ = *p_v_stf;
};

class FunctorTracer2 {
 public:
  FunctorTracer2 (const double &aa) : aa_(aa) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_utf_(iblock, k, j, i) *= oncc_;
    v_vtf_(iblock, k, j, i) *= oncc_;

    v_wkd_(iblock, k, j, i) = aa_ * v_utf_(iblock, k, j, i)
        + (1.0 - aa_) * v_utl_(iblock, k, j, i);
    v_wkb_(iblock, k, j, i) = aa_ * v_vtf_(iblock, k, j, i)
        + (1.0 - aa_) * v_vtl_(iblock, k, j, i);
    return ;
  }

 private:
  const double aa_;
  const double oncc_ = CppPconstMod::oncc;
  const ViewDouble4D v_utf_ = *p_v_utf;
  const ViewDouble4D v_utl_ = *p_v_utl;
  const ViewDouble4D v_vtf_ = *p_v_vtf;
  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_wkb_ = *p_v_wkb;
  const ViewDouble4D v_wkd_ = *p_v_wkd;
};
class FunctorTracer3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    upwell_1 (iblock, j, i, v_stf_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_1(const int &iblock,
      const int &j, const int &i, const ViewDouble3D &v_h0wk)
          const {
    tgrid_to_ugrid (iblock, j, i, v_work_, v_h0wk);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
      const int &iblock, const int &j, const int &i,
          const ViewDouble3D &v_ugrid, 
              const ViewDouble3D &v_tgrid) 
          const {
    if (i >= 1 && j <(NY_BLOCK-1)) {
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
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble3D v_stf_  = *p_v_stf;
  const ViewDouble3D v_work_ = *p_v_work;
};
class FunctorTracer4 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    upwell_2 (iblock, k, j, i, v_wkd_, v_wkb_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_2 (const int &iblock, const int &k, 
      const int &j, const int &i, const ViewDouble4D &v_uwk,
          const ViewDouble4D& v_vwk) const {
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
  const ViewDouble4D v_uk_   = *p_v_uk;
  const ViewDouble4D v_vk_   = *p_v_vk;
  const ViewDouble4D v_wkb_  = *p_v_wkb;
  const ViewDouble4D v_wkd_  = *p_v_wkd;
};
class FunctorTracer5 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    upwell_3 (iblock, k, j, i);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void upwell_3 (const int &iblock,
      const int &k, const int &j, const int &i) const {
    div(iblock, k, j, i, v_wka_, v_uk_, v_vk_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void div(const int &iblock, const int &k, 
      const int &j, const int &i, const ViewDouble4D &v_div_out, 
          const ViewDouble4D &v_ux, const ViewDouble4D &v_uy) const {
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
class FunctorTracer6 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    upwell_4(iblock, j, i, v_stf_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION 
  void upwell_4 (const int &iblock, const int &j, 
      const int &i, const ViewDouble3D &v_h0wk) const {

    double work = 0.0;

    if (i >= 1 && i < (IMT-1) && j >= 1 && j < (JMT-1)) {
      for (int k = 0; k < KM; ++k) {
        work -= v_dzp_(k)
            * v_wka_(iblock, k, j, i) * v_vit_(iblock, k, j, i);
      }
      double ws = v_ws_(iblock, 0, j, i);
      for (int k = 1; k < KM; ++k) {
        ws = v_vit_(iblock, k, j, i) * (ws + v_dzp_(k-1) 
            * (work * v_ohbt_(iblock, j, i)
                + v_wka_(iblock, k-1, j, i)));
        v_ws_(iblock, k, j, i) = ws;
      }

      work = 1.0 / (1.0 + v_h0wk(iblock, j, i)
          * v_ohbt_(iblock, j, i));

      for (int k = 1; k < KM; ++k) {
        v_ws_(iblock, k, j, i) *= work;
      }
    }
    return ;
  }
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_stf_  = *p_v_stf;
  const ViewDouble3D v_ohbt_ = *p_v_ohbt;
  const ViewDouble4D v_ws_   = *p_v_ws;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

// !@@@  COMPUTING DIFFUSION COEFFICIENT
class FunctorTracer7 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
#ifdef ISO
    v_wkc_(iblock, k, j, i) = ahv_ + v_ahisop_(iblock, j, i)
        * v_k3_(iblock, 2, j, k, i);
#else  // ISO
    v_wkc_(iblock, k, j, i) = ahv_;
#endif // ISO

#ifndef CANUTO
    v_akt_(iblock, 0, k, j, i) = v_wkc_(iblock, k, j, i);
    v_akt_(iblock, 1, k, j, i) = v_wkc_(iblock, k, j, i);
#endif // CANUTO
    return ;
  }

 private:
  const double ahv_ = CppPconstMod::ahv;
  const ViewDouble4D v_wkc_ = *p_v_wkc;
#ifdef ISO
  const ViewDouble3D v_ahisop_ = *p_v_ahisop;
  const ViewDouble5D v_k3_     = *p_v_k3;
#endif // ISO

#ifndef CANUTO
  const ViewDouble5D v_akt_ = *p_v_akt;
#endif // CANUTO
};

class FunctorTracer8 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_tf_(iblock, k, j, i) = 0.0;
    v_adv_tt_(k, j, i)     = 0.0;
    return ;
  }
 private:
  const ViewDouble3D v_adv_tt_ = *p_v_adv_tt;
  const ViewDouble4D v_tf_     = *p_v_tf;
};

class FuncAdvTraCenTsp1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advction_tracer_centered_tspas_1 (iblock, k, j, i, v_wkd_, v_wkb_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void advction_tracer_centered_tspas_1 (
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_uuu, const ViewDouble4D &v_vvv) const {
    if (j >= 1 && j < (JMT-1)) {
      v_uv_ws_face_(k, j, i, 0) = (v_uuu(iblock, k, j-1, i) 
          + v_uuu(iblock, k, j, i)) * v_htw_(iblock, j, i) * P25;
      // v_uv_ws_face_(0, k, j, i) = v_u_wface_(k, j, i);
    }
    if (i >= 1 && i < (IMT-1)) {
      v_uv_ws_face_(k, j, i, 1) = (v_vvv(iblock, k, j, i) 
          + v_vvv(iblock, k, j, i+1)) * v_hts_(iblock, j, i) * P25;
      // v_uv_ws_face_(1, k, j, i) = v_v_sface_(k, j, i);
    }
    return ;
  }
 private:
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble4D v_uv_ws_face_ = *p_v_uv_ws_face;
  const ViewDouble4D v_wkb_     = *p_v_wkb;
  const ViewDouble4D v_wkd_     = *p_v_wkd;
};

class FuncAdvTraTsp2 {
 public:
  FuncAdvTraTsp2 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advection_tracer_tspas_2 (iblock, n_, k, j, i, v_ws_, v_at_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void advection_tracer_tspas_2 (
      const int &iblock, const int &n, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_www, const ViewDouble5D &v_ttt) const {

    const double adv_x0 = 
        (v_ttt(iblock, n, k, j, i+1) + v_ttt(iblock, n, k, j, i  ))
            * v_uv_ws_face_(k, j, i+1, 0) * v_tarea_r_(iblock, j, i)
      - (v_ttt(iblock, n, k, j, i  ) + v_ttt(iblock, n, k, j, i-1))
            * v_uv_ws_face_(k, j, i  , 0) * v_tarea_r_(iblock, j, i);
             
    const double adv_y0 = (
        (v_ttt(iblock, n, k, j+1, i) + v_ttt(iblock, n, k, j  , i))
            * v_uv_ws_face_(k, j  , i, 1)
      - (v_ttt(iblock, n, k, j  , i) + v_ttt(iblock, n, k, j-1, i))
            * v_uv_ws_face_(k, j-1, i, 1)) * v_tarea_r_(iblock, j, i);
    
    const double adv_xy1 = - dts_ * (v_ttt(iblock, n, k, j, i+1) - v_ttt(iblock, n, k, j, i))
        * 2.0 * v_tarea_r_(iblock, j, i)
            * v_uv_ws_face_(k, j, i+1, 0) * v_uv_ws_face_(k, j, i+1, 0)
                / (v_h_tu_swen_(iblock, j, i+1, 0, 0) * v_h_tu_swen_(iblock, j, i+1, 0, 1));

    const double adv_xy2 = dts_ * (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j, i-1))
        * 2.0 * v_tarea_r_(iblock, j, i)
            * v_uv_ws_face_(k, j, i  , 0) * v_uv_ws_face_(k, j, i  , 0)
                / (v_h_tu_swen_(iblock, j, i, 0, 0) * v_h_tu_swen_(iblock, j, i, 0, 1));

    const double adv_xy3 = - dts_ * (v_ttt(iblock, n, k, j+1, i) - v_ttt(iblock, n, k, j, i))
        * 2.0 * v_tarea_r_(iblock, j, i)
            * v_uv_ws_face_(k, j  , i, 1) * v_uv_ws_face_(k, j  , i, 1)
                / (v_h_tu_swen_(iblock, j, i, 1, 0) * v_h_tu_swen_(iblock, j, i, 1, 1));

    const double adv_xy4 = dts_ * (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j-1, i))
        * 2.0 * v_tarea_r_(iblock, j, i)
            * v_uv_ws_face_(k, j-1, i, 1) * v_uv_ws_face_(k, j-1, i, 1)
                / (v_h_tu_swen_(iblock, j-1, i, 1, 0) * v_h_tu_swen_(iblock, j-1, i, 1, 1));

    const double adv_c1 = - v_ttt(iblock, n, k, j, i)
        * (v_uv_ws_face_(k, j, i+1, 0) - v_uv_ws_face_(k, j, i, 0))
            * v_tarea_r_(iblock, j, i) * 2.0;
             
    const double adv_c2 = - v_ttt(iblock, n, k, j, i)
        * (v_uv_ws_face_(k, j, i, 1) - v_uv_ws_face_(k, j-1, i, 1))
            * v_tarea_r_(iblock, j, i) * 2.0;

    double adv_za, adv_zc;
    double adv_zb1, adv_zb2;
    const double half_odzp = 0.5 * v_odz_pt_(k, 0);
    if (k == 0) {
      adv_za = - half_odzp * v_www(iblock, 1, j, i)
          * (v_ttt(iblock, n, 1, j, i) + v_ttt(iblock, n, 0, j, i));

      adv_zb1 = 0.0;
      adv_zb2 = half_odzp * v_www(iblock, 1, j, i) * v_www(iblock, 1, j, i) 
          * v_odz_pt_(1, 1) * (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 1, j, i)) * dts_;

      adv_zc = v_odz_pt_(0, 0) * v_ttt(iblock, n, 0, j, i) * v_www(iblock, 1, j, i);
    } else if (k == KM-1){
      adv_za = half_odzp * v_www(iblock, KM-1, j, i)
          * (v_ttt(iblock, n, KM-1, j, i) + v_ttt(iblock, n, KM-2, j, i));
           
      adv_zb1 = - half_odzp * v_www(iblock, KM-1, j, i) * v_www(iblock, KM-1, j, i)
          * v_odz_pt_(KM-1, 1) * (v_ttt(iblock, n, KM-2, j, i) 
              - v_ttt(iblock, n, KM-1, j, i)) * dts_;
      adv_zb2 = 0.0;

      adv_zc = - v_odz_pt_(KM-1, 0) * v_ttt(iblock, n, KM-1, j, i)
          * v_www(iblock, KM-1, j, i);
    } else {
      adv_za = 
          half_odzp * v_www(iblock, k  , j, i)
              * (v_ttt(iblock, n, k, j, i) + v_ttt(iblock, n, k-1, j, i))
        - half_odzp * v_www(iblock, k+1, j, i)
              * (v_ttt(iblock, n, k, j, i) + v_ttt(iblock, n, k+1, j, i));
               
      adv_zb1 = - half_odzp * v_www(iblock, k  , j, i) * v_www(iblock, k  , j, i)
          * v_odz_pt_(k, 1) * (v_ttt(iblock, n, k-1, j, i) 
              - v_ttt(iblock, n, k  , j, i)) * dts_;
      adv_zb2 =   half_odzp * v_www(iblock, k+1, j, i) * v_www(iblock, k+1, j, i)
          * v_odz_pt_(k+1, 1) * (v_ttt(iblock, n, k, j, i) 
              - v_ttt(iblock, n, k+1, j, i)) * dts_;
      
      adv_zc = - v_odz_pt_(k, 0) * v_ttt(iblock, n, k, j, i)
          * (v_www(iblock, k, j, i) - v_www(iblock, k+1, j, i));
    }
    const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
    const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
    const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

    v_at_00_max_min_(k, j, i, 0) = v_ttt(iblock, n, k, j, i)
        + (adv_xx + adv_yy + adv_zz) * dts_;
    return ;
  }

 private:
  const int n_;
  const double dts_ = CppPconstMod::dts;
  const ViewDouble2D v_odz_pt_        = *p_v_odz_pt;
  const ViewDouble3D v_tarea_r_       = *p_v_tarea_r;
  const ViewDouble4D v_ws_            = *p_v_ws;
  const ViewDouble4D v_uv_ws_face_    = *p_v_uv_ws_face;
  const ViewDouble4D v_at_00_max_min_ = *p_v_at_00_max_min;
  const ViewDouble5D v_at_            = *p_v_at;
  const ViewDouble5D v_h_tu_swen_     = *p_v_h_tu_swen;
};

class FuncAdvTraTsp3 {
 public:
  FuncAdvTraTsp3 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advection_tracer_tspas_3 (iblock, n_, k, j, i, v_at_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void advection_tracer_tspas_3 (
      const int &iblock, const int &n, const int &k, const int &j, const int &i,
          const ViewDouble5D &v_ttt) const {

    const double wt1 = - 1.0e10;
    const double wt2 = + 1.0e10;

    if (k == 0) {
      v_at_00_max_min_(k, j, i, 1) = myMax(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt1,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt1,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt1,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt1,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt1,
        v_ttt(iblock, n, k+1, j  , i  ) * v_vit_(iblock, k+1, j  , i  )
            + (1.0 - v_vit_(iblock, k+1, j  , i  )) * wt1);

      v_at_00_max_min_(k, j, i, 2) = myMin(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt2,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt2,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt2,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt2,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt2,
        v_ttt(iblock, n, k+1, j  , i  ) * v_vit_(iblock, k+1, j  , i  )
            + (1.0 - v_vit_(iblock, k+1, j  , i  )) * wt2);
    } else if (k == KM-1) {
      v_at_00_max_min_(k, j, i, 1) = myMax(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt1,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt1,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt1,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt1,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt1,
        v_ttt(iblock, n, k-1, j  , i  ) * v_vit_(iblock, k-1, j  , i  )
            + (1.0 - v_vit_(iblock, k-1, j  , i  )) * wt1);

      v_at_00_max_min_(k, j, i, 2) = myMin(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt2,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt2,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt2,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt2,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt2,
        v_ttt(iblock, n, k-1, j  , i  ) * v_vit_(iblock, k-1, j  , i  )
            + (1.0 - v_vit_(iblock, k-1, j  , i  )) * wt2);
    } else {
      v_at_00_max_min_(k, j, i, 1) = myMax(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt1,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt1,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt1,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt1,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt1,
        v_ttt(iblock, n, k-1, j  , i  ) * v_vit_(iblock, k-1, j  , i  )
            + (1.0 - v_vit_(iblock, k-1, j  , i  )) * wt1,
        v_ttt(iblock, n, k+1, j  , i  ) * v_vit_(iblock, k+1, j  , i  )
            + (1.0 - v_vit_(iblock, k+1, j  , i  )) * wt1);

      v_at_00_max_min_(k, j, i, 2) = myMin(
        v_ttt(iblock, n, k  , j  , i  ) * v_vit_(iblock, k  , j  , i  )
            + (1.0 - v_vit_(iblock, k  , j  , i  )) * wt2,
        v_ttt(iblock, n, k  , j-1, i  ) * v_vit_(iblock, k  , j-1, i  )
            + (1.0 - v_vit_(iblock, k  , j-1, i  )) * wt2,
        v_ttt(iblock, n, k  , j+1, i  ) * v_vit_(iblock, k  , j+1, i  )
            + (1.0 - v_vit_(iblock, k  , j+1, i  )) * wt2,
        v_ttt(iblock, n, k  , j  , i-1) * v_vit_(iblock, k  , j  , i-1)
            + (1.0 - v_vit_(iblock, k  , j  , i-1)) * wt2,
        v_ttt(iblock, n, k  , j  , i+1) * v_vit_(iblock, k  , j  , i+1)
            + (1.0 - v_vit_(iblock, k  , j  , i+1)) * wt2,
        v_ttt(iblock, n, k-1, j  , i  ) * v_vit_(iblock, k-1, j  , i  )
            + (1.0 - v_vit_(iblock, k-1, j  , i  )) * wt2,
        v_ttt(iblock, n, k+1, j  , i  ) * v_vit_(iblock, k+1, j  , i  )
            + (1.0 - v_vit_(iblock, k+1, j  , i  )) * wt2);
    }
    return ;
  }
 private:
  const int n_;
  const ViewDouble4D v_vit_           = *p_v_vit;
  const ViewDouble4D v_at_00_max_min_ = *p_v_at_00_max_min;
  const ViewDouble5D v_at_            = *p_v_at;

  template<typename T> KOKKOS_INLINE_FUNCTION 
      T myMax(const T &val1, const T &val2) const {
    return val1 > val2 ? val1 : val2;
  }

  template<typename T, typename... Args> KOKKOS_INLINE_FUNCTION 
      T myMax(const T &val, const Args &... arg) const {
    T result = myMax(arg...);
    return myMax(val, result);
  }

  template<typename T> KOKKOS_INLINE_FUNCTION 
      T myMin(const T &val1, const T &val2) const {
    return val1 < val2 ? val1 : val2;
  }

  template<typename T, typename... Args> KOKKOS_INLINE_FUNCTION 
      T myMin(const T &val, const Args &... arg) const {
    T result = myMin(arg...);
    return myMin(val, result);
  }
};

class FuncAdvTraTsp4 {
 public:
  FuncAdvTraTsp4 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {

    const int iblock = 0;
    advection_tracer_tspas_4 (iblock, n_, j, i, v_ws_, v_at_, v_adv_tt_);

    return ;
  }

  KOKKOS_INLINE_FUNCTION void advection_tracer_tspas_4 (
      const int &iblock, const int &n, const int &j, const int &i,
          const ViewDouble4D &v_www,
          const ViewDouble5D &v_ttt,
          const ViewDouble3D &v_adv_tt) const {

    double adv_xy1;
    if (v_at_00_max_min_(0, j  , i  , 0) > v_at_00_max_min_(0, j  , i  , 1) ||
        v_at_00_max_min_(0, j  , i  , 0) < v_at_00_max_min_(0, j  , i  , 2) ||
        v_at_00_max_min_(0, j  , i+1, 0) > v_at_00_max_min_(0, j  , i+1, 1) ||
        v_at_00_max_min_(0, j  , i+1, 0) < v_at_00_max_min_(0, j  , i+1, 2)) {
      adv_xy1 = - 
          (v_ttt(iblock, n, 0, j, i+1) - v_ttt(iblock, n, 0, j, i))
              * fabs(v_uv_ws_face_(0, j, i+1, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy1 = - dts_
        * (v_ttt(iblock, n, 0, j, i+1) - v_ttt(iblock, n, 0, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(0, j, i+1, 0) * v_uv_ws_face_(0, j, i+1, 0)
                / (v_h_tu_swen_(iblock, j, i+1, 0, 0) * v_h_tu_swen_(iblock, j, i+1, 0, 1));
                // / (v_htw_(iblock, j, i+1) * v_hun_(iblock, j, i+1));
    }

    double adv_xy2;
    if (v_at_00_max_min_(0, j  , i  , 0) > v_at_00_max_min_(0, j  , i  , 1) ||
        v_at_00_max_min_(0, j  , i  , 0) < v_at_00_max_min_(0, j  , i  , 2) ||
        v_at_00_max_min_(0, j  , i-1, 0) > v_at_00_max_min_(0, j  , i-1, 1) ||
        v_at_00_max_min_(0, j  , i-1, 0) < v_at_00_max_min_(0, j  , i-1, 2)) {
      adv_xy2 =   
          (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 0, j, i-1))
              * fabs(v_uv_ws_face_(0, j, i, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy2 = dts_
        * (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 0, j, i-1))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(0, j, i, 0) * v_uv_ws_face_(0, j, i, 0)
                / (v_h_tu_swen_(iblock, j, i, 0, 0) * v_h_tu_swen_(iblock, j, i, 0, 1));
                // / (v_htw_(iblock, j, i) * v_hun_(iblock, j, i));
    }

    double adv_xy3;
    if (v_at_00_max_min_(0, j  , i  , 0) > v_at_00_max_min_(0, j  , i  , 1) ||
        v_at_00_max_min_(0, j  , i  , 0) < v_at_00_max_min_(0, j  , i  , 2) ||
        v_at_00_max_min_(0, j+1, i  , 0) > v_at_00_max_min_(0, j+1, i  , 1) ||
        v_at_00_max_min_(0, j+1, i  , 0) < v_at_00_max_min_(0, j+1, i  , 2)) {
      adv_xy3 = -
          (v_ttt(iblock, n, 0, j+1, i) - v_ttt(iblock, n, 0, j, i))
              * fabs(v_uv_ws_face_(0, j, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy3 = - dts_
        * (v_ttt(iblock, n, 0, j+1, i) - v_ttt(iblock, n, 0, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(0, j, i, 1) * v_uv_ws_face_(0, j, i, 1)
                / (v_h_tu_swen_(iblock, j, i, 1, 0) * v_h_tu_swen_(iblock, j, i, 1, 1));
                // / (v_hts_(iblock, j, i) * v_hue_(iblock, j, i));
    }

    double adv_xy4;
    if (v_at_00_max_min_(0, j  , i  , 0) > v_at_00_max_min_(0, j  , i  , 1) ||
        v_at_00_max_min_(0, j  , i  , 0) < v_at_00_max_min_(0, j  , i  , 2) ||
        v_at_00_max_min_(0, j-1, i  , 0) > v_at_00_max_min_(0, j-1, i  , 1) ||
        v_at_00_max_min_(0, j-1, i  , 0) < v_at_00_max_min_(0, j-1, i  , 2)) {
      adv_xy4 =  
          (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 0, j-1, i))
              * fabs(v_uv_ws_face_(0, j-1, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy4 = dts_
        * (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 0, j-1, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(0, j-1, i, 1) * v_uv_ws_face_(0, j-1, i, 1)
                / (v_h_tu_swen_(iblock, j-1, i, 1, 0) * v_h_tu_swen_(iblock, j-1, i, 1, 1));
                // / (v_hts_(iblock, j-1, i) * v_hue_(iblock, j-1, i));
    }

    const double adv_zb1 = 0.0;

    double adv_zb2;
    if (v_at_00_max_min_(0, j  , i  , 0) > v_at_00_max_min_(0, j  , i  , 1) ||
        v_at_00_max_min_(0, j  , i  , 0) < v_at_00_max_min_(0, j  , i  , 2) ||
        v_at_00_max_min_(1, j  , i  , 0) > v_at_00_max_min_(1, j  , i  , 1) ||
        v_at_00_max_min_(1, j  , i  , 0) < v_at_00_max_min_(1, j  , i  , 2)) {
      adv_zb2 = 0.5 * fabs(v_ws_(iblock, 1, j, i)) * v_odz_pt_(0, 0)
          * (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 1, j, i));
    } else {
      adv_zb2 = 0.5 * v_odz_pt_(0, 0) * v_www(iblock, 1, j, i) * v_www(iblock, 1, j, i) 
          * v_odz_pt_(1, 1) * (v_ttt(iblock, n, 0, j, i) - v_ttt(iblock, n, 1, j, i)) * dts_;
    }

    const double adv_za = - 0.5 * v_odz_pt_(0, 0) * v_www(iblock, 1, j, i)
        * (v_ttt(iblock, n, 1, j, i) + v_ttt(iblock, n, 0, j, i));

    const double adv_zc = v_odz_pt_(0, 0) * v_ttt(iblock, n, 0, j, i)
        * v_ws_(iblock, 1, j, i);

    const double adv_c1 = - v_ttt(iblock, n, 0, j, i)
        * (v_uv_ws_face_(0, j  , i+1, 0) - v_uv_ws_face_(0, j  , i  , 0))
            * v_tarea_r_(iblock, j, i) * 2.0;
    const double adv_c2 = - v_ttt(iblock, n, 0, j, i)
        * (v_uv_ws_face_(0, j  , i  , 1) - v_uv_ws_face_(0, j-1, i  , 1))
            * v_tarea_r_(iblock, j, i) * 2.0;

    const double adv_x0 =
        (v_ttt(iblock, n, 0, j  , i+1) + v_ttt(iblock, n, 0, j  , i  ))
            * v_uv_ws_face_(0, j  , i+1, 0) * v_tarea_r_(iblock, j, i) 
      - (v_ttt(iblock, n, 0, j  , i  ) + v_ttt(iblock, n, 0, j  , i-1))
            * v_uv_ws_face_(0, j  , i  , 0) * v_tarea_r_(iblock, j, i);
    const double adv_y0 = (
        (v_ttt(iblock, n, 0, j+1, i  ) + v_ttt(iblock, n, 0, j  , i  ))
            * v_uv_ws_face_(0, j  , i  , 1) 
      - (v_ttt(iblock, n, 0, j  , i  ) + v_ttt(iblock, n, 0, j-1, i  ))
            * v_uv_ws_face_(0, j-1, i  , 1)) * v_tarea_r_(iblock, j, i);

    const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
    const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
    const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

    v_adv_tt(0, j, i) = adv_xx + adv_yy + adv_zz;
    return ;
  }

 private:
  const int n_;
  const double dts_ = CppPconstMod::dts;
  const ViewDouble2D v_odz_pt_        = *p_v_odz_pt;
  const ViewDouble3D v_adv_tt_        = *p_v_adv_tt;
  const ViewDouble3D v_tarea_r_       = *p_v_tarea_r;
  const ViewDouble4D v_ws_            = *p_v_ws;
  const ViewDouble4D v_uv_ws_face_    = *p_v_uv_ws_face;
  const ViewDouble4D v_at_00_max_min_ = *p_v_at_00_max_min;
  const ViewDouble5D v_at_            = *p_v_at;
  const ViewDouble5D v_h_tu_swen_     = *p_v_h_tu_swen;
};

class FuncAdvTraTsp5 {
 public:
  FuncAdvTraTsp5 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    advection_tracer_tspas_5 (iblock, n_, k, j, i, v_ws_, v_at_, v_adv_tt_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void advection_tracer_tspas_5 (
      const int &iblock, const int &n, const int &k, const int &j, const int &i,
          const ViewDouble4D &v_www,
          const ViewDouble5D &v_ttt,
          const ViewDouble3D &v_adv_tt) const {

    double adv_xy1;
    if (v_at_00_max_min_(k, j, i  , 0) > v_at_00_max_min_(k, j, i  , 1) ||
        v_at_00_max_min_(k, j, i  , 0) < v_at_00_max_min_(k, j, i  , 2) ||
        v_at_00_max_min_(k, j, i+1, 0) > v_at_00_max_min_(k, j, i+1, 1) ||
        v_at_00_max_min_(k, j, i+1, 0) < v_at_00_max_min_(k, j, i+1, 2)) {
      adv_xy1 = - 
          (v_ttt(iblock, n, k, j, i+1) - v_ttt(iblock, n, k, j, i))
              * fabs(v_uv_ws_face_(k, j, i+1, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy1 = - dts_
        * (v_ttt(iblock, n, k, j, i+1) - v_ttt(iblock, n, k, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(k, j, i+1, 0) * v_uv_ws_face_(k, j, i+1, 0)
                / (v_h_tu_swen_(iblock, j, i+1, 0, 0) * v_h_tu_swen_(iblock, j, i+1, 0, 1));
                // / (v_htw_(iblock, j, i+1) * v_hun_(iblock, j, i+1));
    }
    double adv_xy2;
    if (v_at_00_max_min_(k, j, i  , 0) > v_at_00_max_min_(k, j, i  , 1) ||
        v_at_00_max_min_(k, j, i  , 0) < v_at_00_max_min_(k, j, i  , 2) ||
        v_at_00_max_min_(k, j, i-1, 0) > v_at_00_max_min_(k, j, i-1, 1) ||
        v_at_00_max_min_(k, j, i-1, 0) < v_at_00_max_min_(k, j, i-1, 2)) {
      adv_xy2 =   
          (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j, i-1))
              * fabs(v_uv_ws_face_(k, j, i, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy2 = dts_
        * (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j, i-1))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(k, j, i, 0) * v_uv_ws_face_(k, j, i, 0)
                / (v_h_tu_swen_(iblock, j, i, 0, 0) * v_h_tu_swen_(iblock, j, i, 0, 1));
                // / (v_htw_(iblock, j, i) * v_hun_(iblock, j, i));
    }

    double adv_xy3;
    if (v_at_00_max_min_(k, j  , i, 0) > v_at_00_max_min_(k, j  , i, 1) ||
        v_at_00_max_min_(k, j  , i, 0) < v_at_00_max_min_(k, j  , i, 2) ||
        v_at_00_max_min_(k, j+1, i, 0) > v_at_00_max_min_(k, j+1, i, 1) ||
        v_at_00_max_min_(k, j+1, i, 0) < v_at_00_max_min_(k, j+1, i, 2)) {
      adv_xy3 = -
          (v_ttt(iblock, n, k, j+1, i) - v_ttt(iblock, n, k, j, i))
              * fabs(v_uv_ws_face_(k, j, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy3 = - dts_
        * (v_ttt(iblock, n, k, j+1, i) - v_ttt(iblock, n, k, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(k, j, i, 1) * v_uv_ws_face_(k, j, i, 1)
                / (v_h_tu_swen_(iblock, j, i, 1, 0) * v_h_tu_swen_(iblock, j, i, 1, 1));
                // / (v_hts_(iblock, j, i) * v_hue_(iblock, j, i));
    }

    double adv_xy4;
    if (v_at_00_max_min_(k, j  , i, 0) > v_at_00_max_min_(k, j  , i, 1) ||
        v_at_00_max_min_(k, j  , i, 0) < v_at_00_max_min_(k, j  , i, 2) ||
        v_at_00_max_min_(k, j-1, i, 0) > v_at_00_max_min_(k, j-1, i, 1) ||
        v_at_00_max_min_(k, j-1, i, 0) < v_at_00_max_min_(k, j-1, i, 2)) {
      adv_xy4 =  
          (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j-1, i))
              * fabs(v_uv_ws_face_(k, j-1, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy4 = dts_
        * (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k, j-1, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(k, j-1, i, 1) * v_uv_ws_face_(k, j-1, i, 1)
                / (v_h_tu_swen_(iblock, j-1, i, 1, 0) * v_h_tu_swen_(iblock, j-1, i, 1, 1));
                // / (v_hts_(iblock, j-1, i) * v_hue_(iblock, j-1, i));
    }

    double adv_zb1;
    if (v_at_00_max_min_(k  , j, i, 0) > v_at_00_max_min_(k  , j, i, 1) ||
        v_at_00_max_min_(k  , j, i, 0) < v_at_00_max_min_(k  , j, i, 2) ||
        v_at_00_max_min_(k-1, j, i, 0) > v_at_00_max_min_(k-1, j, i, 1) ||
        v_at_00_max_min_(k-1, j, i, 0) < v_at_00_max_min_(k-1, j, i, 2)) {
      adv_zb1 = - 0.5 * fabs(v_www(iblock, k, j, i)) * v_odz_pt_(k, 0)
          * (v_ttt(iblock, n, k-1, j, i) - v_ttt(iblock, n, k, j, i));
    } else {
      adv_zb1 = - 0.5 * v_odz_pt_(k, 0) * v_www(iblock, k, j, i) * v_www(iblock, k, j, i) 
          * v_odz_pt_(k, 1) * (v_ttt(iblock, n, k-1, j, i) 
              - v_ttt(iblock, n, k, j, i)) * dts_;
    }

    double adv_zb2;
    if (v_at_00_max_min_(k  , j, i, 0) > v_at_00_max_min_(k  , j, i, 1) ||
        v_at_00_max_min_(k  , j, i, 0) < v_at_00_max_min_(k  , j, i, 2) ||
        v_at_00_max_min_(k+1, j, i, 0) > v_at_00_max_min_(k+1, j, i, 1) ||
        v_at_00_max_min_(k+1, j, i, 0) < v_at_00_max_min_(k+1, j, i, 2)) {
      adv_zb2 = 0.5 * fabs(v_www(iblock, k+1, j, i)) * v_odz_pt_(k, 0)
          * (v_ttt(iblock, n, k, j, i) - v_ttt(iblock, n, k+1, j, i));
    } else {
      adv_zb2 = 0.5 * v_odz_pt_(k, 0) * v_www(iblock, k+1, j, i) * v_www(iblock, k+1, j, i)
          * v_odz_pt_(k+1, 1) * (v_ttt(iblock, n, k, j, i) 
              - v_ttt(iblock, n, k+1, j, i)) * dts_;
    }

    const double adv_c1 = - v_ttt(iblock, n, k, j, i)
        * (v_uv_ws_face_(k, j, i+1, 0) - v_uv_ws_face_(k, j  , i, 0))
            * v_tarea_r_(iblock, j, i) * 2.0;
    const double adv_c2 = - v_ttt(iblock, n, k, j, i)
        * (v_uv_ws_face_(k, j, i  , 1) - v_uv_ws_face_(k, j-1, i, 1))
            * v_tarea_r_(iblock, j, i) * 2.0;

    const double adv_za = 
        0.5 * v_odz_pt_(k, 0) * v_www(iblock, k  , j, i)
            * (v_ttt(iblock, n, k, j, i) + v_ttt(iblock, n, k-1, j, i))
      - 0.5 * v_odz_pt_(k, 0) * v_www(iblock, k+1, j, i)
            * (v_ttt(iblock, n, k, j, i) + v_ttt(iblock, n, k+1, j, i));

    const double adv_zc = - v_odz_pt_(k, 0) * v_ttt(iblock, n, k, j, i)
        * (v_www(iblock, k, j, i) - v_www(iblock, k+1, j, i));

    const double adv_x0 =
        (v_ttt(iblock, n, k, j, i+1) + v_ttt(iblock, n, k, j, i  ))
            * v_uv_ws_face_(k, j  , i+1, 0) * v_tarea_r_(iblock, j, i) 
      - (v_ttt(iblock, n, k, j, i  ) + v_ttt(iblock, n, k, j, i-1))
            * v_uv_ws_face_(k, j  , i  , 0) * v_tarea_r_(iblock, j, i);
    const double adv_y0 = (
        (v_ttt(iblock, n, k, j+1, i) + v_ttt(iblock, n, k, j, i  ))
            * v_uv_ws_face_(k, j  , i  , 1) 
      - (v_ttt(iblock, n, k, j, i  ) + v_ttt(iblock, n, k, j-1, i))
            * v_uv_ws_face_(k, j-1, i  , 1)) * v_tarea_r_(iblock, j, i);

    const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
    const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
    const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

    v_adv_tt_(k, j, i) = adv_xx + adv_yy + adv_zz;
    return ;
  }

 private:
  const int n_;
  const double dts_ = CppPconstMod::dts;
  const ViewDouble2D v_odz_pt_        = *p_v_odz_pt;
  const ViewDouble3D v_adv_tt_        = *p_v_adv_tt;
  const ViewDouble3D v_tarea_r_       = *p_v_tarea_r;
  const ViewDouble4D v_ws_            = *p_v_ws;
  const ViewDouble4D v_uv_ws_face_    = *p_v_uv_ws_face;
  const ViewDouble4D v_at_00_max_min_ = *p_v_at_00_max_min;
  const ViewDouble5D v_at_            = *p_v_at;
  const ViewDouble5D v_h_tu_swen_     = *p_v_h_tu_swen;
};

class FuncAdvTraTsp6 {
 public:
  FuncAdvTraTsp6(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {

    const int iblock = 0;
    advection_tracer_tspas_6 (iblock, n_, j, i, v_ws_, v_at_, v_adv_tt_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void advection_tracer_tspas_6 (
      const int &iblock, const int &n, const int &j, const int &i,
          const ViewDouble4D &v_www,
          const ViewDouble5D &v_ttt,
          const ViewDouble3D &v_adv_tt) const {

    double adv_xy1;
    if (v_at_00_max_min_(KM-1, j, i  , 0) > v_at_00_max_min_(KM-1, j, i  , 1) ||
        v_at_00_max_min_(KM-1, j, i  , 0) < v_at_00_max_min_(KM-1, j, i  , 2) ||
        v_at_00_max_min_(KM-1, j, i+1, 0) > v_at_00_max_min_(KM-1, j, i+1, 1) ||
        v_at_00_max_min_(KM-1, j, i+1, 0) < v_at_00_max_min_(KM-1, j, i+1, 2)) {
      adv_xy1 = - 
          (v_ttt(iblock, n, KM-1, j, i+1) - v_ttt(iblock, n, KM-1, j, i))
              * fabs(v_uv_ws_face_(KM-1, j, i+1, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy1 = - dts_
        * (v_ttt(iblock, n, KM-1, j, i+1) - v_ttt(iblock, n, KM-1, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(KM-1, j, i+1, 0) * v_uv_ws_face_(KM-1, j, i+1, 0)
                / (v_h_tu_swen_(iblock, j, i+1, 0, 0) * v_h_tu_swen_(iblock, j, i+1, 0, 1));
                // / (v_htw_(iblock, j, i+1) * v_hun_(iblock, j, i+1));
    }
    double adv_xy2;
    if (v_at_00_max_min_(KM-1, j, i  , 0) > v_at_00_max_min_(KM-1, j, i  , 1) ||
        v_at_00_max_min_(KM-1, j, i  , 0) < v_at_00_max_min_(KM-1, j, i  , 2) ||
        v_at_00_max_min_(KM-1, j, i-1, 0) > v_at_00_max_min_(KM-1, j, i-1, 1) ||
        v_at_00_max_min_(KM-1, j, i-1, 0) < v_at_00_max_min_(KM-1, j, i-1, 2)) {
      adv_xy2 =   
          (v_ttt(iblock, n, KM-1, j, i) - v_ttt(iblock, n, KM-1, j, i-1))
              * fabs(v_uv_ws_face_(KM-1, j, i, 0)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy2 = dts_
        * (v_ttt(iblock, n, KM-1, j, i) - v_ttt(iblock, n, KM-1, j, i-1))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(KM-1, j, i, 0) * v_uv_ws_face_(KM-1, j, i, 0)
                / (v_h_tu_swen_(iblock, j, i, 0, 0) * v_h_tu_swen_(iblock, j, i, 0, 1));
                // / (v_htw_(iblock, j, i) * v_hun_(iblock, j, i));
    }

    double adv_xy3;
    if (v_at_00_max_min_(KM-1, j  , i, 0) > v_at_00_max_min_(KM-1, j  , i, 1) ||
        v_at_00_max_min_(KM-1, j  , i, 0) < v_at_00_max_min_(KM-1, j  , i, 2) ||
        v_at_00_max_min_(KM-1, j+1, i, 0) > v_at_00_max_min_(KM-1, j+1, i, 1) ||
        v_at_00_max_min_(KM-1, j+1, i, 0) < v_at_00_max_min_(KM-1, j+1, i, 2)) {
      adv_xy3 = -
          (v_ttt(iblock, n, KM-1, j+1, i) - v_ttt(iblock, n, KM-1, j, i))
              * fabs(v_uv_ws_face_(KM-1, j, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy3 = - dts_
        * (v_ttt(iblock, n, KM-1, j+1, i) - v_ttt(iblock, n, KM-1, j, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(KM-1, j, i, 1) * v_uv_ws_face_(KM-1, j, i, 1)
                / (v_h_tu_swen_(iblock, j, i, 1, 0) * v_h_tu_swen_(iblock, j, i, 1, 1));
                // / (v_hts_(iblock, j, i) * v_hue_(iblock, j, i));
    }

    double adv_xy4;
    if (v_at_00_max_min_(KM-1,   j, i, 0) > v_at_00_max_min_(KM-1, j  , i, 1) ||
        v_at_00_max_min_(KM-1,   j, i, 0) < v_at_00_max_min_(KM-1, j  , i, 2) ||
        v_at_00_max_min_(KM-1, j-1, i, 0) > v_at_00_max_min_(KM-1, j-1, i, 1) ||
        v_at_00_max_min_(KM-1, j-1, i, 0) < v_at_00_max_min_(KM-1, j-1, i, 2)) {
      adv_xy4 =  
          (v_ttt(iblock, n, KM-1, j, i) - v_ttt(iblock, n, KM-1, j-1, i))
              * fabs(v_uv_ws_face_(KM-1, j-1, i, 1)) * v_tarea_r_(iblock, j, i);
    } else {
      adv_xy4 = dts_
        * (v_ttt(iblock, n, KM-1, j, i) - v_ttt(iblock, n, KM-1, j-1, i))
            * 2.0 * v_tarea_r_(iblock, j, i) * v_uv_ws_face_(KM-1, j-1, i, 1) * v_uv_ws_face_(KM-1, j-1, i, 1)
                / (v_h_tu_swen_(iblock, j-1, i, 1, 0) * v_h_tu_swen_(iblock, j-1, i, 1, 1));
                // / (v_hts_(iblock, j-1, i) * v_hue_(iblock, j-1, i));
    }

    double adv_zb1;
    if (v_at_00_max_min_(KM-1, j, i, 0) > v_at_00_max_min_(KM-1, j, i, 1) ||
        v_at_00_max_min_(KM-1, j, i, 0) < v_at_00_max_min_(KM-1, j, i, 2) ||
        v_at_00_max_min_(KM-2, j, i, 0) > v_at_00_max_min_(KM-2, j, i, 1) ||
        v_at_00_max_min_(KM-2, j, i, 0) < v_at_00_max_min_(KM-2, j, i, 2)) {
      adv_zb1 = - 0.5 * fabs(v_ws_(iblock, KM-1, j, i)) * v_odz_pt_(KM-1, 0)
          * (v_ttt(iblock, n, KM-2, j, i) - v_ttt(iblock, n, KM-1, j, i));
    } else {
      adv_zb1 = - 0.5 * v_odz_pt_(KM-1, 0) * v_www(iblock, KM-1, j, i) * v_www(iblock, KM-1, j, i)
          * v_odz_pt_(KM-1, 1) * (v_ttt(iblock, n, KM-2, j, i) 
              - v_ttt(iblock, n, KM-1, j, i)) * dts_;
    }
    const double adv_zb2 = 0.0;

    const double adv_c1 = - v_ttt(iblock, n, KM-1, j, i)
        * (v_uv_ws_face_(KM-1, j  ,i+1, 0) - v_uv_ws_face_(KM-1, j  ,i  , 0))
            * v_tarea_r_(iblock, j, i) * 2.0;
    const double adv_c2 = - v_ttt(iblock, n, KM-1, j, i)
        * (v_uv_ws_face_(KM-1, j  ,i  , 1) - v_uv_ws_face_(KM-1, j-1,i  , 1))
            * v_tarea_r_(iblock, j, i) * 2.0;

    const double adv_za = 0.5 * v_odz_pt_(KM-1, 0) * v_ws_(iblock, KM-1, j, i)
        * (v_ttt(iblock, n, KM-1, j, i) + v_ttt(iblock, n, KM-2, j, i));

    const double adv_zc = - v_odz_pt_(KM-1, 0) * v_ttt(iblock, n, KM-1, j, i)
        * v_ws_(iblock, KM-1, j, i);


    const double adv_x0 =
        (v_ttt(iblock, n, KM-1, j  , i+1) + v_ttt(iblock, n, KM-1, j  , i  ))
            * v_uv_ws_face_(KM-1, j  ,i+1, 0) * v_tarea_r_(iblock, j, i) 
      - (v_ttt(iblock, n, KM-1, j  , i  ) + v_ttt(iblock, n, KM-1, j  , i-1))
            * v_uv_ws_face_(KM-1, j  ,i  , 0) * v_tarea_r_(iblock, j, i);
    const double adv_y0 = (
        (v_ttt(iblock, n, KM-1, j+1, i  ) + v_ttt(iblock, n, KM-1, j  , i  ))
            * v_uv_ws_face_(KM-1, j  ,i  , 1) 
      - (v_ttt(iblock, n, KM-1, j  , i  ) + v_ttt(iblock, n, KM-1, j-1, i  ))
            * v_uv_ws_face_(KM-1, j-1,i  , 1)) * v_tarea_r_(iblock, j, i);

    const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
    const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
    const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

    v_adv_tt_(KM-1, j, i) = adv_xx + adv_yy + adv_zz;
    return ;
  }

 private:
  const int n_;
  const double dts_ = CppPconstMod::dts;
  const ViewDouble2D v_odz_pt_        = *p_v_odz_pt;
  const ViewDouble3D v_adv_tt_        = *p_v_adv_tt;
  const ViewDouble3D v_tarea_r_       = *p_v_tarea_r;
  const ViewDouble4D v_ws_            = *p_v_ws;
  const ViewDouble4D v_uv_ws_face_    = *p_v_uv_ws_face;
  const ViewDouble4D v_at_00_max_min_ = *p_v_at_00_max_min;
  const ViewDouble5D v_at_            = *p_v_at;
  const ViewDouble5D v_h_tu_swen_     = *p_v_h_tu_swen;
};

class FunctorTracer15 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_tf_(iblock, k, j, i) = v_adv_tt_(k, j, i)
        * v_vit_(iblock, k, j, i);
    return ;
  }
 private:
  const ViewDouble3D v_adv_tt_ = *p_v_adv_tt;
  const ViewDouble4D v_tf_     = *p_v_tf;
  const ViewDouble4D v_vit_    = *p_v_vit;
};

// #ifdef CANUTO
class FunctorTracer16 {
 public:
  FunctorTracer16 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    if (v_akt_(iblock, n_, 0, j, i) < dwndmix_) {
      v_akt_(iblock, n_, 0, j, i) = dwndmix_;
    }
#ifdef ISO
    v_wkc_(iblock, k, j, i) = v_akt_(iblock, n_, k, j, i)
        + v_ahisop_(iblock, j, i) * v_k3_(iblock, 2, j, k+1, i);
#else  // ISO
    v_wkc_(iblock, k, j, i) = v_akt_(iblock, n_, k, j, i);
#endif // ISO
    return ;
  }
 private:
  const int n_;
  const double dwndmix_ = CppPconstMod::dwndmix;

  const ViewDouble4D v_wkc_ = *p_v_wkc;
  const ViewDouble5D v_akt_ = *p_v_akt;

#ifdef ISO
  const ViewDouble3D v_ahisop_ = *p_v_ahisop;
  const ViewDouble5D v_k3_     = *p_v_k3;
#endif // ISO
};
// #endif // CANUTO
class FunctorTracer17{
 public:
  FunctorTracer17 (const int &n) : n_(n) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    hdifft_del4_1 (n_, k, j, i, v_dt2k_, v_atb_);
#endif // BIHAR
    return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdifft_del4_1 (const int &n, 
      const int &k, const int &j, const int &i,
          const ViewDouble3D &v_d2tk, const ViewDouble5D &v_tmix) const {
    const int bid = 0;
    // n s e w c
    int kk = k + 1;
    v_c_cnsew_(k, j, i, 1) = (kk <= v_kmt_nsew_(bid, j, i, 0) && kk <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 0) : C0;                
    v_c_cnsew_(k, j, i, 2) = (kk <= v_kmt_nsew_(bid, j, i, 1) && kk <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 1) : C0;                
    v_c_cnsew_(k, j, i, 3) = (kk <= v_kmt_nsew_(bid, j, i, 2) && kk <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 2) : C0;                
    v_c_cnsew_(k, j, i, 4) = (kk <= v_kmt_nsew_(bid, j, i, 3) && kk <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 3) : C0;

    v_c_cnsew_(k, j, i, 0) = -(v_c_cnsew_(k, j, i, 1) + v_c_cnsew_(k, j, i, 2) 
                             + v_c_cnsew_(k, j, i, 3) + v_c_cnsew_(k, j, i, 4));

    if (i >= (ib_ - 2) && i < (ie_ + 1) && j >= (jb_ - 2) && j < (je_ + 1)) {
      v_d2tk(k, j, i) = v_ahf_(bid, j, i) * 
          (v_c_cnsew_(k, j, i, 0) * v_tmix(bid, n, k+1, j  , i  ) 
         + v_c_cnsew_(k, j, i, 1) * v_tmix(bid, n, k+1, j-1, i  ) 
         + v_c_cnsew_(k, j, i, 2) * v_tmix(bid, n, k+1, j+1, i  ) 
         + v_c_cnsew_(k, j, i, 3) * v_tmix(bid, n, k+1, j  , i+1) 
         + v_c_cnsew_(k, j, i, 4) * v_tmix(bid, n, k+1, j  , i-1));
    }
    return ;
  }
#endif // BIHAR
 private:
  const int n_;
#ifdef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const ViewInt3D v_kmt_        = *p_v_kmt;
  const ViewInt4D v_kmt_nsew_   = *p_v_kmt_nsew;
  const ViewDouble3D v_dt2k_    = *p_v_dt2k;
  const ViewDouble3D v_ahf_     = *p_v_ahf;
  const ViewDouble4D v_c_cnsew_ = *p_v_c_cnsew;
  const ViewDouble4D v_dt_nsew_ = *p_v_dt_nsew;
  const ViewDouble5D v_atb_     = *p_v_atb;
#endif // BIHAR
};
class FunctorTracer18{
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
#ifdef BIHAR
    const int iblock = 0;
    double hdtk;
    hdifft_del4_2(k, j, i, v_dt2k_, hdtk);
    v_tf_(iblock, k, j, i) += hdtk;
#endif // BIHAR
    return ;
  }
#ifdef BIHAR
  KOKKOS_INLINE_FUNCTION void hdifft_del4_2 (const int &k, const int &j, 
      const int &i, const ViewDouble3D &v_d2tk, double &hdtk) const {
    hdtk = C0;
    if (i >= (ib_ - 1) && i < ie_ && j >= (jb_ - 1) && j < je_) {
      hdtk = ah_ * (v_c_cnsew_(k, j, i, 0) * v_d2tk(k, j  , i  ) 
                  + v_c_cnsew_(k, j, i, 1) * v_d2tk(k, j-1, i  ) 
                  + v_c_cnsew_(k, j, i, 2) * v_d2tk(k, j+1, i  ) 
                  + v_c_cnsew_(k, j, i, 3) * v_d2tk(k, j  , i+1) 
                  + v_c_cnsew_(k, j, i, 4) * v_d2tk(k, j  , i-1));
    }
    return;
  }
 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double ah_ = CppHmixDel4::ah;
  const ViewDouble3D v_dt2k_    = *p_v_dt2k;
  const ViewDouble4D v_c_cnsew_ = *p_v_c_cnsew;
  const ViewDouble4D v_tf_      = *p_v_tf;
#endif // BIHAR
};

class FunctorTracer17Del2 {
 public:
  FunctorTracer17Del2 (const int &n) : n_(n) {}
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
#ifndef BIHAR
    const int iblock = 0;
    double hdtk = 0.0;
    hdifft_del2 (n_, k, j, i, hdtk, v_atb_);
    v_tf_(iblock, k, j, i) += hdtk;
#endif // BIHAR
    return;
  }
#ifndef BIHAR
  KOKKOS_INLINE_FUNCTION 
  void hdifft_del2 (const int &n, const int &k, const int &j,const int &i, 
      double &hdtk, const ViewDouble5D &v_tmix) const {
    const int bid = 0;
    const int kk = k + 1;
    const double cn = (kk <= v_kmt_nsew_(bid, j, i, 0) && kk <= v_kmt_(bid, j, i))
        ? v_dtn_(bid, j, i) : C0;
    const double cs = (kk <= v_kmt_nsew_(bid, j, i, 1) && kk <= v_kmt_(bid, j, i))
        ? v_dts_(bid, j, i) : C0;
    const double ce = (kk <= v_kmt_nsew_(bid, j, i, 2) && kk <= v_kmt_(bid, j, i))
        ? v_dte_(bid, j, i) : C0;
    const double cw = (kk <= v_kmt_nsew_(bid, j, i, 3) && kk <= v_kmt_(bid, j, i))
        ? v_dtw_(bid, j, i) : C0;

    const double cc = - (cn + cs + ce + cw);

    if (i >= (ib_ - 1) && i < ie_ && j >= (jb_ - 1) && j < je_) {
      hdtk = ah_ * (cc * v_tmix(bid, n, kk, j    , i    ) 
                  + cn * v_tmix(bid, n, kk, j - 1, i    ) 
                  + cs * v_tmix(bid, n, kk, j + 1, i    ) 
                  + ce * v_tmix(bid, n, kk, j    , i + 1) 
                  + cw * v_tmix(bid, n, kk, j    , i - 1));
    }
    return;
  }
#endif // BIHAR

 private:
  const int n_;
#ifndef BIHAR
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ie;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double ah_ = CppHmixDel2::ah;
  const ViewInt3D v_kmt_      = *p_v_kmt;
  const ViewInt4D v_kmt_nsew_ = *p_v_kmt_nsew;
  const ViewDouble3D v_dtn_   = *p_v_dtn;
  const ViewDouble3D v_dts_   = *p_v_dts;
  const ViewDouble3D v_dte_   = *p_v_dte;
  const ViewDouble3D v_dtw_   = *p_v_dtw;
  const ViewDouble4D v_tf_    = *p_v_tf;
  const ViewDouble5D v_atb_   = *p_v_atb;
#endif // BIHAR
};

// !----------------------------------------
// !     VERTICAL COMPONENT
// !----------------------------------------
class FunctorTracer19 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
#if defined(SOLAR) || defined(SOLARCHLORO)
    const int iblock = 0;
    double wt;
#endif // defined(SOLAR) || defined(SOLARCHLORO)
#ifdef SOLAR
    wt = v_swv_(iblock, j, i) * v_pen_ (0) * v_vit_(iblock, 1, j, i);
#endif // SOLAR
#ifdef SOLARCHLORO
    wt = v_swv_(iblock, j, i) * v_pen_chl_(iblock, 0, j, i) * v_vit_(iblock, 1, j, i);
#endif // SOLARCHLORO
#if defined(SOLAR) || defined(SOLARCHLORO)
    v_tf_ (iblock, 0, j, i) -= v_odzp_(0) * wt;
#endif // defined(SOLAR) || defined(SOLARCHLORO)
    return ;
  }
 private:
  const ViewDouble1D v_odzp_ = *p_v_odzp;
#ifdef SOLAR
  const ViewDouble1D v_pen_  = *p_v_pen;
#endif // SOLAR
  const ViewDouble3D v_swv_  = *p_v_swv;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
#ifdef SOLARCHLORO
  const ViewDouble4D v_pen_chl_ = *p_v_pen_chl;
#endif // SOLARCHLORO
};
class FunctorTracer20 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {

#if defined(SOLAR) || defined(SOLARCHLORO)
    const int iblock = 0;
    double wt1(0.0), wt2(0.0), tmp(0.0);
#endif // defined(SOLAR) || defined(SOLARCHLORO)

#ifdef SOLAR
    wt1 = v_swv_(iblock, j, i) * v_pen_(k-1) * v_vit_(iblock, k  , j, i);
    wt2 = v_swv_(iblock, j, i) * v_pen_(k  ) * v_vit_(iblock, k+1, j, i);
#endif // SOLAR

#ifdef SOLARCHLORO
    wt1 = v_swv_(iblock, j, i) * v_pen_chl_(iblock, k-1, j, i) * v_vit_(iblock, k  , j, i);
    wt2 = v_swv_(iblock, j, i) * v_pen_chl_(iblock, k  , j, i) * v_vit_(iblock, k+1, j, i);
#endif // SOLARCHLORO

#if defined(SOLAR) || defined(SOLARCHLORO)
    tmp = wt1 - wt2;
    v_tf_(iblock, k, j, i) += v_odzp_(k) * tmp;
#endif // defined(SOLAR) || defined(SOLARCHLORO)
    return ;
  }

 private:
  const ViewDouble1D v_odzp_ = *p_v_odzp;
#ifdef SOLAR
  const ViewDouble1D v_pen_  = *p_v_pen;
#endif // SOLAR
  const ViewDouble3D v_swv_  = *p_v_swv;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
#ifdef SOLARCHLORO
  const ViewDouble4D v_pen_chl_ = *p_v_pen_chl;
#endif // SOLARCHLORO
};
class FunctorTracer21 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {

#if defined(SOLAR) || defined(SOLARCHLORO)
    const int iblock = 0;
    double wt(0.0);
#endif // defined(SOLAR) || defined(SOLARCHLORO)

#ifdef SOLAR
    wt = v_swv_(iblock, j, i) * v_pen_(KM-2) * v_vit_(iblock, KM-1, j, i);
#endif // SOLAR

#ifdef SOLARCHLORO
    wt = v_swv_(iblock, j, i) * v_pen_chl_(iblock, KM-2, j, i) * v_vit_(iblock, KM-1, j, i);
#endif // SOLARCHLORO

#if defined(SOLAR) || defined(SOLARCHLORO)
    v_tf_(iblock, KM-1, j, i) += v_odzp_(KM-1) * wt;
#endif // defined(SOLAR) || defined(SOLARCHLORO)
    return ;
  }
 private:
  const ViewDouble1D v_odzp_ = *p_v_odzp;
#ifdef SOLAR
  const ViewDouble1D v_pen_  = *p_v_pen;
#endif // SOLAR
  const ViewDouble3D v_swv_  = *p_v_swv;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
#ifdef SOLARCHLORO
  const ViewDouble4D v_pen_chl_ = *p_v_pen_chl;
#endif // SOLARCHLORO
};

class FunctorTracer22 {
 public:
  FunctorTracer22 (const int &n) : n_(n) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    // const double aidif = 0.5;
    const double one_minus_aidif = 0.5;

    const double wt1 = v_wkc_(iblock, 0, j, i)
        * (v_atb_(iblock, n_, 1, j, i) - v_atb_(iblock, n_, 2, j, i))
            * v_odzt_(1) * v_vit_(iblock, 1, j, i);

    v_tf_(iblock, 0, j, i) -= v_odzp_(0) * wt1 * one_minus_aidif;
    return ;
  }

 private:
  const int n_;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wkc_  = *p_v_wkc;
  const ViewDouble5D v_atb_  = *p_v_atb;
};
class FunctorTracer23 {
 public:
  FunctorTracer23 (const int &n) : n_(n) {}
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    // const double aidif = 0.5;
    const double one_minus_aidif = 0.5;
    const double wt1 = v_wkc_(iblock, k-1, j, i)
        * (v_atb_(iblock, n_, k  , j, i) - v_atb_(iblock, n_, k+1, j, i))
            * v_odzt_(k  ) * v_vit_(iblock, k  , j, i);
    const double wt2 = v_wkc_(iblock, k  , j, i)
        * (v_atb_(iblock, n_, k+1, j, i) - v_atb_(iblock, n_, k+2, j, i))
            * v_odzt_(k+1) * v_vit_(iblock, k+1, j, i);
   
    v_tf_(iblock, k, j, i) += v_odzp_(k) * (wt1 - wt2) * one_minus_aidif;
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wkc_  = *p_v_wkc;
  const ViewDouble5D v_atb_  = *p_v_atb;
};
class FunctorTracer24 {
 public:
  FunctorTracer24 (const int &n) : n_(n) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    // const double aidif = 0.5;
    const double one_minus_aidif = 0.5;
    const double wt2 = v_wkc_(iblock, KM-2, j, i)
        * (v_atb_(iblock, n_, KM-1, j, i) - v_atb_(iblock, n_, KM, j, i))
            * v_odzt_(KM-1) * v_vit_(iblock, KM-1, j, i);
    v_tf_(iblock, KM-1, j, i) += v_odzp_(KM-1) * wt2 * one_minus_aidif;
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_wkc_  = *p_v_wkc;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class FunctorTracer25 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    // const double aidif = 0.5;
    const double one_minus_aidif = 0.5;
    if (v_kmt_(iblock, j, i) > 0) {
#ifdef COUP
      if (boundary_restore_ == 2) {
        v_stf_(iblock, j, i) = v_ssf_(iblock, j, i)
            / v_odzp_(0) + (gamma_ * 50.0 * v_odzp_(0))
                * (v_sss_(iblock, j, i) - v_atb_(iblock, 1, 1, j, i))
                    * v_ifrac_(iblock, j, i) 
            / v_odzp_(0) + (gamma_ * 30.0 / 365.0 / 4.0 * 50.0 * v_odzp_(0))
                * (v_sss_(iblock, j, i) - v_atb_(iblock, 1, 1, j, i))
                    / v_odzp_(0);
      } else {
        v_stf_(iblock, j, i) = v_ssf_(iblock, j, i) / v_odzp_(0);
      }
#else  // COUP
#ifdef FRC_CORE
      v_stf_(iblock, j, i) = (v_fresh_(iblock, j, i) * 34.7 * od0_ * 1.0e-3
          + gamma_ 
              * (v_sss_(iblock, j, i) - v_atb_(iblock, 1, 1, j, i))
                  * v_seaice_(iblock, j, i) / v_odzp_(0)
          + gamma_ * 30.0 / 365.0 / 4.0 * 50.0
              * (v_sss_(iblock, j, i) - v_atb_(iblock, 1, 1, j, i))
                  / v_odzp_(0) * (1.0 - v_seaice_(iblock, j, i)));
#else  // FRC_CORE
      v_stf_(iblock, j, i) = gamma_ * (v_sss_(iblock, j, i) 
          - v_atb_(iblock, 1, 1, j, i)) / v_odzp_(0);
#endif // FRC_CORE
#endif // COUP
      v_tf_(iblock, 0, j, i) += v_stf_(iblock, j, i)
          * one_minus_aidif * v_odzp_(0);

      v_net_(iblock, 1, j, i) = v_stf_(iblock, j, i) * v_odzp_(0);
    }
    return ;
  }
 private:
  const double od0_   = CppPconstMod::od0;
  const double gamma_ = CppPconstMod::gamma;

  const ViewInt3D    v_kmt_  = *p_v_kmt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_sss_  = *p_v_sss;
  const ViewDouble3D v_stf_  = *p_v_stf;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_net_  = *p_v_net;
  const ViewDouble5D v_atb_  = *p_v_atb;

#ifdef COUP                  
  const ViewDouble3D v_ssf_   = *p_v_ssf;
  const ViewDouble3D v_ifrac_ = *p_v_ifrac;
#else  // COUP
#ifdef FRC_CORE
  const ViewDouble3D v_fresh_  = *p_v_fresh;
  const ViewDouble3D v_seaice_ = *p_v_seaice;
#endif // FRC_CORE
#endif // COUP
};

#ifdef SSSNORM
class FunctorTracer26 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i, double &err_norm2) const {
    const int iblock = 0;

    err_norm2 += v_tarea_(iblock, j, i) * v_net_(iblock, 1, j, i)
        * v_vit_(iblock, 0, j, i); 

    return ;
  }
 private:
  const ViewDouble3D v_tarea_ = *p_v_tarea;
  const ViewDouble4D v_vit_   = *p_v_vit;
  const ViewDouble4D v_net_   = *p_v_net;
};

class FunctorTracer27 {
 public:
  FunctorTracer27 (const double &err_norm2) : err_norm2_(err_norm2) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_tf_(iblock, 0, j, i) += err_norm2_ * v_vit_(iblock, 0, j, i);
    return ;
  }
 private:
  const double err_norm2_;

  const ViewDouble4D v_tf_  = *p_v_tf;
  const ViewDouble4D v_vit_ = *p_v_vit;
};
#endif // SSSNORM

class FunctorTracer28 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    if (v_kmt_(iblock, j, i) > 0) {
#ifdef COUP
      v_stf_(iblock, j, i) = v_tsf_(iblock, j, i);
#else  // COUP
#ifdef FRC_CORE
      v_stf_(iblock, j, i) = ((v_swv_(iblock, j, i) + v_nswv_(iblock, j, i))
          * od0cp_ + v_seaice_(iblock, j, i) * gamma_
              * (v_sst_(iblock, j, i) - v_atb_(iblock, 0, 1, j, i)) / v_odzp_(0));
#else  // FRC_CORE
      v_stf_(iblock, j, i) = (v_swv_(iblock, j, i) + v_nswv_(iblock, j, i)
          - v_dqdt_(iblock, j, i) * (v_sst_(iblock, j, i)
              - v_atb_(iblock, 0, 1, j, i))) * od0cp_;
#endif // FRC_CORE
#endif // COUP
    }
    return ;
  }
 private:
  const double od0cp_ = CppPconstMod::od0cp;
  const double gamma_ = CppPconstMod::gamma;

  const ViewInt3D    v_kmt_    = *p_v_kmt;
  const ViewDouble1D v_odzp_   = *p_v_odzp;
  const ViewDouble3D v_stf_    = *p_v_stf;

#ifdef COUP
  const ViewDouble3D v_tsf_    = *p_v_tsf;
#else  // COUP
#ifdef FRC_CORE
  const ViewDouble3D v_nswv_   = *p_v_nswv;
  const ViewDouble3D v_seaice_ = *p_v_seaice;
#else  // FRC_CORE
  const ViewDouble3D v_dqdt_   = *p_v_dqdt;
#endif // FRC_CORE
  const ViewDouble3D v_sst_    = *p_v_sst;
  const ViewDouble3D v_swv_    = *p_v_swv;
  const ViewDouble5D v_atb_    = *p_v_atb;

#endif // COUP
};

class FunctorTracer29 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    // const double aidif = 0.5;
    const double one_minus_aidif = 0.5;
    if (v_kmt_(iblock, j, i) > 0) {
      v_tf_(iblock, 0, j, i) += v_stf_(iblock, j, i) * v_odzp_(0) * one_minus_aidif;

      v_net_(iblock, 0, j, i) = v_stf_(iblock, j, i) * v_odzp_(0);
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmt_  = *p_v_kmt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble3D v_stf_  = *p_v_stf;
  const ViewDouble4D v_tf_   = *p_v_tf;
  const ViewDouble4D v_net_  = *p_v_net;
};

class FunctorTracer30 {
 public:
  FunctorTracer30 (const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_vtl_(iblock, k, j, i) = v_at_(iblock, n_, k, j, i)
        + dts_ * v_tf_(iblock, k, j, i);
    return ;
  }
 private:

  const int n_;

  const double dts_ = CppPconstMod::dts;

  const ViewDouble4D v_tf_  = *p_v_tf;
  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
};

// class FunctorTracer31 {
//  public:
//   KOKKOS_INLINE_FUNCTION void operator () (
//       const int &k, const int &j, const int &i) const {
//     const int iblock = 0;
//     v_vtl_ori_(iblock, k, j, i) = v_vtl_(iblock, k, j, i);

//     return ;
//   }
//  private:
//   const ViewDouble4D v_vtl_     = *p_v_vtl;
//   const ViewDouble4D v_vtl_ori_ = *p_v_vtl_ori;
// };

// invtrit(vtl, stf, wkc, aidif, c2dtts)
class FunctorTracer32 {
 public:
  FunctorTracer32 (const double &c2dtts) : c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const double aidif = 0.5;

    invtrit (j, i, v_vtl_, v_stf_, v_wkc_, aidif, c2dtts_);

    return ;
  }
KOKKOS_INLINE_FUNCTION void invtrit (const int &j, const int &i,
      const ViewDouble4D &v_wk, const ViewDouble3D &v_topbc,
          const ViewDouble4D &v_dcb, const double &aidif, const double &c2dtts) const {

    const int iblock = 0;

    // const double c2dtts_times_aidif = c2dtts * aidif;

    if (v_kmt_(iblock, j, i) > 0) {
      const int kz = v_kmt_(iblock, j, i) - 1;

      for (int k = 1; k <= kz; ++k) {
        v_a8_(j, i, k) = v_dcb(iblock, k-1, j, i) * v_odzt_(k  ) * v_odzp_(k)
            // * c2dtts_times_aidif;
            * c2dtts * aidif;

        v_d8_(j, i, k) = v_wk(iblock, k, j, i);
      }
      for (int k = 1; k <= kz-1; ++k) {
        v_c8_(j, i, k) = v_dcb(iblock, k  , j, i) * v_odzt_(k+1) * v_odzp_(k)
            // * c2dtts_times_aidif;
            * c2dtts * aidif;

        v_b8_(j, i, k) = 1.0 + v_a8_(j, i, k) + v_c8_(j, i, k);
        v_e8_(j, i, k) = 0.0;
        v_f8_(j, i, k) = 0.0;
      }
      // k = 0
      // a8[0] = v_odzp_(0) * c2dtts_times_aidif;
      v_a8_(j, i, 0) = v_odzp_(0) * c2dtts * aidif;
      v_c8_(j, i, 0) = v_dcb(iblock, 0, j, i) * v_odzt_(1) * v_odzp_(0)
            // * c2dtts_times_aidif;
            * c2dtts * aidif;
      v_b8_(j, i, 0) = 1.0 + v_c8_(j, i, 0);
      v_d8_(j, i, 0) = v_wk(iblock, 0, j, i);
      v_e8_(j, i, 0) = 0.0;
      v_f8_(j, i, 0) = 0.0;

      v_b8_(j, i, kz) = 1.0 + v_a8_(j, i, kz);
      // c8[kz] = v_odzp_(kz) * c2dtts_times_aidif;
      v_c8_(j, i, kz) = v_odzp_(kz) * c2dtts * aidif;

      v_e8_(j, i, kz+1) = 0.0;
      v_f8_(j, i, kz+1) = 0.0;

      for (int k = kz; k >= 0; --k) {
        const double g0 = 1.0 / (v_b8_(j, i, k) - v_c8_(j, i, k) * v_e8_(j, i, k+1));

        v_e8_(j, i, k) = v_a8_(j, i, k) * g0;
        v_f8_(j, i, k) = (v_d8_(j, i, k) + v_c8_(j, i, k) * v_f8_(j, i, k+1)) * g0;
      }

      double wk;
      wk = (v_e8_(j, i, 0) * v_topbc(iblock, j, i) + v_f8_(j, i, 0))
          * v_vit_(iblock, 0, j, i);
      v_wk(iblock, 0, j, i) = wk;
      for (int k = 1; k <= kz; ++k) {
        wk = (v_e8_(j, i, k) * wk + v_f8_(j, i, k)) * v_vit_(iblock, k, j, i);
        v_wk(iblock, k, j, i)= wk;
      }
    }
    return ;
  }
 private:
  const double c2dtts_;
  const ViewInt3D    v_kmt_  = *p_v_kmt;
  const ViewDouble1D v_odzp_ = *p_v_odzp;
  const ViewDouble1D v_odzt_ = *p_v_odzt;
  const ViewDouble3D v_a8_   = *p_v_wp1;
  const ViewDouble3D v_b8_   = *p_v_wp2;
  const ViewDouble3D v_c8_   = *p_v_wp3;
  const ViewDouble3D v_d8_   = *p_v_wp4;
  const ViewDouble3D v_e8_   = *p_v_wk1;
  const ViewDouble3D v_f8_   = *p_v_wk2;
  const ViewDouble3D v_stf_  = *p_v_stf;
  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_vtl_  = *p_v_vtl;
  const ViewDouble4D v_wkc_  = *p_v_wkc;
};
// End invtrit

class FunctorTracer33 {
 public:
  FunctorTracer33 (const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int k(0), iblock(0);
    const double aidif = 0.5;
    v_dt_diff_(iblock, n_, k, j, i) =
        (v_vtl_(iblock, k, j, i) - v_vtl_ori_(iblock, k, j, i))
            / c2dtts_ * v_vit_(iblock, k, j, i)
                - v_stf_(iblock, j, i) * aidif * v_odzp_(k);
    return ;
  }

 private:
  const int n_;
  const double c2dtts_;
  const ViewDouble1D v_odzp_    = *p_v_odzp;
  const ViewDouble3D v_stf_     = *p_v_stf;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble4D v_vtl_     = *p_v_vtl;
  const ViewDouble4D v_vtl_ori_ = *p_v_vtl_ori;
  const ViewDouble5D v_dt_diff_ = *p_v_dt_diff;
};

class FunctorTracer34 {
 public:
  FunctorTracer34 (const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_dt_diff_(iblock, n_, k, j, i) =
        (v_vtl_(iblock, k, j, i) - v_vtl_ori_(iblock, k, j, i))
            / c2dtts_ * v_vit_(iblock, k, j, i);
    return ;
  }

 private:
  const int n_;
  const double c2dtts_;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble4D v_vtl_     = *p_v_vtl;
  const ViewDouble4D v_vtl_ori_ = *p_v_vtl_ori;
  const ViewDouble5D v_dt_diff_ = *p_v_dt_diff;
};

// !-----------------------------------------------------------------------
// !     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
// !-----------------------------------------------------------------------
class FunctorTracer35 {
 public:
  FunctorTracer35 (const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    // v_tend_(iblock, n_, k, j, i) = (v_vtl_(iblock, k, j, i)
    //     - v_at_(iblock, n_, k, j, i)) / c2dtts_
    //         * v_vit_(iblock, k, j, i);

    v_at_ (iblock, n_, k  , j, i) = v_vtl_(iblock, k, j, i);

    v_atb_(iblock, n_, k+1, j, i) = v_vtl_(iblock, k, j, i);
    return ;
  }
 private:
  const int n_;
  const double c2dtts_;

  // const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_vtl_  = *p_v_vtl;
  const ViewDouble5D v_at_   = *p_v_at;
  const ViewDouble5D v_atb_  = *p_v_atb;
  // const ViewDouble5D v_tend_ = *p_v_tend;
};

//==================================

// class functor_tracer_10 {
//  public:
//   functor_tracer_10(const int &n, const int &k, const int &iblock) 
//     : n_(n), k_(k), iblock_(iblock) {}

//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &j, const int &i) const {
//     v_tf_(i, j, k_, iblock_)    += v_hdtk_(j, i);
//     v_dx_(i, j, k_, n_, iblock_) = v_hdtk_(j, i);
//     return ;
//   }
//  private:
//   const int n_, k_, iblock_;

//   const ViewDouble2D v_hdtk_ = *p_v_hdtk;
//   const ViewDouble4D v_tf_   = *p_v_tf;
//   const ViewDouble5D v_dx_   = *p_v_dx;
// };


class functor_tracer_21_22 {
 public:
  functor_tracer_21_22(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    const double lamda = 1.0 / (15.0 * 86400.0);

    v_tf_(iblock, k, j, i) += v_vit_(iblock, k, j, i)
        * (v_restore_(iblock, n_, k, j, i) - v_atb_(iblock, n_, k+1, j, i))
            * lamda;
    return ;
  }
 private:
  const int n_;

  const ViewDouble4D v_tf_      = *p_v_tf;
  const ViewDouble4D v_vit_     = *p_v_vit;
  const ViewDouble5D v_restore_ = *p_v_restore;
  const ViewDouble5D v_atb_     = *p_v_atb;
};

 
class functor_tracer_24 {
 public:
  functor_tracer_24(const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_vtl_(iblock, k, j, i) = v_at_(iblock, n_, k, j, i)
        + c2dtts_ * v_tf_(iblock, k, j, i);
    return ;
  }
 private:
  const int n_;
  const double c2dtts_;

  const ViewDouble4D v_tf_  = *p_v_tf;
  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
};



class functor_tracer_27 {
 public:
  functor_tracer_27(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    v_vtl_(iblock, 0, j, i) = v_vtl_(iblock, 0, j, i) 
        - v_at_(iblock, n_, 0, j, i) - v_net_(i, j, n_, iblock) * dts_;

    return ;
  }

 private:
  const int n_;

  const double dts_ = CppPconstMod::dts;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_net_ = *p_v_net;
  const ViewDouble5D v_at_  = *p_v_at;
};

class functor_tracer_28 {
 public:
  functor_tracer_28(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_vtl_(iblock, k, j, i) -= v_at_(iblock, n_, k, j, i);

    return ;
  }

 private:
  const int n_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
};

class functor_tracer_29 {
 public:
  functor_tracer_29(const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    v_vtl_(iblock, 0, j, i) = v_vtl_(iblock, 0, j, i) 
        - v_atb_(iblock, n_, 1, j, i) - v_net_(i, j, n_, iblock) * c2dtts_;

    return ;
  }

 private:
  const int n_;
  const double c2dtts_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_net_ = *p_v_net;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_30 {
 public:
  functor_tracer_30(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_vtl_(iblock, k, j, i) -= v_atb_(iblock, n_, k+1, j, i);

    return ;
  }

 private:
  const int n_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_31 {
 public:
  functor_tracer_31(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    v_vtl_(iblock, 0, j, i) += v_at_(iblock, n_, 0, j, i) 
        + v_net_(i, j, n_, iblock) * dts_;

    return ;
  }

 private:
  const int n_;

  const double dts_ = CppPconstMod::dts;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_net_ = *p_v_net;
  const ViewDouble5D v_at_  = *p_v_at;
};

class functor_tracer_32 {
 public:
  functor_tracer_32(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_vtl_(iblock, k, j, i) += v_at_(iblock, n_, k, j, i);
    return ;
  }

 private:
  const int n_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
};

class functor_tracer_33 {
 public:
  functor_tracer_33(const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    v_vtl_(iblock, 0, j, i) += v_atb_(iblock, n_, 1, j, i)
        + v_net_(i, j, n_, iblock) * c2dtts_;

    return ;
  }

 private:
  const int n_;
  const double c2dtts_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble4D v_net_ = *p_v_net;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_34 {
 public:
  functor_tracer_34(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_vtl_(iblock, k, j, i) += v_atb_(iblock, n_, k+1, j, i);

    return ;
  }

 private:
  const int n_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_35 {
 public:
  functor_tracer_35(const int &n, const double &c2dtts) 
    : n_(n), c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_tend_(iblock, n_, k, j, i) = (v_vtl_(iblock, k, j, i)
        - v_at_(iblock, n_, k, j, i)) / c2dtts_
            * v_vit_(iblock, k, j, i);
    return ;
  }

 private:
  const int n_;
  const double c2dtts_;

  const ViewDouble4D v_vit_  = *p_v_vit;
  const ViewDouble4D v_vtl_  = *p_v_vtl;
  const ViewDouble5D v_at_   = *p_v_at;
  const ViewDouble5D v_tend_ = *p_v_tend;
};

class functor_tracer_37 {
 public:
  functor_tracer_37(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_atb_(iblock, n_, k+1, j, i) = aft2_ * v_at_(iblock, n_, k, j, i)
        + aft1_ * (v_atb_(iblock, n_, k+1, j, i) + v_vtl_(iblock, k, j, i));

    return ;
  }

 private:
  const int n_;
  const double aft1_ = CppPconstMod::aft1;
  const double aft2_ = CppPconstMod::aft2;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_38 {
 public:
  functor_tracer_38(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_at_(iblock, n_, k, j, i) = v_vtl_(iblock, k, j, i);

    return ;
  }

 private:
  const int n_;

  const ViewDouble4D v_vtl_ = *p_v_vtl;
  const ViewDouble5D v_at_  = *p_v_at;
};

#ifdef ISO
// isoflux(n)
class functor_tracer_isoflux_1 {
 public:
  functor_tracer_isoflux_1(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_temp_(iblock, k, j, i) = P25 * v_dzr_(k) * (
        (v_atb_(iblock, n_, k, j+1, i) - v_atb_(iblock, n_, k+2, j+1, i))
      + (v_atb_(iblock, n_, k, j  , i) - v_atb_(iblock, n_, k+2, j  , i)));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_2 {
 public:
  functor_tracer_isoflux_2(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    v_temp_(iblock, 0, j, i) = 0.25 * v_dzr_(0) * (
        (v_atb_(iblock, n_, 1, j+1, i) - v_atb_(iblock, n_, 2, j+1, i))
      + (v_atb_(iblock, n_, 1, j  , i) - v_atb_(iblock, n_, 2, j  , i)));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_3 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    v_temp_(iblock, KM-1, j, i) = C0;
    return ;
  }
 private:
  const ViewDouble4D v_temp_ = *p_v_temp;
};

class functor_tracer_isoflux_4 {
 public:
  functor_tracer_isoflux_4(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    const int k = std::min(v_kmt_(iblock, j, i),
        v_kmt_(iblock, j+1, i)) - 1;

    if (k != -1) {
      const double fxe = v_dzw_(k) + v_dzw_(k+1);

      const double fxa = 0.5 
          * (v_atb_(iblock, n_, k  , j+1, i) + v_atb_(iblock, n_, k  , j, i));
      const double fxb = 0.5 
          * (v_atb_(iblock, n_, k+1, j+1, i) + v_atb_(iblock, n_, k+1, j, i));

      const double fxc = v_dzwr_(k) * (fxb * fxe - fxa * v_dzw_(k+1));

      v_temp_(iblock, k, j, i) = v_dzr_(k) * (0.5 * (fxa + fxb) - fxc);
    }
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble1D v_dzw_  = *p_v_dzw;
  const ViewDouble1D v_dzwr_ = *p_v_dzwr;
  const ViewInt3D    v_kmt_  = *p_v_kmt;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_5 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &i) const {
    const int iblock = 0;
    v_work_1_(iblock, k, 0, i) = C0;
    return ;
  }
 private:
  const ViewDouble4D v_work_1_ = *p_v_work_1;
};

class functor_tracer_isoflux_6 {
 public:
  functor_tracer_isoflux_6(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_work_1_(iblock, k, j, i) = v_hts_(iblock, j, i) * (v_ahisop_(iblock, j, i)
        * (v_atb_(iblock, n_, k+1, j+1, i) - v_atb_(iblock, n_, k+1, j, i))
            / v_hue_(iblock, j, i) + v_ahisop_(iblock, j, i)
                * v_k2_(iblock, 0, j, k+1, i) * v_temp_(iblock, k, j, i))
                    * v_vit_(iblock, k, j, i) * v_vit_(iblock, k, j+1, i);
                     
    v_ddy_iso_(iblock, n_, k, j, i) = v_work_1_(iblock, k, j, i);
    return ;
  }
 private:
  const int n_;
  const ViewDouble3D v_ahisop_  = *p_v_ahisop;
  const ViewDouble3D v_hue_     = *p_v_hue;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble4D v_vit_     = *p_v_vit; 
  const ViewDouble4D v_temp_    = *p_v_temp;
  const ViewDouble4D v_work_1_  = *p_v_work_1;
  const ViewDouble5D v_ddy_iso_ = *p_v_ddy_iso;
  const ViewDouble5D v_atb_     = *p_v_atb;
  const ViewDouble5D v_k2_      = *p_v_k2;
};

class functor_tracer_isoflux_7 {
 public:
  functor_tracer_isoflux_7(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_temp_(iblock, k, j, i) = P25 * v_dzr_(k)
        * (v_atb_(iblock, n_, k, j, i+1) - v_atb_(iblock, n_, k+2, j, i+1)
         + v_atb_(iblock, n_, k , j, i ) - v_atb_(iblock, n_, k+2 , j, i ));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_8 {
 public:
  functor_tracer_isoflux_8(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    v_temp_(iblock, 0, j, i) = P25 * v_dzr_(0)
        * (v_atb_(iblock, n_, 1, j, i+1) - v_atb_(iblock, n_, 1, 2, i+1)
         + v_atb_(iblock, n_, 1, j, i  ) - v_atb_(iblock, n_, 1, 2, i  ));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_9 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    v_temp_(iblock, KM-1, j, i) = C0;
    return ;
  }
 private:
  const ViewDouble4D v_temp_ = *p_v_temp;
};

class functor_tracer_isoflux_10 {
 public:
  functor_tracer_isoflux_10(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const int k = std::min(v_kmt_(iblock, j, i),
        v_kmt_(iblock, j, i+1)) - 1;

    if (k != -1) {
      const double fxe = v_dzw_(k) + v_dzw_(k+1);

      const double fxa = 0.5 
          * (v_atb_(iblock, n_, k  , j, i) + v_atb_(iblock, n_, k  , j, i+1));
      const double fxb = 0.5 
          * (v_atb_(iblock, n_, k+1, j, i) + v_atb_(iblock, n_, k+1, j, i+1));

      const double fxc = v_dzwr_(k) * (fxb * fxe - fxa * v_dzw_(k+1));

      v_temp_(iblock, k, j, i) = v_dzr_(k) * (P5 * (fxa + fxb) - fxc);
    }
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_  = *p_v_dzr;
  const ViewDouble1D v_dzw_  = *p_v_dzw;
  const ViewDouble1D v_dzwr_ = *p_v_dzwr;
  const ViewInt3D    v_kmt_  = *p_v_kmt;
  const ViewDouble4D v_temp_ = *p_v_temp;
  const ViewDouble5D v_atb_  = *p_v_atb;
};

class functor_tracer_isoflux_11 {
 public:
  functor_tracer_isoflux_11(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_work_2_(iblock, k, j, i) = v_htw_(iblock, j, i+1) * (v_ahisop_(iblock, j, i)
        * (v_atb_(iblock, n_, k+1, j, i+1) - v_atb_(iblock, n_, k+1, j, i))
            / v_hun_(iblock, j, i+1) + v_ahisop_(iblock, j, i)
                * v_k1_(iblock, 0, j, k+1, i) * v_temp_(iblock, k, j, i))
                    * v_vit_(iblock, k, j, i+1) * v_vit_(iblock, k, j, i);
    return ;
  }
 private:
  const int n_;
  const ViewDouble3D v_ahisop_ = *p_v_ahisop;
  const ViewDouble3D v_hun_    = *p_v_hun;
  const ViewDouble3D v_htw_    = *p_v_htw;
  const ViewDouble4D v_vit_    = *p_v_vit; 
  const ViewDouble4D v_temp_   = *p_v_temp;
  const ViewDouble4D v_work_2_ = *p_v_work_2;
  const ViewDouble5D v_atb_    = *p_v_atb;
  const ViewDouble5D v_k1_     = *p_v_k1;
};

class functor_tracer_isoflux_12 {
 public:
  functor_tracer_isoflux_12(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_work_3_(iblock, k, j, i) = v_ahisop_(iblock, j, i) * P25 * v_vit_(iblock, k, j, i)
        * (v_k3_(iblock, 0, j, k, i)
           * (v_vit_(iblock, k  , j, k  , i-1) 
                * (v_atb_(iblock, n_, k+1, j, i) - v_atb_(iblock, n_, k+1, j, i-1))
                    / v_hun_(iblock, j, i) 
            + v_vit_(iblock, k-1, j, i-1) 
                * (v_atb_(iblock, n_, k, j, i) - v_atb_(iblock, n_, k, j, i-1))
                    / v_hun_(iblock, j, i) 
            + v_vit_(iblock, k  , j, i+1) 
                * (v_atb_(iblock, n_, k+1, j, i+1) - v_atb_(iblock, n_, k+1, j, i))
                    / v_hun_(iblock, j, i+1) 
            + v_vit_(iblock, k-1, j, i+1) 
                * (v_atb_(iblock, n_, k, j, i+1) - v_atb_(iblock, n_, k, j, i))
                    / v_hun_(iblock, j, i+1))
         + v_k3_(iblock, 1, j, k, i)
           * (v_vit_(iblock, k  , j-1, i) 
                * (v_atb_(iblock, n_, k+1, j, i) - v_atb_(iblock, n_, k+1, j-1, i))
                    / v_hue_(iblock, j-1, i) 
            + v_vit_(iblock, k-1, j-1, i) 
                * (v_atb_(iblock, n_, k, j, i) - v_atb_(iblock, n_, k, j-1, i))
                    / v_hue_(iblock, j-1, i) 
            + v_vit_(iblock, k  , j+1, i) 
                * (v_atb_(iblock, n_, k+1, j+1, i) - v_atb_(iblock, n_, k+1, j, i))
                    / v_hue_(iblock, j, i) 
            + v_vit_(iblock, k-1, j+1, i) 
                * (v_atb_(iblock, n_, k, j+1, i) - v_atb_(iblock, n_, k, j, i))
                    / v_hue_(iblock, j, i)));
    return ;
  }
 private:
  const int n_;
  const ViewDouble3D v_ahisop_ = *p_v_ahisop;
  const ViewDouble3D v_hue_    = *p_v_hue;
  const ViewDouble3D v_hun_    = *p_v_hun;
  const ViewDouble4D v_vit_    = *p_v_vit; 
  const ViewDouble4D v_work_3_ = *p_v_work_3;
  const ViewDouble5D v_atb_    = *p_v_atb;
  const ViewDouble5D v_k3_     = *p_v_k3;
};

class functor_tracer_isoflux_13 {
 public:

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {

    const int iblock = 0;
    v_work_3_(iblock, 0 , j, i) = C0; 
    v_work_3_(iblock, km, j, i) = C0; 
    return ;
  }
 private:
  const int km = CppParamMod::KM;
  const ViewDouble4D v_work_3_ = *p_v_work_3;
};

class functor_tracer_isoflux_14 {
 public:
  functor_tracer_isoflux_14(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_tf_(iblock, k, j, i) +=
        v_tarea_r_(iblock, j, i)
            * (v_work_1_(iblock, k, j, i) - v_work_1_(iblock, k  , j-1, i  ))
      + v_tarea_r_(iblock, j, i)
            * (v_work_2_(iblock, k, j, i) - v_work_2_(iblock, k  , j  , i-1))
      + v_dzr_(k)
            * (v_work_3_(iblock, k, j, i) - v_work_3_(iblock, k+1, j  , i  ));

    v_dx_iso_(iblock, n_, k, j, i) = v_tarea_r_(iblock, j, i)
        * (v_work_2_(iblock, k, j, i) - v_work_2_(iblock, k  , j  , i-1));
    v_dy_iso_(iblock, n_, k, j, i) = v_tarea_r_(iblock, j, i)
        * (v_work_1_(iblock, k, j, i) - v_work_1_(iblock, k  , j-1, i  ));
    v_dz_iso_(iblock, n_, k, j, i) = v_dzr_(k)
        * (v_work_3_(iblock, k, j, i) - v_work_3_(iblock, k+1, j  , i  ));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_     = *p_v_dzr;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_tf_      = *p_v_tf;
  const ViewDouble4D v_work_1_  = *p_v_work_1;
  const ViewDouble4D v_work_2_  = *p_v_work_2;
  const ViewDouble4D v_work_3_  = *p_v_work_3;
  const ViewDouble5D v_dx_iso_  = *p_v_dx_iso;
  const ViewDouble5D v_dy_iso_  = *p_v_dy_iso;
  const ViewDouble5D v_dz_iso_  = *p_v_dz_iso;
};

class functor_tracer_isoflux_15 {
 public:
  functor_tracer_isoflux_15(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_work_1_(iblock, k, j, i) = v_adv_vntiso_(iblock, j, k, i)
        * (v_atb_(iblock, n_, k+1, j+1, i) + v_atb_(iblock, n_, k+1, j, i))
            * v_hts_(iblock, j, i);
    return ;
  }
 private:
  const int n_;
  const ViewDouble3D v_hts_        = *p_v_hts;
  const ViewDouble4D v_work_1_     = *p_v_work_1;
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble5D v_atb_        = *p_v_atb;
};

class functor_tracer_isoflux_16 {
 public:
  functor_tracer_isoflux_16(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_work_2_(iblock, k, j, i) = v_adv_vetiso_(iblock, j, k, i)
        * (v_atb_(iblock, n_, k+1, j, i+1) + v_atb_(iblock, n_, k+1, j, i))
            * v_htw_(iblock, j, i+1);
    return ;
  }
 private:
  const int n_;
  const ViewDouble3D v_htw_        = *p_v_htw;
  const ViewDouble4D v_work_2_     = *p_v_work_2;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
  const ViewDouble5D v_atb_        = *p_v_atb;
};

class functor_tracer_isoflux_17 {
 public:

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    v_work_3_(iblock, 0 , j, i) = C0; 
    v_work_3_(iblock, k_, j, i) = C0; 
    return ;
  }
 private:
  const int k_ = KM;
  const ViewDouble4D v_work_3_ = *p_v_work_3;
};

class functor_tracer_isoflux_18 {
 public:
  functor_tracer_isoflux_18(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {

    const int iblock = 0;

    v_work_3_(iblock, k, j, i) = v_adv_vbtiso_(iblock, j, k, i)
        * (v_atb_(iblock, n_, k+1, j, i) + v_atb_(iblock, n_, k, j, i));
    return ;
  }
 private:
  const int n_;
  const ViewDouble4D v_work_3_     = *p_v_work_3;
  const ViewDouble4D v_adv_vbtiso_ = *p_v_adv_vbtiso;
  const ViewDouble5D v_atb_        = *p_v_atb;
};

class functor_tracer_isoflux_19 {
 public:
  functor_tracer_isoflux_19(const int &n) : n_(n) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {

    const int iblock = 0;

    v_tf_(iblock, k, j, i) = v_tf_(iblock, k, j, i)
      - P5 * (v_work_1_(iblock, k, j, i) - v_work_1_(iblock, k, j-1, i  )
            + v_work_2_(iblock, k, j, i) - v_work_2_(iblock, k, j  , i-1))
                * v_tarea_r_(iblock, j, i)
      - P5 * v_dzr_(k)
           * (v_work_3_(iblock, k, j, i) - v_work_3_(iblock, k+1, j, i));

    v_ax_iso_(iblock, n_, k, j, i) = - P5
        * (v_work_2_(iblock, k, j, i) - v_work_2_(iblock, k, j  , i-1))
            * v_tarea_r_(iblock, j, i);

    v_ay_iso_(iblock, n_, k, j, i) = - P5
        * (v_work_1_(iblock, k, j, i) - v_work_2_(iblock, k, j-1, i  ))
            * v_tarea_r_(iblock, j, i);
    
    v_az_iso_(iblock, n_, k, j, i) = - P5 * v_dzr_(k)
        * (v_work_3_(iblock, k, j, i) - v_work_3_(iblock, k+1, j, i));
    return ;
  }
 private:
  const int n_;
  const ViewDouble1D v_dzr_     = *p_v_dzr;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_tf_      = *p_v_tf;
  const ViewDouble4D v_work_1_  = *p_v_work_1;
  const ViewDouble4D v_work_2_  = *p_v_work_2;
  const ViewDouble4D v_work_3_  = *p_v_work_3;
  const ViewDouble5D v_ax_iso_  = *p_v_ax_iso;
  const ViewDouble5D v_ay_iso_  = *p_v_ay_iso;
  const ViewDouble5D v_az_iso_  = *p_v_az_iso;
  const ViewDouble5D v_atb_     = *p_v_atb;
};
// End isoflux
#else  // ISO
#endif // iSO

#ifdef NODIAG

#ifdef ISO
// isopyc
class functor_tracer_isopyc_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    for (int m = 0; m < NRPL; ++m) {
      v_rhoi_(i, 0, j, m, iblock) = 0.0;
    }
    return ;
  }
 private:
  const ViewDouble5D v_rhoi_ = *p_v_rhoi;
};

class functor_tracer_isopyc_2 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    for (int m = 0; m < NRPL; ++m) {
      const double tref0 = v_to_(v_krplin_(m) - 1);
      const double sref0 = v_so_(v_krplin_(m) - 1);
      const double tq = v_atb_(i, j, k+1, 0, iblock) - tref0;
      const double sq = v_atb_(i, j, k+1, 1, iblock) - sref0;
      const int    kq = v_krplin_(m) - 1;
      v_rhoi_(i, k+1, j, m, iblock) = dens(tq, sq, kq);
    }
    return ;
  }

 private:
  const ViewInt1D    v_krplin_ = *p_v_krplin;
  const ViewDouble1D v_to_     = *p_v_to;
  const ViewDouble1D v_so_     = *p_v_so;
  const ViewDouble2D v_c_      = *p_v_c;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
  const ViewDouble5D v_atb_    = *p_v_atb;

  KOKKOS_INLINE_FUNCTION double
      dens(const double &tq, const double &sq, const int &kk) const {
    double dens;
    dens = (v_c_(kk, 0) + (v_c_(kk, 3) + v_c_(kk, 6) * sq) * sq +
           (v_c_(kk, 2) +  v_c_(kk, 7) * sq + v_c_(kk, 5) * tq) * tq) * tq +
           (v_c_(kk, 1) + (v_c_(kk, 4) + v_c_(kk, 8) * sq) * sq) * sq;
    return dens;
  }
};

// k2_3
class functor_tracer_isopyc_k2_3_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    const double c1e10 = 1.0e10;

    const int m = static_cast<int>(v_kisrpl_(k)) - 1;

    const double fxd = c1e10 * P25 * v_dzr_(k);

    v_e_(iblock, 2, j, k, i) = fxd
        * (v_rhoi_(i, k, j  , m, iblock) - v_rhoi_(i, k+2, j  , m, iblock)
         + v_rhoi_(i, k, j+1, m, iblock) - v_rhoi_(i, k+2, j+1, m, iblock));
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k2_3_2 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const double fxd = c1e10 * v_dzr_(0);
    const double fxe = v_dzw_(0) + v_dzw_(1);

    const int m = static_cast<int>(v_kisrpl_(0)) - 1;

    const double fxa = p5 * 
        (v_rhoi_(i, 2, j, m, iblock) + v_rhoi_(i, 2, j+1, m, iblock));
    const double fxb = p5 * 
        (v_rhoi_(i, 1, j, m, iblock) + v_rhoi_(i, 1, j+1, m, iblock));

    const double fxc = v_dzwr_(1) * (fxb * fxe - fxa * v_dzw_(0));

    v_e_(i, 0, j, 2, iblock) = fxd * (fxc - p5 * (fxa + fxb));
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble1D v_dzw_    = *p_v_dzw;
  const ViewDouble1D v_dzwr_   = *p_v_dzwr;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k2_3_3 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const double c0  = 0.0;

    v_e_(i, KM-1, j, 2, iblock) = c0;
    return ;
  }

 private:
  const ViewDouble5D v_e_ = *p_v_e;
};

class functor_tracer_isopyc_k2_3_4 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {

    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const int k = std::min(v_kmt_(iblock, j, i), v_kmt_(iblock, j+1, i));

    if (k != 0) {
      const double fxe = v_dzw_(k-1) + v_dzw_(k);

      const int m = static_cast<int>(v_kisrpl_(k-1)) - 1;

      const double fxa = p5 * 
          (v_rhoi_(i, k-1, j, m, iblock) + v_rhoi_(i, k-1, j+1, m, iblock));
      const double fxb = p5 * 
          (v_rhoi_(i, k  , j, m, iblock) + v_rhoi_(i, k  , j+1, m, iblock));

      const double fxc = v_dzwr_(k-1) * (fxb * fxe - fxa * v_dzw_(k)); 

      v_e_(i, k-1, j, 2, iblock) = v_dzr_(k-1) * c1e10 
          * (p5 * (fxa + fxb) - fxc);
    }
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble1D v_dzw_    = *p_v_dzw;
  const ViewDouble1D v_dzwr_   = *p_v_dzwr;
  const ViewInt3D    v_kmt_    = *p_v_kmt;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k2_3_5 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {

    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const int m = static_cast<int>(v_kisrpl_(k)) - 1;

    v_e_(iblock, 0, j, k, i) = p5 * c1e10
        / (v_dxu_(iblock, j, i) + v_dxu_(iblock, j, i+1))
        * (v_rhoi_(i+1, k+1, j+1, m, iblock) - v_rhoi_(i-1, k+1, j+1, m, iblock)
         + v_rhoi_(i+1, k+1, j  , m, iblock) - v_rhoi_(i-1, k+1, j  , m, iblock));

    v_e_(iblock, 1, j, k, i) = v_tmask_(iblock, j, k, i)
        * v_tmask_(i, k, j+1, iblock) / v_hue_(iblock, j, i) * c1e10
            * (v_rhoi_(i, k+1, j+1, m, iblock) - v_rhoi_(i, k+1, j, m, iblock));
    return ;
  }

 private:
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble3D v_dxu_    = *p_v_dxu;
  const ViewDouble3D v_hue_    = *p_v_hue;
  const ViewDouble4D v_tmask_  = *p_v_tmask;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k2_3_6 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {

    const int iblock = 0;

    const double c0 = 0.0;
    const double c1 = 1.0;

    const double olmask =
        v_tmask_(i-1, k, j+1, iblock) * v_tmask_(i+1, k, j+1, iblock)
      * v_tmask_(i-1, k, j  , iblock) * v_tmask_(i+1, k, j  , iblock);

    if (olmask < c1) {
      v_e_(iblock, 0, j, k, i) = c0;
    }
    return ;
  }

 private:
  const ViewDouble4D v_tmask_ = *p_v_tmask;
  const ViewDouble5D v_e_     = *p_v_e;
};
#ifdef LDD97
class functor_tracer_isopyc_k2_3_7 {
 public:
  functor_tracer_isopyc_k2_3_7(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_f1_(iblock, k, j, i) = 0.0;
    v_f2_(iblock, k, j, i) = 0.0;
    return ;
  }

 private:
  const ViewDouble4D v_f1_;
  const ViewDouble4D v_f2_;
};

class functor_tracer_isopyc_k2_3_8 {
 public:
  functor_tracer_isopyc_k2_3_8(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;
    
    const double chkslp = - sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2)) * slmxr_;

    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }

    const double slopemod = sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2))
            / fabs(v_e_(iblock, 2, j, k, i) + eps);

    v_f1_(iblock, k, j, i) = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));

    const double nondimr = - v_zkt_(k) / (v_rrd2_(iblock, j, i)
        * (slopemod + eps));

    if (nondimr >= 1.0) {
      f2(iblock, k, j, i) = 1.0;
    } else {
      f2(iblock, k, j, i) = 0.5 * (1.0 + sin(PI * (nondimr - 0.5)));
    }

    v_k2_(iblock, 0, j, k+1, i) = (- v_e_(iblock, 1, j, k, i)
        * v_e_(iblock, 2, j, k, i) * v_f1_(iblock, k, j, i)
            * v_f2_(iblock, k, j, i)) / (pow(v_e_(iblock, 2, j, k, i), 2) + eps);
    return ;
  }

 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble4D v_f1_;
  const ViewDouble4D v_f2_;
  const ViewDouble1D v_zkt_  = *p_v_zkt;
  const ViewDouble3D v_rrd2_ = *p_v_rrd2;
  const ViewDouble5D v_e_    = *p_v_e;
  const ViewDouble5D v_k2_   = *p_v_k2;
};
#else  // LDD97
class functor_tracer_isopyc_k2_3_9 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;

    const double chkslp = - sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2)) * slmxr_;
    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }
    v_k2_(iblock, 0, j, k+1, i) = (- v_e_(iblock, 1, j, k, i)
        * v_e_(iblock, 2, j, k, i) * v_fzisop_(k))
            / (pow(v_e_(iblock, 2, j, k, i), 2) + eps);
    return ;
  }

 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble1D v_fzisop_ = *p_v_fzisop;
  const ViewDouble5D v_e_    = *p_v_e;
  const ViewDouble5D v_k2_   = *p_v_k2;
};

#endif // LDD97
// End k2_3

// k1_3
class functor_tracer_isopyc_k1_3_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double p25   = 0.25;
    const double c1e10 = 1.0e10;

    const int m = static_cast<int>(v_kisrpl_(k)) - 1;

    const double fxd = c1e10 * p25 * v_dzr_(k);

    v_e_(iblock, 2, j, k, i) = fxd
        * (v_rhoi_(i  , k, j, m, iblock) - v_rhoi_(i  , k+2, j, m, iblock)
         + v_rhoi_(i+1, k, j, m, iblock) - v_rhoi_(i+1, k+2, j, m, iblock));
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k1_3_2 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const double fxd = c1e10 * v_dzr_(0);

    const double fxe = v_dzw_(0) + v_dzw_(1);

    const int m = static_cast<int>(v_kisrpl_(0)) - 1;

    const double fxa = p5 * 
        (v_rhoi_(i, 2, j, m, iblock) + v_rhoi_(i+1, 2, j, m, iblock));
    const double fxb = p5 * 
        (v_rhoi_(i, 1, j, m, iblock) + v_rhoi_(i+1, 1, j, m, iblock));

    const double fxc = v_dzwr_(1) * (fxb * fxe - fxa * v_dzw_(0));

    v_e_(i, 0, j, 2, iblock) = fxd * (fxc - p5 * (fxa + fxb));
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble1D v_dzw_    = *p_v_dzw;
  const ViewDouble1D v_dzwr_   = *p_v_dzwr;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k1_3_3 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;

    const double c0 = 0.0;

    v_e_(i, KM-1, j, 2, iblock) = c0;
    return ;
  }

 private:
  const ViewDouble5D v_e_ = *p_v_e;
};

class functor_tracer_isopyc_k1_3_4 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {

    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const int k = std::min(v_kmt_(iblock, j, i), v_kmt_(iblock, j, i+1)) - 1;

    if (k != -1) {
      const double fxe = v_dzw_(k) + v_dzw_(k+1);

      const int m = static_cast<int>(v_kisrpl_(k)) - 1;

      const double fxa = p5 * 
          (v_rhoi_(i, k  , j, m, iblock) + v_rhoi_(i+1, k  , j, m, iblock));
      const double fxb = p5 * 
          (v_rhoi_(i, k+1, j, m, iblock) + v_rhoi_(i+1, k+1, j, m, iblock));

      const double fxc = v_dzwr_(k) * (fxb * fxe - fxa * v_dzw_(k+1)); 

      v_e_(iblock, 2, j, k, i) = v_dzr_(k) * c1e10 
          * (p5 * (fxa + fxb) - fxc);
    }
    return ;
  }

 private:
  const ViewDouble1D v_dzr_    = *p_v_dzr;
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble1D v_dzw_    = *p_v_dzw;
  const ViewDouble1D v_dzwr_   = *p_v_dzwr;
  const ViewInt3D    v_kmt_    = *p_v_kmt;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k1_3_5 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {

    const int iblock = 0;

    const double p5    = 0.5;
    const double c1e10 = 1.0e10;

    const int m = static_cast<int>(v_kisrpl_(k)) - 1;

    v_e_(iblock, 0, j, k, i) = v_tmask_(iblock, j, k, i)
        * v_tmask_(i+1, k, j, iblock) / v_hun_(iblock, j, i+1) * c1e10
            * (v_rhoi_(i+1, k+1, j, m, iblock) - v_rhoi_(i, k+1, j, m, iblock));

    v_e_(iblock, 1, j, k, i) = p5 * c1e10
        / (v_dyu_(i+1, j-1, iblock) + v_dyu_(iblock, j, i+1))
        * (v_rhoi_(i  , k+1, j+1, m, iblock) - v_rhoi_(i  , k+1, j-1, m, iblock)
         + v_rhoi_(i+1, k+1, j+1, m, iblock) - v_rhoi_(i+1, k+1, j-1, m, iblock));
    return ;
  }

 private:
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble3D v_dyu_    = *p_v_dyu;
  const ViewDouble3D v_hun_    = *p_v_hun;
  const ViewDouble4D v_tmask_  = *p_v_tmask;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k1_3_6 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {

    const int iblock = 0;

    const double c0 = 0.0;
    const double c1 = 1.0;

    const double olmask =
        v_tmask_(i  , k, j-1, iblock) * v_tmask_(i  , k, j+1, iblock)
      * v_tmask_(i+1, k, j-1, iblock) * v_tmask_(i+1, k, j+1, iblock);

    if (olmask < c1) {
      v_e_(iblock, 1, j, k, i) = c0;
    }
    return ;
  }

 private:
  const ViewDouble4D v_tmask_ = *p_v_tmask;
  const ViewDouble5D v_e_     = *p_v_e;
};

#ifdef LDD97
class functor_tracer_isopyc_k1_3_7 {
 public:
  functor_tracer_isopyc_k1_3_7(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_f1_(iblock, k, j, i) = 1.0;
    v_f2_(iblock, k, j, i) = 1.0;
    return ;
  }

 private:
  const ViewDouble4D v_f1_;
  const ViewDouble4D v_f2_;
};

class functor_tracer_isopyc_k1_3_8 {
 public:
  functor_tracer_isopyc_k1_3_8(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;
    
    const double chkslp = - sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2)) * slmxr_;

    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }

    const double slopemod = sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2))
            / fabs(v_e_(iblock, 2, j, k, i) + eps);

    v_f1_(iblock, k, j, i) = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));

    const double nondimr = - v_zkt_(k) / (v_rrd1_(iblock, j, i)
        * (slopemod + eps));

    if (nondimr >= 1.0) {
      f2(iblock, k, j, i) = 1.0;
    } else {
      f2(iblock, k, j, i) = 0.5 * (1.0 + sin(PI * (nondimr - 0.5)));
    }

    v_k1_(iblock, 0, j, k+1, i) = (- v_e_(iblock, 1, j, k, i)
        * v_e_(iblock, 2, j, k, i) * v_f1_(iblock, k, j, i)
            * v_f2_(iblock, k, j, i)) / (pow(v_e_(iblock, 2, j, k, i), 2) + eps);
    return ;
  }
 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble4D v_f1_, v_f2_;
  const ViewDouble1D v_zkt_  = *p_v_zkt;
  const ViewDouble3D v_rrd1_ = *p_v_rrd1;
  const ViewDouble5D v_e_    = *p_v_e;
  const ViewDouble5D v_k1_   = *p_v_k1;
};

#else  // LDD97
class functor_tracer_isopyc_k1_3_9 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;

    const double chkslp = - sqrt(pow(v_e_(iblock, 0, j, k, i), 2)
        + pow(v_e_(iblock, 1, j, k, i), 2)) * slmxr_;
    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }
    v_k1_(iblock, 0, j, k+1, i) = (- v_e_(iblock, 0, j, k, i)
        * v_e_(iblock, 2, j, k, i) * v_fzisop_(k))
            / (pow(v_e_(iblock, 2, j, k, i), 2) + eps);
    return ;
  }

 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble1D v_fzisop_ = *p_v_fzisop;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_k1_     = *p_v_k1;
};

#endif // LDD97
// End k1_3

// k3_123
class functor_tracer_isopyc_k3_123_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double p25   = 0.25;
    const double c1e10 = 1.0e10;

    const int m = static_cast<int>(v_kisrpl_(k)) - 1;

    v_e_(i, k-1, j, 0, iblock) = p25 * c1e10 * (
        v_tmask_(i-1, k-1, j, iblock) * v_tmask_(i  , k-1, j, iblock)
            * (v_rhoi_(i  , k  , j, m, iblock) - v_rhoi_(i-1, k  , j, m, iblock))
                / v_hun_(i  , j, iblock)
      + v_tmask_(i  , k-1, j, iblock) * v_tmask_(i+1, k-1, j, iblock)
            * (v_rhoi_(i+1, k  , j, m, iblock) - v_rhoi_(i  , k  , j, m, iblock))
                / v_hun_(iblock, j, i+1)
      + v_tmask_(i-1, k  , j, iblock) * v_tmask_(i  , k  , j, iblock)
            * (v_rhoi_(i  , k+1, j, m, iblock) - v_rhoi_(i-1, k+1, j, m, iblock))
                / v_hun_(i  , j, iblock)
      + v_tmask_(i  , k  , j, iblock) * v_tmask_(i+1, k  , j, iblock)
            * (v_rhoi_(i+1, k+1, j, m, iblock) - v_rhoi_(i  , k+1, j, m, iblock))
                / v_hun_(iblock, j, i+1));

    v_e_(i, k-1, j, 1, iblock) = p25 * c1e10 * (
        v_tmask_(i, k-1, j-1, iblock) * v_tmask_(i, k-1, j  , iblock)
            * (v_rhoi_(i, k  , j  , m, iblock) - v_rhoi_(i, k  , j-1, m, iblock))
                / v_hue_(iblock, j-1, i)
      + v_tmask_(i, k-1, j  , iblock) * v_tmask_(i, k-1, j+1, iblock)
            * (v_rhoi_(i, k  , j+1, m, iblock) - v_rhoi_(i, k  , j  , m, iblock))
                / v_hue_(i, j  , iblock)
      + v_tmask_(i, k  , j-1, iblock) * v_tmask_(i, k  , j  , iblock)
            * (v_rhoi_(i, k+1, j  , m, iblock) - v_rhoi_(i, k+1, j-1, m, iblock))
                / v_hue_(iblock, j-1, i)
      + v_tmask_(i, k  , j  , iblock) * v_tmask_(i, k  , j+1, iblock)
            * (v_rhoi_(i, k+1, j+1, m, iblock) - v_rhoi_(i, k+1, j  , m, iblock))
                / v_hue_(i, j  , iblock));

    v_e_(i, k-1, j, 2, iblock) = v_dzwr_(k)
        * v_tmask_(i, k-1, j, iblock) * v_tmask_(iblock, j, k, i) * c1e10
            * (v_rhoi_(i, k, j, m, iblock) - v_rhoi_(i, k+1, j, m, iblock)); 
    return ;
  }

 private:
  const ViewDouble1D v_kisrpl_ = *p_v_kisrpl;
  const ViewDouble1D v_dzwr_   = *p_v_dzwr;
  const ViewDouble3D v_hue_    = *p_v_hue;
  const ViewDouble3D v_hun_    = *p_v_hun;
  const ViewDouble4D v_tmask_  = *p_v_tmask;
  const ViewDouble5D v_e_      = *p_v_e;
  const ViewDouble5D v_rhoi_   = *p_v_rhoi;
};

class functor_tracer_isopyc_k3_123_2 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) 
          const {
    const int iblock = 0;
    const double c0 = 0.0;
    v_e_(i, KM-1, j, 0, iblock) = c0;
    v_e_(i, KM-1, j, 1, iblock) = c0;
    v_e_(i, KM-1, j, 2, iblock) = c0;
    return ;
  }
 private:
  const ViewDouble5D v_e_ = *p_v_e;
};
#ifdef LDD97
class functor_tracer_isopyc_k3_123_3 {
 public:
  functor_tracer_isopyc_k3_123_3(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;
    v_f1_(iblock, k, j, i) = 0.0;
    v_f2_(iblock, k, j, i) = 0.0;
    return ;
  }

 private:
  const ViewDouble4D v_f1_, v_f2_;
};

class functor_tracer_isopyc_k3_123_4 {
 public:
  functor_tracer_isopyc_k3_123_4(
      const ViewDouble4D &v_f1,
      const ViewDouble4D &v_f2)
    : v_f1_(v_f1), v_f2_(v_f2) {}
 
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;

    const double chkslp = - sqrt(
        pow(v_e_(iblock, 0, j, k, i), 2) + pow(v_e_(iblock, 1, j, k, i), 2))
            * slmxr_;
    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }
    const double slopemod = sqrt(
        pow(v_e_(iblock, 0, j, k, i), 2) + pow(v_e_(iblock, 1, j, k, i), 2))
            / fabs(v_e_(iblock, 2, j, k, i) + eps);

    v_f1_(iblock, k, j, i) = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));

    const double nondimr = - v_zkp_(k) / (v_rrd1_(iblock, j, i)
        * (slopemod + eps));

    if (nondimr >= 1.0) {
      v_f2_(iblock, k, j, i) = 1.0;
    } else {
      v_f2_(iblock, k, j, i) = 0.5 * (1.0 + sin(PI * (nondimr - 0.5)));
    }

    const double ahfctr = 1.0 / (pow(v_e_(iblock, 2, j, k, i), 2) + eps)
        * v_f1_(iblock, k, j, i) * v_f2_(iblock, k, j, i);

    v_k3_(iblock, 0, j, k+1, i) = - v_e_(iblock, 2, j, k, i)
                                  * v_e_(iblock, 0, j, k, i) * ahfctr;
                                   
    v_k3_(i, k+1, j, 1, iblock) = - v_e_(iblock, 2, j, k, i)
                                  * v_e_(iblock, 1, j, k, i) * ahfctr;
                                   
    v_k3_(iblock, 2, j, k+1, i) = (pow(v_e_(iblock, 0, j, k, i), 2)
                                 + pow(v_e_(iblock, 1, j, k, i), 2));
    return ;
  }

 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble4D v_f1_;
  const ViewDouble4D v_f2_;
  const ViewDouble1D v_zkp_  = *p_v_zkp;
  const ViewDouble3D v_rrd1_ = *p_v_rrd1;
  const ViewDouble5D v_e_    = *p_v_e;
  const ViewDouble5D v_k3_   = *p_v_k3;
};
#else //  LDD97
class functor_tracer_isopyc_k3_123_5 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;

    const double eps = 1.0e-25;

    const double chkslp = - sqrt(
        pow(v_e_(iblock, 0, j, k, i), 2) + pow(v_e_(iblock, 1, j, k, i), 2))
            * slmxr_;

    if (v_e_(iblock, 2, j, k, i) > chkslp) {
      v_e_(iblock, 2, j, k, i) = chkslp;
    }

    const double ahfctr = 1.0 / (pow(v_e_(iblock, 2, j, k, i), 2) + eps);

    v_k3_(iblock, 0, j, k+1, i) = - v_e_(iblock, 2, j, k, i)
                                  * v_e_(iblock, 0, j, k, i) * ahfctr;
                                   
    v_k3_(i, k+1, j, 1, iblock) = - v_e_(iblock, 2, j, k, i)
                                  * v_e_(iblock, 1, j, k, i) * ahfctr;
                                   
    v_k3_(iblock, 2, j, k+1, i) = (pow(v_e_(iblock, 0, j, k, i), 2)
                                 + pow(v_e_(iblock, 1, j, k, i), 2)) 
                                     * ahfctr;
    return ;
  }

 private:
  const double slmxr_ = CppIsopycMod::slmxr;

  const ViewDouble5D v_e_  = *p_v_e;
  const ViewDouble5D v_k3_ = *p_v_k3;
};

#endif // LDD97
// End k3_123 

// isoadv
class functor_tracer_isopyc_isoadv_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    const double c0 = 0.0;
    v_adv_vntiso_(iblock, j, k, i) = c0;
    v_adv_vetiso_(iblock, j, k, i) = c0;
    return ;
  }
 private:
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
};
class functor_tracer_isopyc_isoadv_2 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    const double c0 = 0.0;
    v_adv_vbtiso_(iblock, j, k, i) = c0;
    return ;
  }
 private:
  const ViewDouble4D v_adv_vbtiso_ = *p_v_adv_vbtiso;
};

class functor_tracer_isopyc_isoadv_3 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    const double p5 = 0.5;
    const double fxa = - p5 * v_dzr_(k) * v_athkdf_(iblock, j, i);
    v_adv_vntiso_(iblock, j, k, i) = fxa
        * v_tmask_(iblock, j, k, i) * v_tmask_(iblock, j+1, k, i)
            * (v_k2_(iblock, 0, j, k, i) - v_k2_(iblock, 0, j, k+2, i));
    return ;
  }

 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble5D v_k2_         = *p_v_k2;
};

class functor_tracer_isopyc_isoadv_4 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const double p5 = 0.5;
    const double fxa = - p5 * v_dzr_(0) * v_athkdf_(iblock, j, i);
    v_adv_vntiso_(iblock, j, 0, i) = - fxa
        * v_tmask_(iblock, j, 0, i) * v_tmask_(iblock, j+1, 0, i)
            * (v_k2_(iblock, 0, j, 1, i) + v_k2_(iblock, 0, j, 2, i));
    return ;
  }
 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble5D v_k2_         = *p_v_k2;
};

class functor_tracer_isopyc_isoadv_5 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const double p5 = 0.5;

    const int k = std::min(v_kmt_(iblock, j, i), v_kmt_(iblock, j+1, i)) - 1;

    if (k != -1) {
      v_adv_vntiso_(iblock, j, k, i) = - p5 * v_dzr_(k) * v_athkdf_(iblock, j, i) 
          * v_tmask_(iblock, j, k, i) * v_tmask_(iblock, j+1, k, i)
              * (v_k2_(iblock, 0, j, k+1, i) + v_k2_(iblock, 0, j, k, i));
    }
    return ;
  }

 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble5D v_k2_         = *p_v_k2;
};

class functor_tracer_isopyc_isoadv_6 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    const double p5 = 0.5;
    const double fxa = - p5 * v_dzr_(k) * v_athkdf_(iblock, j, i);

    v_adv_vetiso_(iblock, j, k, i) = fxa
        * v_tmask_(iblock, j, k, i) * v_tmask_(iblock, j, k, i+1)
            * (v_k1_(iblock, 0, j, k, i) - v_k1_(iblock, 0, j, k+2, i));
    return ;
  }
 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
  const ViewDouble5D v_k1_         = *p_v_k1;
};

class functor_tracer_isopyc_isoadv_7 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const double p5 = 0.5;
    const double fxa = - p5 * v_dzr_(0) * v_athkdf_(iblock, j, i);

    v_adv_vetiso_(iblock, j, 0, i) = - fxa
        * v_tmask_(iblock, j, 0, i) * v_tmask_(iblock, j, 0, i+1)
            * (v_k1_(iblock, 0, j, 1, i) + v_k1_(iblock, 0, j, 2, i));
    return ;
  }
 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
  const ViewDouble5D v_k1_         = *p_v_k1;
};

class functor_tracer_isopyc_isoadv_8 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) const {
    const int iblock = 0;
    const double p5 = 0.5;
    const int k = std::min(v_kmt_(iblock, j, i), v_kmt_(iblock, j, i+1)) - 1;

    if (k != -1) {
      v_adv_vetiso_(iblock, j, k, i) = - p5 * v_dzr_(k) * v_athkdf_(iblock, j, i)
          * v_tmask_(iblock, j, k, i) * v_tmask_(iblock, j, k, i+1)
              * (v_k1_(iblock, 0, j, k+1, i) + v_k1_(iblock, 0, j, k, i));
    }
    return ;
  }
 private:
  const ViewDouble1D v_dzr_        = *p_v_dzr;
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble3D v_athkdf_     = *p_v_athkdf;
  const ViewDouble4D v_tmask_      = *p_v_tmask;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
  const ViewDouble5D v_k1_         = *p_v_k1;
};

class functor_tracer_isopyc_isoadv_9 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &k, const int &i) 
          const {
    const int iblock = 0;
    v_adv_vbtiso_(i, k+1, j, iblock) = v_dzp_(k)
        * ((v_adv_vetiso_(iblock, j  , k, i  )
          - v_adv_vetiso_(iblock, j  , k, i-1) * v_htw_(iblock, j  , i  ))
              * v_tarea_r_(iblock, j, i)
         + (v_adv_vntiso_(iblock, j  , k, i  ) * v_hts_(iblock, j  , i  )
          - v_adv_vntiso_(iblock, j-1, k, i  ) * v_hts_(iblock, j-1, i  ))
              * v_tarea_r_(iblock, j, i));
    return ;
  }
 private:
  const ViewDouble1D v_dzp_        = *p_v_dzp;
  const ViewDouble3D v_hts_        = *p_v_hts;
  const ViewDouble3D v_htw_        = *p_v_htw;
  const ViewDouble3D v_tarea_r_    = *p_v_tarea_r;
  const ViewDouble4D v_adv_vetiso_ = *p_v_adv_vetiso;
  const ViewDouble4D v_adv_vntiso_ = *p_v_adv_vntiso;
  const ViewDouble4D v_adv_vbtiso_ = *p_v_adv_vbtiso;
};

class functor_tracer_isopyc_isoadv_10 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) 
          const {
    const int iblock = 0;

    for (int k = 1; k <= KM-1; ++k) {
      v_adv_vbtiso_(iblock, j, k, i) +=
          v_adv_vbtiso_(iblock, j, k-1, i);
    }
    return ;
  }

 private:
  const ViewDouble4D v_adv_vbtiso_ = *p_v_adv_vbtiso;
};

class functor_tracer_isopyc_isoadv_11 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j, const int &i) 
          const {
    const int iblock = 0;
    const double c0 = 0.0;
    v_adv_vbtiso_(i, v_kmt_(iblock, j, i), j, iblock) = c0;
    return ;
  }

 private:
  const ViewInt3D    v_kmt_        = *p_v_kmt;
  const ViewDouble4D v_adv_vbtiso_ = *p_v_adv_vbtiso;
};
// End isoadv
// End isopyc
#endif // ISO

// advection_tracer(wkd[iblock], wkb[iblock], ws[iblock], at[iblock][n], adv_tt, iblock, n);


// class functor_tracer_advection_tracer_3 {
//  public:
//   functor_tracer_advection_tracer_3(const int &iblock)
//     : iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &k, const int &j, const int &i) 
//           const {
//     v_v_sface_(k, j, i) = (
//         v_wkb_(iblock_, k, j, i  ) * v_dxu_(iblock_, j, i  )
//       + v_wkb_(iblock_, k, j, i+1) * v_dxu_(iblock_, j, i+1)) * P25;
//     return ;
//   }

//  private:
//   const int iblock_;
//   const ViewDouble3D v_dxu_     = *p_v_dxu;
//   const ViewDouble3D v_v_sface_ = *p_v_v_sface;
//   const ViewDouble4D v_wkb_     = *p_v_wkb;
// };

// class functor_tracer_advection_tracer_4 {
//  public:
//   functor_tracer_advection_tracer_4(const int &iblock)
//     : iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &k, const int &j, const int &i) 
//           const {
//     v_u_wface_(k, j, i) = (
//         v_wkd_(iblock_, k, j-1, i  ) * v_dyu_(iblock_, j-1, i  )
//       + v_wkd_(iblock_, k, j  , i  ) * v_dyu_(iblock_, j  , i  )) * P25;
//     return ;
//   }

//  private:
//   const int iblock_;
//   const ViewDouble3D v_dyu_     = *p_v_dyu;
//   const ViewDouble3D v_u_wface_ = *p_v_u_wface;
//   const ViewDouble4D v_wkd_     = *p_v_wkd;
// };

// class functor_tracer_advection_tracer_5 {
//  public:
//   functor_tracer_advection_tracer_5(const int &n, const int &iblock)
//     : n_(n), iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &j, const int &i) const {
//     v_adv_tt_(0, j, i) = (
//         - v_u_wface_(0, j  , i  ) 
//             * (v_at_(iblock_, n_, 0, j  , i  ) - v_at_(iblock_, n_, 0, j  , i-1))
//         - v_u_wface_(0, j  , i+1) 
//             * (v_at_(iblock_, n_, 0, j  , i+1) - v_at_(iblock_, n_, 0, j  , i  ))
//         - v_v_sface_(0, j  , i  ) 
//             * (v_at_(iblock_, n_, 0, j+1, i  ) - v_at_(iblock_, n_, 0, j  , i  ))
//         - v_v_sface_(0, j-1, i  ) 
//             * (v_at_(iblock_, n_, 0, j  , i  ) - v_at_(iblock_, n_, 0, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i);

//     v_ax_(iblock_, n_, 0, j, i) += (
//         - v_u_wface_(0, j  , i  ) 
//             * (v_at_(iblock_, n_, 0, j  , i  ) - v_at_(iblock_, n_, 0, j  , i-1))
//         - v_u_wface_(0, j  , i+1) 
//             * (v_at_(iblock_, n_, 0, j  , i+1) - v_at_(iblock_, n_, 0, j  , i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);
                 
//     v_ay_(iblock_, n_, 0, j, i) += (
//         - v_v_sface_(0, j  , i  ) 
//             * (v_at_(iblock_, n_, 0, j+1, i  ) - v_at_(iblock_, n_, 0, j  , i  ))
//         - v_v_sface_(0, j-1, i  ) 
//             * (v_at_(iblock_, n_, 0, j  , i  ) - v_at_(iblock_, n_, 0, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);

//     const double adv_z2 = v_ws_(iblock_, 1, j, i)
//         * (v_at_(iblock_, n_, 0, j, i) - v_at_(iblock_, n_, 1, j, i));

//     v_adv_tt_(0, j, i) -= P5 * v_odzp_(0) * adv_z2;

//     v_az_(iblock_, n_, 0, j, i) -= P5 * v_odzp_(0) * adv_z2
//         / static_cast<double>(nss_);
//     return ;
//   }

//  private:
//   const int n_, iblock_;

//   const double nss_ = CppPconstMod::nss;

//   const ViewDouble1D v_odzp_    = *p_v_odzp;
//   const ViewDouble3D v_adv_tt_  = *p_v_adv_tt;
//   const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
//   const ViewDouble3D v_u_wface_ = *p_v_u_wface;
//   const ViewDouble3D v_v_sface_ = *p_v_v_sface;
//   const ViewDouble4D v_ws_      = *p_v_ws;
//   const ViewDouble5D v_ax_      = *p_v_ax;
//   const ViewDouble5D v_ay_      = *p_v_ay;
//   const ViewDouble5D v_az_      = *p_v_az;
//   const ViewDouble5D v_at_      = *p_v_at;
// };

// class functor_tracer_advection_tracer_6 {
//  public:
//   functor_tracer_advection_tracer_6(const int &n, const int &iblock)
//     : n_(n), iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &k, const int &j, const int &i) 
//           const {
//     v_adv_tt_(k, j, i) = (
//         - v_u_wface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j  , i-1))
//         - v_u_wface_(k, j  , i+1) 
//             * (v_at_(iblock_, n_, k, j  , i+1) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j+1, i  ) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j-1, i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i);

//     v_ax_(iblock_, n_, k, j, i) += (
//         - v_u_wface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j  , i-1))
//         - v_u_wface_(k, j  , i+1) 
//             * (v_at_(iblock_, n_, k, j  , i+1) - v_at_(iblock_, n_, k, j  , i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);
                 
//     v_ay_(iblock_, n_, k, j, i) += (
//         - v_v_sface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j+1, i  ) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j-1, i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);

//     const double adv_z1 = v_ws_(iblock_, k  , j, i)
//         * (v_at_(iblock_, n_, k-1, j, i) - v_at_(iblock_, n_, k  , j, i));
         
//     const double adv_z2 = v_ws_(i, j, k+1, iblock_)
//         * (v_at_(iblock_, n_, k  , j, i) - v_at_(iblock_, n_, k+1, j, i));

//     v_adv_tt_(k, j, i) -= P5 * v_odzp_(k) * (adv_z1 + adv_z2);

//     v_az_(iblock_, n_, k, j, i) -= P5 * v_odzp_(k) * (adv_z1 + adv_z2)
//         / static_cast<double>(nss_);
//     return ;
//   }

//  private:
//   const int n_, iblock_;

//   const double nss_ = CppPconstMod::nss;
//   const ViewDouble1D v_odzp_    = *p_v_odzp;
//   const ViewDouble3D v_adv_tt_  = *p_v_adv_tt;
//   const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
//   const ViewDouble3D v_u_wface_ = *p_v_u_wface;
//   const ViewDouble3D v_v_sface_ = *p_v_v_sface;
//   const ViewDouble4D v_ws_      = *p_v_ws;
//   const ViewDouble5D v_ax_      = *p_v_ax;
//   const ViewDouble5D v_ay_      = *p_v_ay;
//   const ViewDouble5D v_az_      = *p_v_az;
//   const ViewDouble5D v_at_      = *p_v_at;
// };

// class functor_tracer_advection_tracer_7 {
//  public:
//   functor_tracer_advection_tracer_7(const int &n, const int &iblock)
//     : n_(n), iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &j, const int &i) const {
//     const int k = KM - 1;
//     v_adv_tt_(k, j, i) = (
//         - v_u_wface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j  , i-1))
//         - v_u_wface_(k, j  , i+1) 
//             * (v_at_(iblock_, n_, k, j  , i+1) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j+1, i  ) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j-1, i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i);

//     v_ax_(iblock_, n_, k, j, i) += (
//         - v_u_wface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j  , i-1))
//         - v_u_wface_(k, j  , i+1) 
//             * (v_at_(iblock_, n_, k, j  , i+1) - v_at_(iblock_, n_, k, j  , i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);
                 
//     v_ay_(iblock_, n_, k, j, i) += (
//         - v_v_sface_(k, j  , i  ) 
//             * (v_at_(iblock_, n_, k, j+1, i  ) - v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j-1, i  ) 
//             * (v_at_(iblock_, n_, k, j  , i  ) - v_at_(iblock_, n_, k, j-1, i  )))
//                 * v_tarea_r_(iblock_, j, i) / static_cast<double>(nss_);

//     const double adv_z1 = v_ws_(iblock_, k, j, i)
//         * (v_at_(iblock_, n_, k-1, j, i) - v_at_(iblock_, n_, k, j, i));
         
//     v_adv_tt_(k, j, i) -= P5 * v_odzp_(k) * adv_z1;

//     v_az_(iblock_, n_, k, j, i) -= P5 * v_odzp_(k) * adv_z1
//         / static_cast<double>(nss_);
//     return ;
//   }

//  private:
//   const int n_, iblock_;
//   const double nss_ = CppPconstMod::nss;

//   const ViewDouble1D v_odzp_    = *p_v_odzp;
//   const ViewDouble3D v_adv_tt_  = *p_v_adv_tt;
//   const ViewDouble3D v_u_wface_ = *p_v_u_wface;
//   const ViewDouble3D v_v_sface_ = *p_v_v_sface;
//   const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
//   const ViewDouble4D v_ws_      = *p_v_ws;
//   const ViewDouble5D v_ax_      = *p_v_ax;
//   const ViewDouble5D v_ay_      = *p_v_ay;
//   const ViewDouble5D v_az_      = *p_v_az;
//   const ViewDouble5D v_at_      = *p_v_at;
// };

// class functor_tracer_advection_tracer_8 {
//  public:
//   functor_tracer_advection_tracer_8(const int &n, const int &iblock)
//     : n_(n), iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void
//       operator () (const int &k, const int &j, const int &i) 
//           const {
//     v_adv_tt_(k, j, i) = (
//         - v_u_wface_(k, j  , i  )
//             * (v_at_(iblock_, n_, k, j  , i  ) + v_at_(iblock_, n_, k, j  , i-1))
//         + v_u_wface_(k, j  , i+1)                                    
//             * (v_at_(iblock_, n_, k, j  , i+1) + v_at_(iblock_, n_, k, j  , i  ))
//         - v_v_sface_(k, j-1, i  )                                    
//             * (v_at_(iblock_, n_, k, j  , i  ) + v_at_(iblock_, n_, k, j-1, i  ))
//         + v_v_sface_(k, j  , i  )                                    
//             * (v_at_(iblock_, n_, k, j+1, i  ) + v_at_(iblock_, n_, k, j  , i  )))
//                 * v_tarea_r_(iblock_, j, i);
//     double adv_z1, adv_z2;
//     if (k == 0) {
//       adv_z1 = 0.0; 
//     } else {
//       adv_z1 = v_ws_(iblock_, k, j, i) 
//           * (v_at_(iblock_, n_, k-1, j, i) + v_at_(iblock_, n_, k, j, i)) * P5;
//     }
//     if (k == KM-1) {
//       adv_z2 = 0.0;
//     } else {
//       adv_z2 = v_ws_(iblock_, k+1, j, i) 
//           * (v_at_(iblock_, n_, k, j, i) + v_at_(iblock_, n_, k+1, j, i)) * P5;
//     }
//     v_adv_tt_(k, j, i) -= v_odzp_(k) * (adv_z2 - adv_z1);
//     return ;
//   }

//  private:
//   const int n_, iblock_;
//   const ViewDouble1D v_odzp_    = *p_v_odzp;
//   const ViewDouble3D v_adv_tt_  = *p_v_adv_tt;
//   const ViewDouble3D v_u_wface_ = *p_v_u_wface;
//   const ViewDouble3D v_v_sface_ = *p_v_v_sface;
//   const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
//   const ViewDouble4D v_ws_      = *p_v_ws;
//   const ViewDouble5D v_at_      = *p_v_at;
// };
// End advection tracer


// smts(vtl, vit, fil_lat)
class functor_tracer_smts_1 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j) const {
    v_nn_(j) = 0;
    return ;
  }
 private:
  const ViewInt1D v_nn_ = *p_v_nn;
};

class functor_tracer_smts_2 {
 public:
  functor_tracer_smts_2(const double &fil_lat) 
    : fil_lat_(fil_lat) {}

  KOKKOS_INLINE_FUNCTION void
      operator() (const int &j) const {
    const int max_nn = 10;
    if (cos(v_tlat_(0, j, 0)) <= cos(DEGtoRAD(fil_lat_))) {
      v_nn_(j) = std::min(max_nn,
          abs(static_cast<int>(cos(DEGtoRAD(fil_lat_)) 
              / cos(v_tlat_(0, j, 0)) * 1.2)));
    }
    if (v_nn_(j) < 0) {
      v_nn_(j) = 0;
    }
    return ;
  }
 private:
  const double fil_lat_;

  const ViewInt1D v_nn_ = *p_v_nn;

  const ViewDouble3D v_tlat_ = *p_v_tlat;

  KOKKOS_INLINE_FUNCTION 
      double DEGtoRAD(const double &degree) const {
    return (degree * (PI / 180));
  }
};

class functor_tracer_smts_3 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j) const {

    const int iblock = 0;
    const int max_nn = 10;

    for (int ncy = 0; ncy < max_nn; ++ncy) {
      if (v_nn_(j) >= ncy) {
        for (int i = 0; i < IMT; ++i) {
          v_xs_(i) = v_vtl_(iblock, k, j, i) * v_vit_(iblock, k, j, i);
        }
        for (int i = 2; i < IMT-2; ++i) {
          v_vtl_(iblock, k, j, i) = (v_xs_(i) * (1.0 
              - 0.25 * v_vit_(iblock, k, j, i-1)
              - 0.25 * v_vit_(iblock, k, j, i+1))
              + 0.25 * (v_xs_(i-1) + v_xs_(i+1))) 
                  * v_vit_(iblock, k, j, i);
        }
      }
    }
    return ;
  }
 private:
  const ViewInt1D    v_nn_  = *p_v_nn;
  const ViewDouble1D v_xs_  = *p_v_xs;

  const ViewDouble4D v_vit_ = *p_v_vit;
  const ViewDouble4D v_vtl_ = *p_v_vtl;
};

// End smts

#else  // NODIAG

class functor_tracer_39 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    v_atb_(iblock, 0, k, j, i) = 12.0;
    v_atb_(iblock, 1, k, j, i) = C0;
    return ;
  }

 private:
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class functor_tracer_40 {
 public:
  KOKKOS_INLINE_FUNCTION void
      operator() (const int &k, const int &j, const int &i) 
          const {
    const int iblock = 0;

    for (int n = 0; n < NTRA; ++n) {
      v_at_(iblock, n, k, j, i) = v_atb_(iblock, n, k+1, j, i);
    }
    return ;
  }

 private:
  const ViewDouble5D v_at_  = *p_v_at;
  const ViewDouble5D v_atb_ = *p_v_atb;
};
#endif // NODIAG


KOKKOS_REGISTER_FOR_2D(FunctorTracer1, FunctorTracer1)
KOKKOS_REGISTER_FOR_3D(FunctorTracer2, FunctorTracer2)
KOKKOS_REGISTER_FOR_2D(FunctorTracer3, FunctorTracer3)
KOKKOS_REGISTER_FOR_3D(FunctorTracer4, FunctorTracer4)
KOKKOS_REGISTER_FOR_3D(FunctorTracer5, FunctorTracer5)
KOKKOS_REGISTER_FOR_2D(FunctorTracer6, FunctorTracer6)
KOKKOS_REGISTER_FOR_3D(FunctorTracer7, FunctorTracer7)
KOKKOS_REGISTER_FOR_3D(FunctorTracer8, FunctorTracer8)
KOKKOS_REGISTER_FOR_3D(FuncAdvTraCenTsp1, FuncAdvTraCenTsp1)
KOKKOS_REGISTER_FOR_3D(FuncAdvTraTsp2, FuncAdvTraTsp2)
KOKKOS_REGISTER_FOR_3D(FuncAdvTraTsp3, FuncAdvTraTsp3)
KOKKOS_REGISTER_FOR_2D(FuncAdvTraTsp4, FuncAdvTraTsp4)
KOKKOS_REGISTER_FOR_3D(FuncAdvTraTsp5, FuncAdvTraTsp5)
KOKKOS_REGISTER_FOR_2D(FuncAdvTraTsp6, FuncAdvTraTsp6)
KOKKOS_REGISTER_FOR_3D(FunctorTracer15, FunctorTracer15)
KOKKOS_REGISTER_FOR_3D(FunctorTracer16, FunctorTracer16)
KOKKOS_REGISTER_FOR_3D(FunctorTracer17, FunctorTracer17)
KOKKOS_REGISTER_FOR_3D(FunctorTracer17Del2, FunctorTracer17Del2)
KOKKOS_REGISTER_FOR_3D(FunctorTracer18, FunctorTracer18)
KOKKOS_REGISTER_FOR_2D(FunctorTracer19, FunctorTracer19)
KOKKOS_REGISTER_FOR_3D(FunctorTracer20, FunctorTracer20)
KOKKOS_REGISTER_FOR_2D(FunctorTracer21, FunctorTracer21)
KOKKOS_REGISTER_FOR_2D(FunctorTracer22, FunctorTracer22)
KOKKOS_REGISTER_FOR_3D(FunctorTracer23, FunctorTracer23)
KOKKOS_REGISTER_FOR_2D(FunctorTracer24, FunctorTracer24)
KOKKOS_REGISTER_FOR_2D(FunctorTracer25, FunctorTracer25)
KOKKOS_REGISTER_REDUCE_2D(FunctorTracer26, FunctorTracer26)
KOKKOS_REGISTER_FOR_2D(FunctorTracer27, FunctorTracer27)
KOKKOS_REGISTER_FOR_2D(FunctorTracer28, FunctorTracer28)
KOKKOS_REGISTER_FOR_2D(FunctorTracer29, FunctorTracer29)
KOKKOS_REGISTER_FOR_3D(FunctorTracer30, FunctorTracer30)
// KOKKOS_REGISTER_FOR_3D(FunctorTracer31, FunctorTracer31)
KOKKOS_REGISTER_FOR_2D(FunctorTracer32, FunctorTracer32)
// KOKKOS_REGISTER_FOR_2D(FunctorTracer33, FunctorTracer33)
// KOKKOS_REGISTER_FOR_3D(FunctorTracer34, FunctorTracer34)
KOKKOS_REGISTER_FOR_3D(FunctorTracer35, FunctorTracer35)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_TRACER_HPP_
