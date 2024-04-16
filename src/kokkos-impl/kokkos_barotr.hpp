#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BAROTR_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BAROTR_HPP_
#include "../head/def-undef.h"

#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_extern_functions.h"
#ifndef BIHAR
#include "../head/cpp_hmix_del2.h"
#else  // BIHAR
#include "../head/cpp_hmix_del4.h"
#endif // BIHAR
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/fortran_extern_functions.h"

#include "../head/cpp_extern_functions.h"

#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_grid.h"

#ifndef BIHAR
#include "../head/kokkos_hmix_del2.h"
#else  // BIHAR
#include "../head/kokkos_hmix_del4.h"
#endif // BIHAR

#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_work_mod.h"
#include "../head/kokkos_tmp_var.h"

#include "../head/kokkos_config.hpp"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include "Kokkos_Core.hpp"

using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::JST;
using CppParamMod::JET;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;

using CppConstantMod::G;
using CppConstantMod::C0;
using CppConstantMod::P5;
using CppConstantMod::P25;

using KokkosDynMod   ::p_v_h0;
using KokkosDynMod   ::p_v_h0f;
using KokkosDynMod   ::p_v_h0bf;
using KokkosDynMod   ::p_v_h0p;
using KokkosDynMod   ::p_v_vb;
using KokkosDynMod   ::p_v_ub;
using KokkosDynMod   ::p_v_vbp;
using KokkosDynMod   ::p_v_ubp;
using KokkosDynMod   ::p_v_dlub;
using KokkosDynMod   ::p_v_dlvb;
using KokkosGrid     ::p_v_au0;
using KokkosGrid     ::p_v_aus;
using KokkosGrid     ::p_v_auw;
using KokkosGrid     ::p_v_ausw;
using KokkosGrid     ::p_v_dxu;
using KokkosGrid     ::p_v_dyu;
using KokkosGrid     ::p_v_dxur;
using KokkosGrid     ::p_v_dyur;
using KokkosGrid     ::p_v_dxyur;
using KokkosGrid     ::p_v_fcor;
using KokkosGrid     ::p_v_htw;
using KokkosGrid     ::p_v_hts;
using KokkosGrid     ::p_v_kmt;
using KokkosGrid     ::p_v_kmu;
using KokkosGrid     ::p_v_kmtn;
using KokkosGrid     ::p_v_kmts;
using KokkosGrid     ::p_v_kmte;
using KokkosGrid     ::p_v_kmtw;
using KokkosGrid     ::p_v_kmt_nsew;
using KokkosGrid     ::p_v_uarea;
using KokkosGrid     ::p_v_tarea_r;
using KokkosPconstMod::p_v_vit;
using KokkosPconstMod::p_v_viv;
using KokkosPconstMod::p_v_dzph;
using KokkosPconstMod::p_v_ebea;
using KokkosPconstMod::p_v_ebeb;
using KokkosWorkMod  ::p_v_pax;
using KokkosWorkMod  ::p_v_pay;
using KokkosWorkMod  ::p_v_pxb;
using KokkosWorkMod  ::p_v_pyb;
using KokkosWorkMod  ::p_v_whx;
using KokkosWorkMod  ::p_v_why;
using KokkosWorkMod  ::p_v_wka;
using KokkosWorkMod  ::p_v_wgp;
using KokkosWorkMod  ::p_v_work;
#ifdef BIHAR
using KokkosHmixDel4::  p_v_ahf;
using KokkosHmixDel4::  p_v_amf;
using KokkosHmixDel4::  p_v_duc;
using KokkosHmixDel4::  p_v_dum;
using KokkosHmixDel4::  p_v_dun;
using KokkosHmixDel4::  p_v_dus;
using KokkosHmixDel4::  p_v_due;
using KokkosHmixDel4::  p_v_duw;
using KokkosHmixDel4::  p_v_dmc;
using KokkosHmixDel4::  p_v_dmn;
using KokkosHmixDel4::  p_v_dms;
using KokkosHmixDel4::  p_v_dme;
using KokkosHmixDel4::  p_v_dmw;
using KokkosHmixDel4::  p_v_dtn;
using KokkosHmixDel4::  p_v_dts;
using KokkosHmixDel4::  p_v_dte;
using KokkosHmixDel4::  p_v_dtw;
using KokkosHmixDel4::  p_v_dt_nsew;
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
using KokkosHmixDel2::  p_v_dtn;
using KokkosHmixDel2::  p_v_dts;
using KokkosHmixDel2::  p_v_dte;
using KokkosHmixDel2::  p_v_dtw;
#endif // BIHAR

#ifdef BIHAR

using KokkosTmpVar::p_v_c_cnsew;

using KokkosTmpVar::p_v_curl;
using KokkosTmpVar::p_v_d2uk;
using KokkosTmpVar::p_v_d2vk;
using KokkosTmpVar::p_v_dt2k;
using KokkosTmpVar::p_v_gradx;
using KokkosTmpVar::p_v_grady;
#endif // BIHAR
using KokkosTmpVar::p_v_div_out;
// using KokkosTmpVar::p_v_am_factor;

class FunctorBarotr1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    v_wka_(iblock, k, j, i) = 0.0;
    return ;
  }
 private:
  const ViewDouble4D v_wka_ = *p_v_wka;
};

class FunctorBarotr2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    tgrid_to_ugrid(iblock, j, i, v_work_, v_h0_);
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
      v_ugrid(iblock, j, i) = 0.0;
    }
    return ;
  }
 private:
  // const ViewDouble3D v_au0_  = *p_v_au0;
  // const ViewDouble3D v_aus_  = *p_v_aus;
  // const ViewDouble3D v_auw_  = *p_v_auw;
  // const ViewDouble3D v_ausw_ = *p_v_ausw;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_work_ = *p_v_work;
};

class FunctorBarotr3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    // if (i < (IMT-1) && j >= 1) {
      v_wka_(iblock, 0, j, i) = v_ub_(iblock, j, i)
          * (v_dzph_(iblock, j, i) + v_work_(iblock, j, i));
      v_wka_(iblock, 1, j, i) = v_vb_(iblock, j, i)
          * (v_dzph_(iblock, j, i) + v_work_(iblock, j, i));
    // }
    return ;
  }
 private:
  const ViewDouble3D v_ub_   = *p_v_ub;
  const ViewDouble3D v_vb_   = *p_v_vb;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble3D v_dzph_ = *p_v_dzph;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

class FunctorBarotr4 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    div(iblock, 0, j, i, v_div_out_, v_wka_, v_wka_);
    v_work_(iblock, j, i) = - v_vit_(iblock, 0, j, i)
        * v_div_out_(j, i) * P25; 
    return ;
  }
  KOKKOS_INLINE_FUNCTION void div (const int &iblock, const int &k, const int &j, 
      const int &i, const ViewDouble2D &v_div_out, const ViewDouble4D &v_ux, 
          const ViewDouble4D &v_uy) const {
    v_div_out(j, i) = C0;
    if (i < (NX_BLOCK-1) && j >= 1) {
      const int bid = 0;
      if (k <= (v_kmt_(bid, j, i) - 1)) {
        v_div_out(j, i) = P5 * (
            (v_ux(iblock, 0, j  , i+1) + v_ux(iblock, 0, j-1, i+1)) * v_htw_(bid, j  , i+1)
          - (v_ux(iblock, 0, j  , i  ) + v_ux(iblock, 0, j-1, i  )) * v_htw_(bid, j  , i  )
          + (v_uy(iblock, 1, j  , i+1) + v_uy(iblock, 1, j  , i  )) * v_hts_(bid, j  , i  )
          - (v_uy(iblock, 1, j-1, i+1) + v_uy(iblock, 1, j-1, i  )) * v_hts_(bid, j-1, i  )) 
              * v_tarea_r_(bid, j, i);
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_work_    = *p_v_work;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_wka_     = *p_v_wka;
  const ViewDouble4D v_vit_     = *p_v_vit;
};

class FunctorBarotr5 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_h0_(iblock, j, i) = v_h0p_(iblock, j, i) + v_work_(iblock, j, i) * dtb_;
    return;
  }
 private:
  const double dtb_ = CppPconstMod::dtb;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_h0p_  = *p_v_h0p;
  const ViewDouble3D v_work_ = *p_v_work;
};

class FunctorBarotr6 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    hdiffu_del4_1 (iblock, 0, j, i, v_ubp_, v_vbp_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_1(
      const int &iblock, const int &k, const int &j, const int &i,
          const ViewDouble3D &v_umixk, const ViewDouble3D &v_vmixk)
              const {
    div   (j, i, v_div_out_, v_umixk, v_vmixk);
    zcurl (j, i, v_curl_,    v_umixk, v_vmixk);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void div (const int &j, const int &i, 
      const ViewDouble2D &v_div_out, const ViewDouble3D &v_ux, const ViewDouble3D &v_uy) const {
    v_div_out(j, i) = C0;
    if (i < (NX_BLOCK-1) && j >= 1) {
      const int bid = 0;
      // const int k = 0;
      if (0 <= (v_kmt_(bid, j, i) - 1)) {
        v_div_out(j, i) = P5 * (
            (v_ux(bid, j  , i+1) + v_ux(bid, j-1, i+1)) * v_htw_(bid, j  , i+1)
          - (v_ux(bid, j  , i  ) + v_ux(bid, j-1, i  )) * v_htw_(bid, j  , i  )
          + (v_uy(bid, j  , i+1) + v_uy(bid, j  , i  )) * v_hts_(bid, j  , i  )
          - (v_uy(bid, j-1, i+1) + v_uy(bid, j-1, i  )) * v_hts_(bid, j-1, i  )) 
              * v_tarea_r_(bid, j, i);
      }
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void zcurl (const int &j, const int &i,
      const ViewDouble2D &v_curl, const ViewDouble3D &v_ux, const ViewDouble3D &v_uy) const {
    v_curl(j, i) = C0;
    if (i >= 1 && j >= 1) {
      const int bid = 0;
      // const int k = 0;
      if (0 <= v_kmt_(bid, j, i) - 1) {
        v_curl(j, i) = P5 * (
            v_uy(bid, j  , i  ) * v_dyu_(bid, j  , i  )
          + v_uy(bid, j-1, i  ) * v_dyu_(bid, j-1, i  )
          - v_uy(bid, j  , i-1) * v_dyu_(bid, j  , i-1)
          - v_uy(bid, j-1, i-1) * v_dyu_(bid, j-1, i-1)
          - v_ux(bid, j  , i  ) * v_dxu_(bid, j  , i  )
          - v_ux(bid, j  , i-1) * v_dxu_(bid, j  , i-1)
          + v_ux(bid, j-1, i  ) * v_dxu_(bid, j-1, i  )
          + v_ux(bid, j-1, i-1) * v_dxu_(bid, j-1, i-1));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble2D v_curl_    = *p_v_curl;
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_dxu_     = *p_v_dxu;
  const ViewDouble3D v_dyu_     = *p_v_dyu;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble3D v_ubp_     = *p_v_ubp;
  const ViewDouble3D v_vbp_     = *p_v_vbp;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
};

class FunctorBarotr7 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    hdiffu_del4_2 (j, i, v_ubp_, v_vbp_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_2 (
      const int &j, const int &i,
          const ViewDouble3D &v_umixk, const ViewDouble3D &v_vmixk)
              const {
    const int bid = 0;

    if (i >= (ib_-2) && i < (ie_+1) && j >= (jb_-2) && j < (je_+1)) {
      const double cc = v_du_cnsewm_(bid, j, i, 0) + v_du_cnsewm_(bid, j, i, 5);
      v_d2uk_(0, j, i) = (           cc * v_umixk(bid, j  , i  )
           + v_du_cnsewm_(bid, j, i, 1) * v_umixk(bid, j-1, i  )
           + v_du_cnsewm_(bid, j, i, 2) * v_umixk(bid, j+1, i  )
           + v_du_cnsewm_(bid, j, i, 3) * v_umixk(bid, j  , i+1)
           + v_du_cnsewm_(bid, j, i, 4) * v_umixk(bid, j  , i-1))
          + (v_dm_cnsew_ (bid, j, i, 0) * v_vmixk(bid, j  , i  )
           + v_dm_cnsew_ (bid, j, i, 1) * v_vmixk(bid, j-1, i  )
           + v_dm_cnsew_ (bid, j, i, 2) * v_vmixk(bid, j+1, i  )
           + v_dm_cnsew_ (bid, j, i, 3) * v_vmixk(bid, j  , i+1)
           + v_dm_cnsew_ (bid, j, i, 4) * v_vmixk(bid, j  , i-1));
      
      v_d2vk_(0, j, i) = (           cc * v_vmixk(bid, j  , i  )
           + v_du_cnsewm_(bid, j, i, 1) * v_vmixk(bid, j-1, i  )
           + v_du_cnsewm_(bid, j, i, 2) * v_vmixk(bid, j+1, i  )
           + v_du_cnsewm_(bid, j, i, 3) * v_vmixk(bid, j  , i+1)
           + v_du_cnsewm_(bid, j, i, 4) * v_vmixk(bid, j  , i-1))
          + (v_dm_cnsew_ (bid, j, i, 0) * v_umixk(bid, j  , i  )
           + v_dm_cnsew_ (bid, j, i, 1) * v_umixk(bid, j-1, i  )
           + v_dm_cnsew_ (bid, j, i, 2) * v_umixk(bid, j+1, i  )
           + v_dm_cnsew_ (bid, j, i, 3) * v_umixk(bid, j  , i+1)
           + v_dm_cnsew_ (bid, j, i, 4) * v_umixk(bid, j  , i-1));
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
  const ViewDouble3D v_ubp_       = *p_v_ubp;
  const ViewDouble3D v_vbp_       = *p_v_vbp;
  const ViewDouble4D v_dm_cnsew_  = *p_v_dm_cnsew;
  const ViewDouble4D v_du_cnsewm_ = *p_v_du_cnsewm;
};

class FunctorBarotr8 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    hdiffu_del4_3 (j, i);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_3 (const int &j, const int &i) const {
    const int bid = 0;
    double am_factor = 1.0;
    const double amf = v_amf_(bid, j, i);
    if (i >= (ib_-2) && i < (ie_+1) && j >= (jb_-2) && j < (je_+1)) {
      const double sqrt_uarea = sqrt (v_uarea_(bid, j, i));
      const double dxdy = sqrt_uarea * sqrt_uarea * sqrt_uarea * sqrt_uarea * sqrt_uarea * 45.0;
      double gradx1, grady1;
      grad (j, i, gradx1, grady1, v_curl_);
      gradx1 *= gradx1;
      grady1 *= grady1;
      double gradx2, grady2;
      grad (j, i, gradx2, grady2, v_div_out_);
      gradx2 *= gradx2;
      grady2 *= grady2;
      am_factor = sqrt(gradx1 + gradx2 + grady1 + grady2)  
          * dxdy / fabs(am_ * amf);
    }
    am_factor = fmin(40.0, am_factor);
    am_factor = fmax(1.0,  am_factor);
    // const int k = 0;
    if (0 <= v_kmu_(bid, j, i) - 1) {
      v_d2uk_(0, j, i) *= am_factor * amf;
      v_d2vk_(0, j, i) *= am_factor * amf;
    } else {
      v_d2uk_(0, j, i) = C0; 
      v_d2vk_(0, j, i) = C0;
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &j, const int &i, 
      double &gradx, double &grady, const ViewDouble2D &v_f) const {
    const int bid = 0;
    gradx = C0;
    grady = C0;
    if (i >= 1 && j < (NY_BLOCK-1)) {
      // const int k = 0;
      if (0 <= v_kmu_(bid, j, i) - 1) {
        gradx = v_dxyur_(bid, j, i, 0) * P5 * 
            (v_f(j+1, i  ) - v_f(j, i-1) - 
             v_f(j+1, i-1) + v_f(j, i  ));
      
        grady = v_dxyur_(bid, j, i, 1) * P5 * 
            (v_f(j+1, i  ) - v_f(j, i-1) + 
             v_f(j+1, i-1) - v_f(j, i  ));
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
  const ViewInt3D    v_kmu_     = *p_v_kmu;
  const ViewDouble2D v_curl_    = *p_v_curl;
  const ViewDouble3D v_d2uk_    = *p_v_d2uk;
  const ViewDouble3D v_d2vk_    = *p_v_d2vk;
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_amf_     = *p_v_amf;
  const ViewDouble3D v_uarea_   = *p_v_uarea;
  const ViewDouble4D v_dxyur_   = *p_v_dxyur;
};

class FunctorBarotr9 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      double hduk, hdvk;
      hdiffu_del4_4 (j, i, hduk, hdvk);
      v_wka_(0, 4, j, i) += hduk;
      v_wka_(0, 5, j, i) += hdvk;
    }

    return ;
  }
  KOKKOS_INLINE_FUNCTION void hdiffu_del4_4 (const int &j, const int &i,
      double &hduk, double &hdvk) const {
    const int bid = 0;
    hduk = C0;
    hdvk = C0;
    if (i >= (ib_-1) && i < (ie_) && j >= (jb_-1) && j < (je_)) {
      const double cc = v_du_cnsewm_(bid, j, i, 0) 
          + v_du_cnsewm_(bid, j, i, 5);
      hduk = am_ * ((               cc * v_d2uk_(0, j  , i  )
          + v_du_cnsewm_(bid, j, i, 1) * v_d2uk_(0, j-1, i  )
          + v_du_cnsewm_(bid, j, i, 2) * v_d2uk_(0, j+1, i  )
          + v_du_cnsewm_(bid, j, i, 3) * v_d2uk_(0, j  , i+1)
          + v_du_cnsewm_(bid, j, i, 4) * v_d2uk_(0, j  , i-1))
          + (v_dm_cnsew_(bid, j, i, 0) * v_d2vk_(0, j  , i  )
          +  v_dm_cnsew_(bid, j, i, 1) * v_d2vk_(0, j-1, i  )
          +  v_dm_cnsew_(bid, j, i, 2) * v_d2vk_(0, j+1, i  )
          +  v_dm_cnsew_(bid, j, i, 3) * v_d2vk_(0, j  , i+1)
          +  v_dm_cnsew_(bid, j, i, 4) * v_d2vk_(0, j  , i-1)));
   
      hdvk = am_ * ((               cc * v_d2vk_(0, j  , i  )
          + v_du_cnsewm_(bid, j, i, 1) * v_d2vk_(0, j-1, i  )
          + v_du_cnsewm_(bid, j, i, 2) * v_d2vk_(0, j+1, i  )
          + v_du_cnsewm_(bid, j, i, 3) * v_d2vk_(0, j  , i+1)
          + v_du_cnsewm_(bid, j, i, 4) * v_d2vk_(0, j  , i-1))
          + (v_dm_cnsew_(bid, j, i, 0) * v_d2uk_(0, j  , i  )
          +  v_dm_cnsew_(bid, j, i, 1) * v_d2uk_(0, j-1, i  )
          +  v_dm_cnsew_(bid, j, i, 2) * v_d2uk_(0, j+1, i  )
          +  v_dm_cnsew_(bid, j, i, 3) * v_d2uk_(0, j  , i+1)
          +  v_dm_cnsew_(bid, j, i, 4) * v_d2uk_(0, j  , i-1)));
    }
    // const int k = 0;
    if (0 > v_kmu_(bid, j, i) - 1) {
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
  const ViewDouble4D v_wka_       = *p_v_wka;
  const ViewDouble4D v_dm_cnsew_  = *p_v_dm_cnsew;
  const ViewDouble4D v_du_cnsewm_ = *p_v_du_cnsewm;
};

class FunctorBarotr11 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    grad(iblock, 0, j, i, v_gradx_, v_grady_, v_h0_);

    if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      const double gstar = (v_wgp_(iblock, j, i) - 1.0) * G;
      v_wka_(iblock, 0, j, i) = v_wka_(iblock, 4, j, i)
          + gstar * v_gradx_(j, i);
      v_wka_(iblock, 1, j, i) = v_wka_(iblock, 5, j, i)
          + gstar * v_grady_(j, i);
    }
    return ;
  }
  KOKKOS_INLINE_FUNCTION void grad (const int &iblock, const int &k,
      const int &j, const int &i, const ViewDouble2D &v_gradx,
          const ViewDouble2D &v_grady, const ViewDouble3D &v_f) 
              const {
    const int bid = 0;
    v_gradx(j, i) = C0;
    v_grady(j, i) = C0;
    if (i >= 1 && j < (NY_BLOCK-1)) {
      if (k <= v_kmu_(bid, j, i) - 1) {
        v_gradx(j, i) = v_dxyur_(bid, j, i, 0) * P5 * 
            (v_f(iblock, j+1, i  ) - v_f(iblock, j, i-1) - 
             v_f(iblock, j+1, i-1) + v_f(iblock, j, i  ));
      
        v_grady(j, i) = v_dxyur_(bid, j, i, 1) * P5 * 
            (v_f(iblock, j+1, i  ) - v_f(iblock, j, i-1) + 
             v_f(iblock, j+1, i-1) - v_f(iblock, j, i  ));
      }
    }
    return ;
  }
 private:
  const ViewInt3D    v_kmu_   = *p_v_kmu;
  const ViewDouble2D v_gradx_ = *p_v_gradx;
  const ViewDouble2D v_grady_ = *p_v_grady;
  const ViewDouble3D v_h0_    = *p_v_h0;
  const ViewDouble3D v_wgp_   = *p_v_wgp;
  const ViewDouble4D v_wka_   = *p_v_wka;
  const ViewDouble4D v_dxyur_ = *p_v_dxyur;
};

class FunctorBarotr12 {
 public:
  KOKKOS_INLINE_FUNCTION void operator() (
      const int &j, const int &i) const {
    const int iblock = 0;
    tgrid_to_ugrid(iblock, j, i, v_work_, v_h0_);
    return ;
  }
  KOKKOS_INLINE_FUNCTION void tgrid_to_ugrid (
      const int &iblock, const int &j, const int &i,
          const ViewDouble3D &v_ugrid, const ViewDouble3D &v_tgrid) const {
    if (i >= 1 && j < (NY_BLOCK-1)) {
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

class FunctorBarotr13 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) 
      const {
    const int iblock = 0;

    // if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      v_wka_(iblock, 0, j, i) = v_viv_(iblock, 0, j, i)
          * (v_wka_(iblock, 0, j, i) + v_dlub_(iblock, j, i)
              - v_fcor_(iblock, j, i) * v_vbp_(iblock, j, i)
                  + v_pax_(iblock, j, i) + v_pxb_(iblock, j, i)
                      - v_work_(iblock, j, i) * v_whx_(iblock, j, i));
    // }
    return ;
  }
 private:
  const ViewDouble3D v_pax_  = *p_v_pax;
  const ViewDouble3D v_pxb_  = *p_v_pxb;
  const ViewDouble3D v_vbp_  = *p_v_vbp;
  const ViewDouble3D v_whx_  = *p_v_whx;
  const ViewDouble3D v_fcor_ = *p_v_fcor;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble3D v_dlub_ = *p_v_dlub;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_viv_  = *p_v_viv;
};
class FunctorBarotr14 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;

    // if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {
      v_wka_(iblock, 1, j, i) = v_viv_(iblock, 0, j, i)
          * (v_wka_(iblock, 1, j, i) + v_dlvb_(iblock, j, i)
              + v_fcor_(iblock, j, i) * v_ubp_(iblock, j, i)
                  + v_pay_(iblock, j, i) + v_pyb_(iblock, j, i)
                      - v_work_(iblock, j, i) * v_why_(iblock, j, i));
    // }
    return ;
  }
 private:
  const ViewDouble3D v_pay_  = *p_v_pay;
  const ViewDouble3D v_pyb_  = *p_v_pyb;
  const ViewDouble3D v_ubp_  = *p_v_ubp;
  const ViewDouble3D v_why_  = *p_v_why;
  const ViewDouble3D v_fcor_ = *p_v_fcor;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble3D v_dlvb_ = *p_v_dlvb;
  const ViewDouble4D v_wka_  = *p_v_wka;
  const ViewDouble4D v_viv_  = *p_v_viv;
};

class FunctorBarotr15 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;

    // if (i >= 2 && i < (IMT-2) && j >= 2 && j < (JMT-2)) {

      v_wka_(iblock, 2, j ,i) = v_ebea_(iblock, j, i) * v_wka_(iblock, 0, j ,i)
                              - v_ebeb_(iblock, j, i) * v_wka_(iblock, 1, j ,i);
      v_wka_(iblock, 3, j ,i) = v_ebea_(iblock, j, i) * v_wka_(iblock, 1, j ,i)
                              + v_ebeb_(iblock, j, i) * v_wka_(iblock, 0, j ,i);
    // }
    return ;
  }
 private:
  const ViewDouble3D v_ebea_ = *p_v_ebea;
  const ViewDouble3D v_ebeb_ = *p_v_ebeb;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

class FunctorBarotr16 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_ub_(iblock, j, i) = v_ubp_(iblock, j, i)
        + v_wka_(iblock, 2, j ,i) * dtb_;
    v_vb_(iblock, j, i) = v_vbp_(iblock, j, i)
        + v_wka_(iblock, 3, j ,i) * dtb_;
    
    if (i < (IMT-1) && j >= 1) {
      v_wka_(iblock, 0, j ,i) = v_ub_(iblock, j, i)
          * (v_dzph_(iblock, j, i) + v_work_(iblock, j, i));
      v_wka_(iblock, 1, j ,i) = v_vb_(iblock, j, i)
          * (v_dzph_(iblock, j, i) + v_work_(iblock, j, i));
    }
    return;
  }

 private:
  const double dtb_ = CppPconstMod::dtb;
  const ViewDouble3D v_ub_   = *p_v_ub;
  const ViewDouble3D v_vb_   = *p_v_vb;
  const ViewDouble3D v_ubp_  = *p_v_ubp;
  const ViewDouble3D v_vbp_  = *p_v_vbp;
  const ViewDouble3D v_dzph_ = *p_v_dzph;
  const ViewDouble3D v_work_ = *p_v_work;
  const ViewDouble4D v_wka_  = *p_v_wka;
};

class FunctorBarotr17 {
 public:
  FunctorBarotr17 (const int &iblock) : iblock_(iblock) {}
  KOKKOS_INLINE_FUNCTION void operator() (
      const int &j, const int &i) const {
    div(iblock_, 0, j, i, v_div_out_, v_wka_, v_wka_);
    return;
  }
  KOKKOS_INLINE_FUNCTION void div (const int &iblock, const int &k, const int &j,
      const int &i, const ViewDouble2D &v_div_out, const ViewDouble4D &v_ux,
          const ViewDouble4D &v_uy) const {
    v_div_out(j, i) = C0;
    if (i < (NX_BLOCK - 1) && j >= 1) {
      const int bid = 0;
      if (k <= (v_kmt_(bid, j, i) - 1)) {
        v_div_out(j, i) = P5 * (
            (v_ux(iblock, 0, j  , i+1) + v_ux(iblock, 0, j-1, i+1)) * v_htw_(bid, j  , i+1)
          - (v_ux(iblock, 0, j  , i  ) + v_ux(iblock, 0, j-1, i  )) * v_htw_(bid, j  , i  )
          + (v_uy(iblock, 1, j  , i+1) + v_uy(iblock, 1, j  , i  )) * v_hts_(bid, j  , i  )
          - (v_uy(iblock, 1, j-1, i+1) + v_uy(iblock, 1, j-1, i  )) * v_hts_(bid, j-1, i  ))
              * v_tarea_r_(bid, j, i);
      }
    }
    return ;
  }

 private:
  const int iblock_;
  const ViewInt3D    v_kmt_     = *p_v_kmt;
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_htw_     = *p_v_htw;
  const ViewDouble3D v_hts_     = *p_v_hts;
  const ViewDouble3D v_tarea_r_ = *p_v_tarea_r;
  const ViewDouble4D v_wka_     = *p_v_wka;
};

class FunctorBarotr18 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {
    const int iblock = 0;
    hdifft_del4_1 (iblock, 1, j, i, v_dt2k_, v_h0p_);
    return;
  }
  KOKKOS_INLINE_FUNCTION void hdifft_del4_1 (
      const int &iblock, const int &k, const int &j,
      const int &i, const ViewDouble3D &v_d2tk,
      const ViewDouble3D &v_tmix) const {
    const int bid = 0;
    // n s e w c
    v_c_cnsew_(0, j, i, 1) = (k <= v_kmt_nsew_(bid, j, i, 0) && k <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 0) : C0;
    v_c_cnsew_(0, j, i, 2) = (k <= v_kmt_nsew_(bid, j, i, 1) && k <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 1) : C0;
    v_c_cnsew_(0, j, i, 3) = (k <= v_kmt_nsew_(bid, j, i, 2) && k <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 2) : C0;
    v_c_cnsew_(0, j, i, 4) = (k <= v_kmt_nsew_(bid, j, i, 3) && k <= v_kmt_(bid, j, i))
        ? v_dt_nsew_(bid, j, i, 3) : C0;

    v_c_cnsew_(0, j, i, 0) = -(v_c_cnsew_(0, j, i, 1) + v_c_cnsew_(0, j, i, 2) 
                             + v_c_cnsew_(0, j, i, 3) + v_c_cnsew_(0, j, i, 4));

    if (i >= (ib_ - 2) && i < (ie_ + 1) && j >= (jb_ - 2) && j < (je_ + 1)) {
      v_d2tk(0, j, i) = v_ahf_(bid, j, i) * 
          (v_c_cnsew_(0, j, i, 0) * v_tmix(iblock, j    , i    ) 
         + v_c_cnsew_(0, j, i, 1) * v_tmix(iblock, j - 1, i    ) 
         + v_c_cnsew_(0, j, i, 2) * v_tmix(iblock, j + 1, i    ) 
         + v_c_cnsew_(0, j, i, 3) * v_tmix(iblock, j    , i + 1) 
         + v_c_cnsew_(0, j, i, 4) * v_tmix(iblock, j    , i - 1));
    }
    return;
  }

 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ib;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const ViewInt3D v_kmt_        = *p_v_kmt;
  const ViewInt4D v_kmt_nsew_   = *p_v_kmt_nsew;
  const ViewDouble3D v_dt2k_    = *p_v_dt2k;
  const ViewDouble3D v_ahf_     = *p_v_ahf;
  const ViewDouble3D v_h0p_     = *p_v_h0p;
  const ViewDouble4D v_c_cnsew_ = *p_v_c_cnsew;
  const ViewDouble4D v_dt_nsew_ = *p_v_dt_nsew;
};

class FunctorBarotr19 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;

    double hdtk;

    hdifft_del4_2 (j, i, v_dt2k_, hdtk);

    if (i >= 2 && i < (IMT - 2) && j >= 1 && j < (JMT - 2)) {
      v_work_(iblock, j, i) = v_vit_(iblock, 0, j, i) * (hdtk - v_div_out_(j, i));
    }

    return;
  }
  KOKKOS_INLINE_FUNCTION void hdifft_del4_2 (const int &j, const int &i,
      const ViewDouble3D &v_d2tk, double &hdtk) const {
    hdtk = C0;
    if (i >= (ib_ - 1) && i < ie_ && j >= (jb_ - 1) && j < je_) {
      hdtk = ah_ * (v_c_cnsew_(0, j, i, 0) * v_d2tk(0, j    , i    ) 
                  + v_c_cnsew_(0, j, i, 1) * v_d2tk(0, j - 1, i    ) 
                  + v_c_cnsew_(0, j, i, 2) * v_d2tk(0, j + 1, i    ) 
                  + v_c_cnsew_(0, j, i, 3) * v_d2tk(0, j    , i + 1) 
                  + v_c_cnsew_(0, j, i, 4) * v_d2tk(0, j    , i - 1));
    }
    return;
  }

 private:
  const int ib_ = CppBlocks::ib;
  const int ie_ = CppBlocks::ib;
  const int jb_ = CppBlocks::jb;
  const int je_ = CppBlocks::je;
  const double ah_ = CppHmixDel4::ah;
  const ViewDouble3D v_dt2k_    = *p_v_dt2k;
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_work_    = *p_v_work;
  const ViewDouble4D v_c_cnsew_ = *p_v_c_cnsew;
  const ViewDouble4D v_vit_     = *p_v_vit;
};

class FunctorBarotr20 {
 public:
  KOKKOS_INLINE_FUNCTION void operator() (
      const int &j, const int &i) const {
    const int iblock = 0;
    // const double hdtk = C0;
    // v_work_(iblock, j, i) = v_vit_(iblock, 0, j, i) * (hdtk - v_div_out_(j, i));
    v_work_(iblock, j, i) = - v_vit_(iblock, 0, j, i) * v_div_out_(j, i);
    return;
  }

 private:
  const ViewDouble2D v_div_out_ = *p_v_div_out;
  const ViewDouble3D v_work_    = *p_v_work;
  const ViewDouble4D v_vit_     = *p_v_vit;
};

class FunctorBarotr21 {
public:
  KOKKOS_INLINE_FUNCTION void operator() (
      const int &j, const int &i) const {
    const int iblock = 0;
    v_h0_(iblock, j, i) = v_h0p_(iblock, j, i) + v_work_(iblock, j, i) * dtb_;

    v_ubp_(iblock, j, i) = v_ub_(iblock, j, i);
    v_vbp_(iblock, j, i) = v_vb_(iblock, j, i);
    v_h0p_(iblock, j, i) = v_h0_(iblock, j, i);

    if (j >= (JST - 1) && j < JET) {
      v_h0f_(iblock, j, i)  += v_h0_(iblock, j, i);
      v_h0bf_(iblock, j, i) += v_h0_(iblock, j, i);
    }
    return;
  }
 private:
  const double dtb_ = CppPconstMod::dtb;
  const ViewDouble3D v_ub_   = *p_v_ub;
  const ViewDouble3D v_vb_   = *p_v_vb;
  const ViewDouble3D v_ubp_  = *p_v_ubp;
  const ViewDouble3D v_vbp_  = *p_v_vbp;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_h0p_  = *p_v_h0p;
  const ViewDouble3D v_h0f_  = *p_v_h0f;
  const ViewDouble3D v_h0bf_ = *p_v_h0bf;
  const ViewDouble3D v_work_ = *p_v_work;
};
//===========================

// class functor_barotr_6 {
//  public:
//   functor_barotr_6 (const int &iblock) : iblock_(iblock) {}

//   KOKKOS_INLINE_FUNCTION void
//       operator() (const int &j, const int &i) const {
//     v_wka_(i, j, 4, iblock_) = v_hduk_(j, i);
//     v_wka_(i, j, 5, iblock_) = v_hdvk_(j, i);
//     return;
//   }

//  private:
//   const int iblock_;
//   const ViewDouble2D v_hduk_ = *p_v_hduk;
//   const ViewDouble2D v_hdvk_ = *p_v_hdvk;
//   const ViewDouble4D v_wka_ = *p_v_wka;
// };

// class functor_barotr_13 {
//  public:
//   functor_barotr_13 (const int &iblock) : iblock_(iblock) {}
//   KOKKOS_INLINE_FUNCTION void operator() (const int &j, const int &i) const {
//     v_work_(iblock_, j, i) = v_vit_(iblock_, 0, j, i) * (v_hdtk_(j, i) * 1.0 - v_div_out_(j, i));
//     return;
//   }

//  private:
//   const int iblock_;
//   const ViewDouble2D v_hdtk_ = *p_v_hdtk;
//   const ViewDouble2D v_div_out_ = *p_v_div_out;
//   const ViewDouble3D v_work_ = *p_v_work;
//   const ViewDouble4D v_vit_ = *p_v_vit;
// };

#ifdef SMAG1
#else // SMAG1
#ifdef BIHAR
// hdiffu_del4
#else  // BIHAR
// hdiffu_del2
class functor_barotr_hdiffu_del2_1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(const int &j, const int &i) const {
    v_hduk_(j, i) = C0;
    v_hdvk_(j, i) = C0;
    return;
  };

 private:
  const ViewDouble2D v_hduk_ = *p_v_hduk;
  const ViewDouble2D v_hdvk_ = *p_v_hduk;
};
class functor_barotr_hdiffu_del2_2 {
 public:
  functor_barotr_hdiffu_del2_2(const int &k, const int &iblock)
      : k_(k), iblock_(iblock), v_hduk_(v_hduk), v_hdvk_(v_hdvk) {}
  KOKKOS_INLINE_FUNCTION void operator()(const int &j, const int &i) const {
    const int bid = 0;
    const double cc = v_duc_(bid, j, i) + v_dum_(bid, j, i);
    v_hduk_(j, i) = am_ * ((cc * v_ubp_(iblock, j, i_) + v_dun_(bid, j, i) * v_ubp_(i, j - 1, iblock_) + v_dus_(bid, j, i) * v_ubp_(i, j + 1, iblock_) + v_due_(bid, j, i) * v_ubp_(i + 1, j, iblock_) + v_duw_(bid, j, i) * v_ubp_(i - 1, j, iblock_)) + (v_dmc_(bid, j, i) * v_vbp_(iblock, j, i_) + v_dmn_(bid, j, i) * v_vbp_(i, j - 1, iblock_) + v_dms_(bid, j, i) * v_vbp_(i, j + 1, iblock_) + v_dme_(bid, j, i) * v_vbp_(i + 1, j, iblock_) + v_dmw_(bid, j, i) * v_vbp_(i - 1, j, iblock_)) * v_viv_(i, j, k_, bid));
    v_hdvk_(j, i) = am_ * ((cc * v_vp_(iblock, j, i_) + v_dun_(bid, j, i) * v_vp_(i, j - 1, iblock_) + v_dus_(bid, j, i) * v_vp_(i, j + 1, iblock_) + v_due_(bid, j, i) * v_vp_(i + 1, j, iblock_) + v_duw_(bid, j, i) * v_vp_(i - 1, j, iblock_)) + (v_dmc_(bid, j, i) * v_up_(iblock, j, i_) + v_dmn_(bid, j, i) * v_up_(i, j - 1, iblock_) + v_dms_(bid, j, i) * v_up_(i, j + 1, iblock_) + v_dme_(bid, j, i) * v_up_(i + 1, j, iblock_) + v_dmw_(bid, j, i) * v_up_(i - 1, j, iblock_)) * v_viv_(i, j, k_, bid));
    return;
  };

private:
  const int k_, iblock_;
  const double am_ = CppHmixDel2::am;
  const ViewDouble2D v_hduk_ = *p_v_hduk;
  const ViewDouble2D v_hdvk_ = *p_v_hdvk;
  const ViewDouble3D v_dun_ = *p_v_dun;
  const ViewDouble3D v_dus_ = *p_v_dus;
  const ViewDouble3D v_due_ = *p_v_due;
  const ViewDouble3D v_duw_ = *p_v_duw;
  const ViewDouble3D v_dmc_ = *p_v_dmc;
  const ViewDouble3D v_dmn_ = *p_v_dmn;
  const ViewDouble3D v_dms_ = *p_v_dms;
  const ViewDouble3D v_dme_ = *p_v_dme;
  const ViewDouble3D v_dmw_ = *p_v_dmw;
  const ViewDouble3D v_ubp_ = *p_v_ubp;
  const ViewDouble3D v_vbp_ = *p_v_vbp;
  const ViewDouble4D v_viv_ = *p_v_viv;
};
class functor_barotr_hdiffu_del2_3 {
 public:
  functor_barotr_hdiffu_del2_3(const int &k) : k_(k) {}
  KOKKOS_INLINE_FUNCTION void operator()(const int &j, const int &i) const {
    const int bid = 0;
    if (k_ > v_kmu_(bid, j, i) - 1) {
      v_hduk_(j, i) = C0;
      v_hdvk_(j, i) = C0;
    }
    return ;
  }

 private:
  const ViewDouble2D v_hduk_ = *p_v_hduk;
  const ViewDouble2D v_hdvk_ = *p_v_hdvk;
};
  // hdiffu_del2
  // hdifft_del2
class functor_barotr_hdifft_del2_1 {
 public:
  functor_barotr_hdifft_del2_1(const int &k) : k_(k) {}
  KOKKOS_INLINE_FUNCTION void operator()(const int &j, const int &i) const {
    const int bid = 0;
    v_cn_(j, i) = (k_ <= v_kmtn_(bid, j, i) && k_ <= v_kmt_(bid, j, i))
                      ? v_dtn_(bid, j, i) : C0;
    v_cs_(j, i) = (k_ <= v_kmts_(bid, j, i) && k_ <= v_kmt_(bid, j, i))
                      ? v_dts_(bid, j, i) : C0;
    v_ce_(j, i) = (k_ <= v_kmte_(bid, j, i) && k_ <= v_kmt_(bid, j, i))
                      ? v_dte_(bid, j, i) : C0;
    v_cw_(j, i) = (k_ <= v_kmtw_(bid, j, i) && k_ <= v_kmt_(bid, j, i))
                      ? v_dtw_(bid, j, i) : C0;
    v_cc_(j, i) = -(v_cn_(j, i) + v_cs_(j, i) + v_ce_(j, i) + v_cw_(j, i));
    return;
  }

 private:
    const int k_;
    const ViewDouble2D v_cc_ = *p_v_cc;
    const ViewDouble2D v_cn_ = *p_v_cn;
    const ViewDouble2D v_cs_ = *p_v_cs;
    const ViewDouble2D v_ce_ = *p_v_ce;
    const ViewDouble2D v_cw_ = *p_v_cw;
    const ViewDouble2D v_hdtk_ = *p_v_hdtk;
    const ViewDouble3D v_dtn_ = *p_v_dtn;
    const ViewDouble3D v_dts_ = *p_v_dts;
    const ViewDouble3D v_dte_ = *p_v_dte;
    const ViewDouble3D v_dtw_ = *p_v_dtw;
    const ViewInt3D v_kmt_ = *p_v_kmt;
    const ViewInt3D v_kmtn_ = *p_v_kmtn;
    const ViewInt3D v_kmts_ = *p_v_kmts;
    const ViewInt3D v_kmte_ = *p_v_kmte;
    const ViewInt3D v_kmtw_ = *p_v_kmtw;
  };

  class functor_barotr_hdifft_del2_2 {
  public:
    KOKKOS_INLINE_FUNCTION void operator()(const int &j, const int &i) const {
      v_hdtk_(j, i) = C0;
      return;
    }

 private:
    const ViewDouble2D v_hdtk_ = *p_v_hdtk;
  };

class functor_barotr_hdifft_del2_3 {
 public:
  functor_barotr_hdifft_del2_3 (const int &iblock) : iblock_(iblock) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int &j, const int &i) const {
    v_hdtk_(j, i) = ah_ * (v_cc_(j, i) * v_h0p_(i    , j    , iblock_) 
                         + v_cn_(j, i) * v_h0p_(i    , j - 1, iblock_) 
                         + v_cs_(j, i) * v_h0p_(i    , j + 1, iblock_) 
                         + v_ce_(j, i) * v_h0p_(i - 1, j    , iblock_) 
                         + v_cw_(j, i) * v_h0p_(i - 1, j    , iblock_));
    return;
  }

 private:
  const int iblock_;
  const double ah_ = CppHmixDel2::ah;
  const ViewDouble2D v_cc_ = *p_v_cc;
  const ViewDouble2D v_cn_ = *p_v_cn;
  const ViewDouble2D v_cs_ = *p_v_cs;
  const ViewDouble2D v_ce_ = *p_v_ce;
  const ViewDouble2D v_cw_ = *p_v_cw;
  const ViewDouble2D v_hdtk_ = *p_v_hdtk;
  const ViewDouble3D v_h0p_ = *p_v_h0p;
};
// End hdifft_del4
#endif // BIHAR
#endif // SMAG1

//=========================
KOKKOS_REGISTER_FOR_3D(FunctorBarotr1,  FunctorBarotr1)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr2,  FunctorBarotr2)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr3,  FunctorBarotr3)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr4,  FunctorBarotr4)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr5,  FunctorBarotr5)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr6,  FunctorBarotr6)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr7,  FunctorBarotr7)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr8,  FunctorBarotr8)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr9,  FunctorBarotr9)
// KOKKOS_REGISTER_FOR_2D(FunctorBarotr10, FunctorBarotr10)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr11, FunctorBarotr11)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr12, FunctorBarotr12)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr13, FunctorBarotr13)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr14, FunctorBarotr14)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr15, FunctorBarotr15)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr16, FunctorBarotr16)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr17, FunctorBarotr17)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr18, FunctorBarotr18)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr19, FunctorBarotr19)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr20, FunctorBarotr20)
KOKKOS_REGISTER_FOR_2D(FunctorBarotr21, FunctorBarotr21)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_BAROTR_HPP_
