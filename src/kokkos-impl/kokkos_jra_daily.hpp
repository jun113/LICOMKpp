#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_JRA_DAILY_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_JRA_DAILY_HPP_

#include "../head/def-undef.h"

#include "../head/cpp_param_mod.h"
#include "../head/cpp_output_mod.h"

#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_tracer_mod.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

#include <cstdlib>

using CppParamMod::IMM;
using CppParamMod::JSM;
using CppParamMod::JEM;
using CppParamMod::S_IMT;
using CppParamMod::S_JMT;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::JMT_GLOBAL;

class FunctorInterplationNearest {
 public:
  FunctorInterplationNearest(const int& irec,
      const ViewFloat3D &v_source) : irec_(irec),
          v_source_(v_source) {}
  // KOKKOS_INLINE_FUNCTION void operator () (
  //     const int &j, const int &i) const {
  //   const int iwk = S_IMT + 2;
  //   const int jwk = S_JMT + 2;

  //   // both dx and dy are 0.5625
  //   const double dx_dy = 0.5625;

  //   // wx
  //   // j=0; i=0
  //   if (i == 0 && j == 0) {
  //     v_s_wx_(0, 0) = v_s_lon_(0, 0);
  //   }
  //   // j=0; i=1:iwk-1
  //   if (j == 0) {
  //     v_s_wx_(0, i+1) = v_s_lon_(0, i);
  //   }
  //   // j=0; i=iwk-1
  //   if (i == (S_IMT-1) && j == 0) {
  //     v_s_wx_(0, iwk-1) = v_s_lon_(0, iwk-3) + dx_dy;
  //   }
  //   // j=1:jwk-2; i=1:iwk-2
  //   v_s_wx_(j+1, i+1) = v_s_lon_(j, i); 
  //   // j=1:jwk-1; i=0
  //   if (i == 0) {
  //     v_s_wx_(j+1, 0) = v_s_lon_(j, 0) - dx_dy;
  //   }
  //   // j=jwk-1; i=0
  //   if (i == 0 && j == (S_JMT-1)) {
  //     v_s_wx_(jwk-1, 0) = v_s_lon_(jwk-3, 0);
  //   }
  //   // j=1:jwk-1; i=iwk-1
  //   if (i == (S_IMT-1)) {
  //     v_s_wx_(j+1, iwk-1) = v_s_lon_(j, iwk-3) + dx_dy;
  //   }
  //   // j=jwk-1; i:1,iwk-1
  //   if (j == (S_JMT-1)) {
  //     v_s_wx_(jwk-1, i+1) = v_s_lon_(0, i);
  //   }
  //   // j=jwk-1; i=iwk-1
  //   if (i == (S_IMT-1) && j == (S_JMT-1)) {
  //     v_s_wx_(jwk-1, iwk-1) = v_s_lon_(jwk-3, iwk-3) + dx_dy;
  //   }

  //   // wy
  //   // j=0; i=0
  //   if (i == 0 && j == 0) {
  //     v_s_wy_(0, 0) = v_s_lat_(0, 0) + dx_dy;
  //   }
  //   // j=0; i=1:iwk-1
  //   if (j == 0) {
  //     v_s_wy_(0, i+1) = v_s_lat_(0, i) + dx_dy;
  //   }
  //   // j=0; i=iwk-1
  //   if (i == (S_IMT-1) && j == 0) {
  //     v_s_wy_(0, iwk-1) = v_s_lat_(0, iwk-3) + dx_dy;
  //   }
  //   // j=1:jwk-2; i=1:iwk-2
  //   v_s_wy_(j+1, i+1) = v_s_lat_(j, i); 
  //   // j=1:jwk-1; i=0
  //   if (i == 0) {
  //     v_s_wy_(j+1, 0)   = v_s_lat_(j, 0);
  //   }
  //   // j=jwk-1; i=0
  //   if (i == 0 && j == (S_JMT-1)) {
  //     v_s_wy_(jwk-1, 0) = v_s_lat_(jwk-3, 0) - dx_dy;
  //   }
  //   // j=1:jwk-1; i=iwk-1
  //   if (i == (S_IMT-1)) {
  //     v_s_wy_(j+1, iwk-1)   = v_s_lon_(j, iwk-3);
  //   }
  //   // j=jwk-1; i:1,iwk-1
  //   if (j == (S_JMT-1)) {
  //     v_s_wy_(jwk-1, i+1) = v_s_lat_(jwk-3, i) - dx_dy;
  //   }
  //   // j=jwk-1; i=iwk-1
  //   if (i == (S_IMT-1) && j == (S_JMT-1)) {
  //     v_s_wy_(jwk-1, iwk-1) = v_s_lat_(jwk-3, iwk-3) - dx_dy;
  //   }

  //   // work
  //   if (j == 0) {
  //     v_s_work_(0, i) = v_source_(irec_-1, 0, i);
  //   }
  //   if (i == (S_IMT-1) && j == 0) {
  //     v_s_work_(0, S_IMT-2) = v_source_(irec_-1, 0, S_IMT-2);
  //     v_s_work_(0, S_IMT-1) = v_source_(irec_-1, 0, S_IMT-1);
  //   }
  //   // j=1:jwk-2; i=1:iwk-2
  //   v_s_work_(j+1, i+1) = v_source_(irec_-1, j, i);
  //   // j=jwk-1; i=0:iwk-1
  //   if (j == (S_JMT-1)) {
  //     v_s_work_(jwk-1, i) = v_source_(irec_-1, jwk-3, i);
  //   }

  //   if (i == 0) {
  //     v_s_work_(j+1, 0) = v_source_(irec_-1, j, 0);
  //   }
  //   // iwk - 1
  //   if (i == (S_IMT-1)) {
  //     v_s_work_(j+1, iwk-1) = v_source_(irec_-1, j, iwk-3);
  //   }
  //   if (i == (S_IMT-1) && j == (S_JMT-1)) {
  //     v_s_work_(jwk-1, S_IMT-2) = v_source_(irec_-1, jwk-3, S_IMT-2);
  //     v_s_work_(jwk-1, S_IMT-1) = v_source_(irec_-1, jwk-3, S_IMT-1);
  //   }
  //   return ;
  // }
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iwk = S_IMT + 2;
    const int jwk = S_JMT + 2;

    // both dx and dy are 0.5625
    const double dx_dy = 0.5625;

    v_s_wx_(j+1, i+1)   = v_s_lon_(j, i); 
    v_s_wy_(j+1, i+1)   = v_s_lat_(j, i); 
    v_s_work_(j+1, i+1) = v_source_(irec_-1, j, i);

    if (i == 0) {
      v_s_wx_(j+1, 0)   = v_s_lon_(j, 0) - dx_dy;
      v_s_wy_(j+1, 0)   = v_s_lat_(j, 0);
      v_s_work_(j+1, 0) = v_source_(irec_-1, j, 0);
    }
    // iwk - 1
    if (i == (S_IMT-1)) {
      v_s_wx_(j+1, iwk-1)   = v_s_lon_(j, iwk-3) + dx_dy;
      v_s_wy_(j+1, iwk-1)   = v_s_lon_(j, iwk-3);
      v_s_work_(j+1, iwk-1) = v_source_(irec_-1, j, iwk-3);
    }
    if (j == 0) {
      v_s_wx_(0, i+1) = v_s_lon_(0, i);
      v_s_wy_(0, i+1) = v_s_lat_(0, i) + dx_dy;
      v_s_work_(0, i) = v_source_(irec_-1, 0, i);
    }
    // jwk - 1;
    if (j == (S_JMT-1)) {
      v_s_wx_(jwk-1, i+1) = v_s_lon_(0, i);
      v_s_wy_(jwk-1, i+1) = v_s_lat_(jwk-3, i) - dx_dy;
      v_s_work_(jwk-1, i) = v_source_(irec_-1, jwk-3, i);
    }

    if (i == 0 && j == 0) {
      v_s_wx_(0, 0) = v_s_lon_(0, 0);
      v_s_wy_(0, 0) = v_s_lat_(0, 0) + dx_dy;
    }
    if (i == (S_IMT-1) && j == 0) {
      v_s_wx_(0, iwk-1) = v_s_lon_(0, iwk-3) + dx_dy;
      v_s_wy_(0, iwk-1) = v_s_lat_(0, iwk-3) + dx_dy;
      v_s_work_(0, S_IMT-2) = v_source_(irec_-1, 0, S_IMT-2);
      v_s_work_(0, S_IMT-1) = v_source_(irec_-1, 0, S_IMT-1);
    }
    if (i == 0 && j == (S_JMT-1)) {
      v_s_wx_(jwk-1, 0) = v_s_lon_(jwk-3, 0);
      v_s_wy_(jwk-1, 0) = v_s_lat_(jwk-3, 0) - dx_dy;
    }
    if (i == (S_IMT-1) && j == (S_JMT-1)) {
      v_s_wx_(jwk-1, iwk-1) = v_s_lon_(jwk-3, iwk-3) + dx_dy;
      v_s_wy_(jwk-1, iwk-1) = v_s_lat_(jwk-3, iwk-3) - dx_dy;
      v_s_work_(jwk-1, S_IMT-2) = v_source_(irec_-1, jwk-3, S_IMT-2);
      v_s_work_(jwk-1, S_IMT-1) = v_source_(irec_-1, jwk-3, S_IMT-1);
    }
    return ;
  }
 private:
  const int irec_;
  ViewFloat3D  v_source_;
  ViewDouble2D v_s_wx_   = *KokkosForcMod::p_v_s_wx; 
  ViewDouble2D v_s_wy_   = *KokkosForcMod::p_v_s_wy;
  ViewDouble2D v_s_work_ = *KokkosForcMod::p_v_s_work; 
  ViewDouble1D v_buffer_ = *KokkosForcMod::p_v_buffer;
  ViewDouble2D v_s_lon_  = *KokkosPconstMod::p_v_s_lon;
  ViewDouble2D v_s_lat_  = *KokkosPconstMod::p_v_s_lat;
};

class FunctorNear1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int mx = S_IMT + 2;
    const int my = S_JMT + 2;
    v_buffer_(j*IMT_GLOBAL + i) = SPVAL_;
    int ic = INT_MAX;
    int jc = INT_MAX;
    double lon_tmp;
    // ---find adjacent two grids on x direction (ic)
    if (v_lon_o_(j, i) < 0.0) {
      lon_tmp = static_cast<double>(v_lon_o_(j, i)) + 360.0;
    } else {
      lon_tmp = static_cast<double>(v_lon_o_(j, i));
    }
    for (int ip = 1; ip < mx; ++ip) {
      if ((v_alon_(0, ip-1) <= lon_tmp) &&
          (v_alon_(0, ip)   >= lon_tmp)) {
        ic = ip;
        break;
      }
    }
    // ---find adjacent two grids on y direction (jc)
    double lat_tmp = v_lat_o_(j, i);
    for (int jp = 1; jp < my; ++jp) {
      if ((v_alat_(jp,   0) <= lat_tmp) &&
          (v_alat_(jp-1, 0) >= lat_tmp)) {
        jc = jp;
        break;
      }
    }
    if ((ic == INT_MAX) || (jc == INT_MAX)) {
      printf("break adjacent grids has no found\n");
      // TODO exit
			// exit(EXIT_FAILURE);
      return ;
    }
    // Bilinear interpolater
    v_buffer_(j*IMT_GLOBAL + i) = v_a_(jc, ic);

    return ;
  }
 private:
  const float SPVAL_     = CppOutputMod::SPVAL;
  ViewDouble2D v_a_      = *KokkosForcMod::p_v_s_work;
  ViewDouble2D v_alon_   = *KokkosForcMod::p_v_s_wx;
  ViewDouble2D v_alat_   = *KokkosForcMod::p_v_s_wy;
  ViewDouble1D v_buffer_ = *KokkosForcMod::p_v_buffer;
  ViewFloat2D v_lon_o_   = *KokkosPconstMod::p_v_lon_o;
  ViewFloat2D v_lat_o_   = *KokkosPconstMod::p_v_lat_o;
};

class FunctorNear2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j) const {
    v_buffer_(j*IMT_GLOBAL + IMT_GLOBAL-2) = v_buffer_(j*IMT_GLOBAL    );
    v_buffer_(j*IMT_GLOBAL + IMT_GLOBAL-1) = v_buffer_(j*IMT_GLOBAL + 1);
    return ;
  }
 private:
  ViewDouble1D v_buffer_ = *KokkosForcMod::p_v_buffer;
};

class FunctorJRADaily1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    double uu(0.0), vv(0.0);
    if ((j >= JSM - 1) && (j < JEM) && (i >= 1) && (i < IMM)) {
      const double epsln = 1e-25;
      const double tmp = v_vit_(0, 0, j, i) / (
          v_viv_(0, 0, j  , i) + v_viv_(0, 0, j  , i+1)
        + v_viv_(0, 0, j-1, i) + v_viv_(0, 0, j-1, i+1) + epsln);

      uu = tmp * (v_u_(0, 0, j  , i) + v_u_(0, 0, j  , i+1)
                + v_u_(0, 0, j-1, i) + v_u_(0, 0, j-1, i+1));
      vv = tmp * (v_v_(0, 0, j  , i) + v_v_(0, 0, j  , i+1)
                + v_v_(0, 0, j-1, i) + v_v_(0, 0, j-1, i+1));
    }
    // relative speed to surface currents
    v_windx_(j, i) = (v_wspdu3_(j, i) - uu) * v_vit_(0, 0, j, i);
    v_windy_(j, i) = (v_wspdv3_(j, i) + vv) * v_vit_(0, 0, j, i);
    // 1.0 is from mom4
    v_wspd3_(j, i) = std::sqrt(v_windx_(j, i) * v_windx_(j, i)
                             + v_windy_(j, i) * v_windy_(j, i) + 1.0) * v_vit_(0, 0, j, i);
    // using a transient temperature, not daily mean
    const double tok = 273.15;
    v_model_sst_(j, i) = (v_at_(0, 0, 0, j, i) + tok) * v_vit_(0, 0, j, i);

    v_zz_(j, i) = 10.0;

    v_qs_(j, i) = (0.98 * 640380 * std::exp(-5107.4 / v_model_sst_(j, i))
        / 1.22) * v_vit_(0, 0, j, i);
    // temperature to potential temperature
    v_theta_(j, i) = v_tsa3_(j, i) * std::pow(100000.0 / v_psa3_(j, i), 0.286)
        * v_vit_(0, 0, j, i);
    v_runoff_(0, j, i) = v_runoff3_(j, i) * v_vit_(0, 0, j, i);
    v_seaice_(0, j, i) = v_seaice3_(j, i) * v_vit_(0, 0, j, i);
    return ;
  }
 private:
  ViewDouble4D v_u_         = *KokkosDynMod::p_v_u;
  ViewDouble4D v_v_         = *KokkosDynMod::p_v_v;
  ViewDouble2D v_zz_        = *KokkosForcMod::p_v_zz;
  ViewDouble2D v_qs_        = *KokkosForcMod::p_v_qs;
  ViewDouble2D v_windx_     = *KokkosForcMod::p_v_windx;
  ViewDouble2D v_windy_     = *KokkosForcMod::p_v_windy;
  ViewDouble2D v_tsa3_      = *KokkosForcMod::p_v_tsa3;
  ViewDouble2D v_psa3_      = *KokkosForcMod::p_v_psa3;
  ViewDouble2D v_theta_     = *KokkosForcMod::p_v_theta;
  ViewDouble2D v_wspd3_     = *KokkosForcMod::p_v_wspd3;
  ViewDouble2D v_wspdu3_    = *KokkosForcMod::p_v_wspdu3;
  ViewDouble2D v_wspdv3_    = *KokkosForcMod::p_v_wspdv3;
  ViewDouble2D v_runoff3_   = *KokkosForcMod::p_v_runoff3;
  ViewDouble2D v_seaice3_   = *KokkosForcMod::p_v_seaice3;
  ViewDouble2D v_model_sst_ = *KokkosForcMod::p_v_model_sst;
  ViewDouble3D v_runoff_    = *KokkosForcMod::p_v_runoff;
  ViewDouble3D v_seaice_    = *KokkosForcMod::p_v_seaice;
  ViewDouble4D v_viv_       = *KokkosPconstMod::p_v_viv;
  ViewDouble4D v_vit_       = *KokkosPconstMod::p_v_vit;
  ViewDouble5D v_at_        = *KokkosTracerMod::p_v_at;
};

class FunctorJRADaily2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const double vit_times_one_minus_seaice = 
        v_vit_(0, 0, j, i) * (1.0 - v_seaice_(0, j, i));

    v_sshf_(0, j, i) = v_core_sensible_(j, i) * vit_times_one_minus_seaice;
    v_lthf_(0, j, i) = (v_core_latent_(j, i) - v_snow3_(j, i) * 3.335e5) 
        * vit_times_one_minus_seaice;

    double tmp_sst = v_model_sst_(j, i) * v_model_sst_(j, i);
    tmp_sst *= tmp_sst;
    v_lwv_(0, j, i) = (0.95*v_lwv3_(j, i) - 0.95*5.67e-8*tmp_sst)
        * vit_times_one_minus_seaice;

    v_qs_(j, i) =   v_core_tau_(j, i) * v_windx_(j, i) * v_vit_(0, 0, j, i);
    v_zz_(j, i) = - v_core_tau_(j, i) * v_windy_(j, i) * v_vit_(0, 0, j, i);

    v_nswv_(0, j, i) = v_lwv_(0, j, i) + v_sshf_(0, j, i) + v_lthf_(0, j, i);
    v_swv_(0, j, i) = 0.934 * v_swv3_(j, i) * vit_times_one_minus_seaice;

    v_fresh_(0, j, i) = - (v_core_latent_(j, i) / (2.5e+6) 
        + v_rain3_(j, i) + v_snow3_(j, i) + v_runoff3_(j, i))
            * vit_times_one_minus_seaice;

  v_ustar_(0, j, i) = v_ustar_(0, j, i) * v_vit_(0, 0, j, i);
  return ;
  }
 private:
  ViewDouble2D v_qs_  = *KokkosForcMod::p_v_qs;
  ViewDouble2D v_zz_  = *KokkosForcMod::p_v_zz;

  ViewDouble2D v_lwv3_    = *KokkosForcMod::p_v_lwv3;
  ViewDouble2D v_swv3_    = *KokkosForcMod::p_v_swv3;
  ViewDouble2D v_rain3_   = *KokkosForcMod::p_v_rain3;
  ViewDouble2D v_snow3_   = *KokkosForcMod::p_v_snow3;
  ViewDouble2D v_runoff3_ = *KokkosForcMod::p_v_runoff3;

  ViewDouble2D v_windx_ = *KokkosForcMod::p_v_windx;
  ViewDouble2D v_windy_ = *KokkosForcMod::p_v_windy;

  ViewDouble2D v_core_tau_      = *KokkosForcMod::p_v_core_tau;
  ViewDouble2D v_core_latent_   = *KokkosForcMod::p_v_core_latent;
  ViewDouble2D v_core_sensible_ = *KokkosForcMod::p_v_core_sensible;

  ViewDouble2D v_model_sst_ = *KokkosForcMod::p_v_model_sst;

  ViewDouble3D v_lwv_  = *KokkosForcMod::p_v_lwv;
  ViewDouble3D v_swv_  = *KokkosForcMod::p_v_swv;
  ViewDouble3D v_nswv_ = *KokkosForcMod::p_v_nswv;

  ViewDouble3D v_sshf_ = *KokkosForcMod::p_v_sshf;
  ViewDouble3D v_lthf_ = *KokkosForcMod::p_v_lthf;

  ViewDouble3D v_ustar_  = *KokkosForcMod::p_v_ustar;
  ViewDouble3D v_fresh_  = *KokkosForcMod::p_v_fresh;
  ViewDouble3D v_seaice_ = *KokkosForcMod::p_v_seaice;

  ViewDouble4D v_vit_ = *KokkosPconstMod::p_v_vit;
};

class FunctorJRADaily3 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    v_su_(0, j, i) = 0.25 * (
      v_qs_(j  , i) + v_qs_(j  , i-1) 
    + v_qs_(j+1, i) + v_qs_(j+1, i-1)) * v_viv_(0, 0, j, i);
    v_sv_(0, j, i) = 0.25 * (
      v_zz_(j  , i) + v_zz_(j  , i-1) 
    + v_zz_(j+1, i) + v_zz_(j+1, i-1)) * v_viv_(0, 0, j, i);
    // using CppParamMod::JST;
    // if (i >= 1 && i < JEM && j >= (JST-1) && j < JEM) {
    //   v_su_(0, j, i) = 0.25 * (
    //     v_qs_(j  , i) + v_qs_(j  , i-1) 
    //   + v_qs_(j+1, i) + v_qs_(j+1, i-1)) * v_viv_(0, 0, j, i);
    //   v_sv_(0, j, i) = 0.25 * (
    //     v_zz_(j  , i) + v_zz_(j  , i-1) 
    //   + v_zz_(j+1, i) + v_zz_(j+1, i-1)) * v_viv_(0, 0, j, i);
    // } else {
    //   v_su_(0, j, i) = 0.0;
    //   v_sv_(0, j, i) = 0.0;
    // }
  return ;
  }
 private:
  ViewDouble2D v_qs_  = *KokkosForcMod::p_v_qs;
  ViewDouble2D v_zz_  = *KokkosForcMod::p_v_zz;
  ViewDouble3D v_su_  = *KokkosForcMod::p_v_su;
  ViewDouble3D v_sv_  = *KokkosForcMod::p_v_sv;
  ViewDouble4D v_viv_ = *KokkosPconstMod::p_v_viv;
};

class FunctorNcarOceanFluxesJra {
 public:
  FunctorNcarOceanFluxesJra (
      const ViewDouble2D &v_u_del,
      const ViewDouble2D &v_t,
      const ViewDouble2D &v_ts,
      const ViewDouble2D &v_q,
      const ViewDouble2D &v_qs,
      const ViewDouble2D &v_z,
      const ViewDouble4D &v_avail,
      const ViewDouble2D &v_sh,
      const ViewDouble2D &v_lh,
      const ViewDouble2D &v_tau,
      const ViewDouble3D &v_ustar) :
    v_u_del_(v_u_del), v_t_(v_t), v_ts_(v_ts), v_q_(v_q), 
        v_qs_(v_qs), v_z_(v_z), v_avail_(v_avail), v_sh_(v_sh), 
            v_lh_(v_lh), v_tau_(v_tau), v_ustar_(v_ustar) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    double cd = 0.0;
    double ch = 0.0;
    double ce = 0.0;
            
    if (v_avail_(0, 0, j, i) > 0.5) {
      const double grav(9.80), vonkarm(0.40);
      const double reciprocal_vonkarm = 2.5;

      const double tv = v_t_(j, i) * (1.0 + 0.608 * v_q_(j, i));
      const double u  = std::max(v_u_del_(j, i), 0.5);
      double u10      = u;

      double cd_n10    = (2.7/u10 + 0.142 + 0.0764*u10) * 0.001;
      double cd_n10_rt = std::sqrt(cd_n10);
      double ce_n10    = 34.6 * cd_n10_rt * 0.001;
      double stab      = 0.5 + sign(0.5, v_t_(j, i) - v_ts_(j, i));
      double ch_n10    = (18.0 * stab + 32.7 * (1.0 - stab)) * cd_n10_rt * 0.001;
 
      cd = cd_n10;
      ch = ch_n10;
      ce = ce_n10;

      double bstar;
      for (int jj = 0; jj < 2; ++jj) {
        // for: 1, n_itts. n_itts = 2
      // Monin-Obukhov iteration
      const double cd_rt = std::sqrt(cd);
      v_ustar_(0, j, i) = cd_rt * u;
      const double tstar = (ch / cd_rt) * (v_t_(j, i) - v_ts_(j, i));
      const double qstar = (ce / cd_rt) * (v_q_(j, i) - v_qs_(j, i));
      bstar = grav * (tstar/tv + qstar/(v_q_(j, i) + 1.0/0.608));
 
      double zeta = vonkarm * bstar * v_z_(j, i) / 
          (v_ustar_(0, j, i) * v_ustar_(0, j, i));
      zeta = sign(std::min(std::abs(zeta), 10.0), zeta);
 
      double x2 = std::sqrt(std::abs(1 - 16.0 * zeta));
      x2 = std::max(x2, 1.0);
      const double x = std::sqrt(x2);
      double psi_m;
      double psi_h;
      if (zeta > 0.0) {
        psi_m = - 5.0 * zeta;
        psi_h = psi_m;
      } else {
        psi_m = std::log((1.0 + 2.0*x + x2) * (1.0+x2)*0.125) - 2.0 * (std::atan(x) - std::atan(1.0));
        psi_h = 2.0 * std::log((1.0 + x2) * 0.5);
      }
      u10 = u / (1.0 + cd_n10_rt * (std::log(v_z_(j, i) * 0.1) - psi_m) * reciprocal_vonkarm);
 
      cd_n10    = (2.7/u10 + 0.142 + 0.0764*u10) * 0.001;
      cd_n10_rt = std::sqrt(cd_n10);
      ce_n10    = 34.6 * cd_n10_rt * 0.001;
      stab      = 0.5 + sign(0.5, zeta);
      ch_n10    = (18.0*stab + 32.7*(1.0-stab)) * cd_n10_rt * 0.001;
      // diagnostic
      // z0     = 10 * exp(-vonkarm / cd_n10_rt);
      double xx = (std::log(v_z_(j, i) * 0.1) - psi_m) * reciprocal_vonkarm;
      const double tmp = (1.0 + cd_n10_rt * xx);
      cd        = cd_n10 / (tmp * tmp);
      xx = (std::log(v_z_(j, i) * 0.1) - psi_h) * reciprocal_vonkarm;
 
      ch = ch_n10 / (1.0 + ch_n10 * xx / cd_n10_rt) * std::sqrt(cd / cd_n10);
      ce = ce_n10 / (1.0 + ce_n10 * xx / cd_n10_rt) * std::sqrt(cd / cd_n10);
      }
    }
    const double L  = 2.5e6;
    const double R0 = 1.22;
    const double CP = 1000.5;
    v_sh_(j, i)  = R0 * CP * ch * 
        (v_t_(j, i) - v_ts_(j, i)) * v_u_del_(j, i);
    v_lh_(j, i)  = R0 * ce * L * 
      (v_q_(j, i) - v_qs_(j, i)) * v_u_del_(j, i);
    v_tau_(j, i) = R0 * cd * v_u_del_(j, i);
    return ;
  }
 private:
  const ViewDouble2D v_u_del_;
  const ViewDouble2D v_t_;
  const ViewDouble2D v_ts_;
  const ViewDouble2D v_q_;
  const ViewDouble2D v_qs_;
  const ViewDouble2D v_z_;
  const ViewDouble4D v_avail_;
  const ViewDouble2D v_sh_;
  const ViewDouble2D v_lh_;
  const ViewDouble2D v_tau_;
  const ViewDouble3D v_ustar_;

  template<typename T>
  KOKKOS_INLINE_FUNCTION T sign (const T &x, const T &y) const {
    return y >= static_cast<T>(0) ? std::abs(x) : -std::abs(x);
  }
};
KOKKOS_REGISTER_FOR_2D(FunctorInterplationNearest, FunctorInterplationNearest)
KOKKOS_REGISTER_FOR_2D(FunctorNear1, FunctorNear1)
KOKKOS_REGISTER_FOR_1D(FunctorNear2, FunctorNear2)
KOKKOS_REGISTER_FOR_2D(FunctorJRADaily1, FunctorJRADaily1)
KOKKOS_REGISTER_FOR_2D(FunctorJRADaily2, FunctorJRADaily2)
KOKKOS_REGISTER_FOR_2D(FunctorJRADaily3, FunctorJRADaily3)
KOKKOS_REGISTER_FOR_2D(FunctorNcarOceanFluxesJra, FunctorNcarOceanFluxesJra)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_JRA_DAILY_HPP_
