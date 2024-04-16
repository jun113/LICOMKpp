#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_NEXTSTEP_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_NEXTSTEP_HPP_

#include "../head/def-undef.h"

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_param_mod.h"

#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_tracer_mod.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

using KokkosDynMod::p_v_u;
using KokkosDynMod::p_v_v;
using KokkosDynMod::p_v_h0;
using KokkosDynMod::p_v_h0p;
using KokkosDynMod::p_v_h0f;
using KokkosDynMod::p_v_h0bf;
using KokkosDynMod::p_v_up;
using KokkosDynMod::p_v_vp;
using KokkosDynMod::p_v_up;
using KokkosDynMod::p_v_ub;
using KokkosDynMod::p_v_vb;
using KokkosDynMod::p_v_ubp;
using KokkosDynMod::p_v_vbp;
using KokkosDynMod::p_v_utf;
using KokkosDynMod::p_v_vtf;
using KokkosTracerMod::p_v_at;
using KokkosTracerMod::p_v_atb;

using KokkosPconstMod::p_v_dzp;
using KokkosPconstMod::p_v_ohbu;
using KokkosPconstMod::p_v_viv;

class FunctorNextStep1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    v_up_(iblock,  k, j, i) = v_u_(iblock, k, j, i);
    v_utf_(iblock, k, j, i) = v_u_(iblock, k, j, i);
    v_vp_(iblock,  k, j, i) = v_v_(iblock, k, j, i);
    v_vtf_(iblock, k, j, i) = v_v_(iblock, k, j, i);

    // NOTE: ATB[0:KM]
    v_atb_(iblock, 0, k+1, j, i) = v_at_(iblock, 0, k, j, i);
    v_atb_(iblock, 1, k+1, j, i) = v_at_(iblock, 1, k, j, i);
    return ;
  };
 private:
  const ViewDouble4D v_u_   = *p_v_u;
  const ViewDouble4D v_v_   = *p_v_v;
  const ViewDouble4D v_up_  = *p_v_up;
  const ViewDouble4D v_vp_  = *p_v_vp;
  const ViewDouble4D v_utf_ = *p_v_utf;
  const ViewDouble4D v_vtf_ = *p_v_vtf;
  const ViewDouble5D v_at_  = *p_v_at;
  const ViewDouble5D v_atb_ = *p_v_atb;
};

class FunctorNextStep2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int iblock = 0;
    vinteg (j, i, v_u_, v_ub_);
    vinteg (j, i, v_v_, v_vb_);

    v_h0p_(iblock,  j, i) = v_h0_(iblock, j, i);
    v_h0f_(iblock,  j, i) = v_h0_(iblock, j, i);
    v_h0bf_(iblock, j, i) = v_h0_(iblock, j, i);
    v_ubp_(iblock,  j, i) = v_ub_(iblock, j, i);
    v_vbp_(iblock,  j, i) = v_vb_(iblock, j, i);
    return ;
  };
 private:
  const ViewDouble1D v_dzp_  = *p_v_dzp;
  const ViewDouble3D v_ub_   = *p_v_ub;
  const ViewDouble3D v_vb_   = *p_v_vb;
  const ViewDouble3D v_ubp_  = *p_v_ubp;
  const ViewDouble3D v_vbp_  = *p_v_vbp;
  const ViewDouble3D v_h0_   = *p_v_h0;
  const ViewDouble3D v_h0p_  = *p_v_h0p;
  const ViewDouble3D v_h0f_  = *p_v_h0f;
  const ViewDouble3D v_h0bf_ = *p_v_h0bf;
  const ViewDouble3D v_ohbu_ = *p_v_ohbu;
  const ViewDouble4D v_u_    = *p_v_u;
  const ViewDouble4D v_v_    = *p_v_v;
  const ViewDouble4D v_viv_  = *p_v_viv;

  KOKKOS_INLINE_FUNCTION void vinteg (const int &j, const int &i,
      const ViewDouble4D &v_wk3, const ViewDouble3D &v_wk2) 
          const {
    const int iblock = 0;
    v_wk2(iblock, j, i) = CppConstantMod::C0;
    for (int k = 0; k < CppParamMod::KM; ++k) {
      v_wk2(iblock, j, i) += v_dzp_(k) * v_ohbu_(iblock, j, i) 
          * v_wk3(iblock, k, j, i) *v_viv_(iblock, k, j, i);
    }
    return ;
  }
};

KOKKOS_REGISTER_FOR_3D(FunctorNextStep1,  FunctorNextStep1)
KOKKOS_REGISTER_FOR_2D(FunctorNextStep2,  FunctorNextStep2)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_NEXTSTEP_HPP_
