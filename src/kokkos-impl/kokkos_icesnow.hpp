#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_ICESNOW_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_ICESNOW_HPP_
#include "../head/def-undef.h"

#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_shr_const_mod.h"

#include "../head/kokkos_grid.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_tracer_mod.h"
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

#include <cmath>

using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;
using CppParamMod::MAX_BLOCKS_CLINIC;
//==================


class FunctorIcesnow1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &j, const int &i) const {
    const int k(0), iblock(0);
    const double sal_ocn(34.7), sal_ice(4.0);

    if (v_at_(iblock, 0, k, j, i) < tbice_ 
        && v_kmt_(iblock, j, i) > 0) {
  
      const double tdiff = tbice_ - v_at_(iblock, 0, k, j, i);
      v_licomqice_(iblock, j, i) += tdiff * v_dzp_(k) / v_dzp_(0);
  
      v_at_(iblock, 1, k, j, i) += tdiff * (sal_ocn - sal_ice) 
          * cp_ / heat_ice_fusion_ * 0.001;
  
      v_at_(iblock, 0, k, j, i) = tbice_;       
    } //end if

    if (v_licomqice_(iblock, j, i) > 0.0 && 
        v_at_(iblock, 0, k, j, i) > tbice_) {

      const double tdiff = fmin((v_at_(iblock, 0, 0, j, i) - tbice_), 
          v_licomqice_(iblock, j, i));

      v_licomqice_(iblock, j, i) -= tdiff;

      v_at_(iblock, 0, k, j, i)  -= tdiff;
      v_at_(iblock, 1, k, j, i)  -= tdiff * (sal_ocn - sal_ice) 
          * cp_ / heat_ice_fusion_ * 0.001;
    }//end if
    return ;
  }
 private:
  const double cp_              = CppPconstMod::cp;
  const double tbice_           = CppPconstMod::tbice;
  const double heat_ice_fusion_ = CppShrConstMod::SHR_CONST_LATICE;

  const ViewInt3D    v_kmt_       = *KokkosGrid     ::p_v_kmt;
  const ViewDouble1D v_dzp_       = *KokkosPconstMod::p_v_dzp;
  const ViewDouble3D v_licomqice_ = *KokkosTracerMod::p_v_licomqice;
  const ViewDouble5D v_at_        = *KokkosTracerMod::p_v_at;
};
class FunctorIcesnow2 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (
      const int &k, const int &j, const int &i) const {
    const int iblock = 0;
    if (v_at_(iblock, 0, k, j, i) < tbice_) {
        v_at_(iblock, 0, k, j, i) = tbice_ ;
    } //end if
    return ;
  }
 private:
  const double tbice_      = CppPconstMod::tbice;
  const ViewDouble5D v_at_ = *KokkosTracerMod::p_v_at;
};

KOKKOS_REGISTER_FOR_2D(FunctorIcesnow1, FunctorIcesnow1)
KOKKOS_REGISTER_FOR_3D(FunctorIcesnow2, FunctorIcesnow2)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_ICESNOW_HPP_
