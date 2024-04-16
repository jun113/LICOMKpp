#ifndef LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_CONVADJ_HPP_
#define LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_CONVADJ_HPP_

#include "../head/def-undef.h"

#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

#include "../head/kokkos_grid.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_tracer_mod.h"
#include "../head/kokkos_output_mod.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

#include <cstdio>

using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;
using CppParamMod::MAX_BLOCKS_CLINIC;

class FunctorConvadj1 {
 public:
  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {
    const int iblock = 0;
    const int kcon = v_kmt_(iblock, j, i);
    if (kcon == 0){
      return ;
    }  
    int lctot =   0;
#ifdef LOWRES   
    int lcven =   1;
#endif // LOWRES   
    int lcon  = - 1;
    double rhoup[KM] = {0.0};
    double rholo[KM] = {0.0};
    for (int l = 0; l < KM - 1; ++l) {
      int l1 = l + 1;
      const double tup = v_at_(iblock, 0, l1, j, i) - v_to_(l1);
      const double sup = v_at_(iblock, 1, l1, j, i) - v_so_(l1);
      const double tlo = v_at_(iblock, 0, l, j, i) - v_to_(l1);
      const double slo = v_at_(iblock, 1, l, j, i) - v_so_(l1);
      rhoup[l1] = dens(tup, sup, l1);
      rholo[l]  = dens(tlo, slo, l1);
    }
    for (int k = kcon-1; k >= 1; --k) {
      if (rholo[k-1] > rhoup[k]) {
        lcon = k - 1;
      }
    }
    if (lcon < 0) {
      return ;
    }
    int lcona;
    int lconb;
    double dztsum;
    double tramix;
    double trasum[2];
    for (;;) { //start conv_1
      lcona = lcon;
      lconb = lcon + 1;
      dztsum = v_dzp_(lcona) + v_dzp_(lconb);
      for (int n = 0; n < 2; ++n) {
        trasum[n] = v_at_(iblock, n, lcona, j, i) * v_dzp_(lcona)
                  + v_at_(iblock, n, lconb, j, i) * v_dzp_(lconb);

        tramix = trasum[n] / dztsum;
        v_at_(iblock, n, lcona, j, i) = tramix;
        v_at_(iblock, n, lconb, j, i) = tramix;
      } //end for
      for (;;) { //start conv_2
        if (lconb != (kcon - 1)) {
          int l1 = lconb + 1;
          rholo[lconb] = dens(
              v_at_(iblock, 0, lconb, j, i) - v_to_(l1),
              v_at_(iblock, 1, lconb, j, i) - v_so_(l1), l1);

          if (rholo[lconb] > rhoup[l1]) {
            ++lconb;
            dztsum += v_dzp_(lconb);  
            for (int n = 0; n < 2; ++n) {
              trasum[n] += v_at_(iblock, n, lconb, j, i) 
                  * v_dzp_(lconb);
                   
              tramix = trasum[n] / dztsum;
              for (int lmix = lcona; lmix <= lconb; ++lmix) {
                v_at_(iblock, n, lmix, j, i) = tramix;
              }
            } //end for
            continue ;
          } //end if
        }   //end if
        if (lcona > 0) {
          int l1 = lcona - 1;
          rholo[l1] = dens(
              v_at_(iblock, 0, l1, j, i) - v_to_(lcona),
              v_at_(iblock, 1, l1, j, i) - v_so_(lcona),
                  lcona);
          rhoup[lcona] = dens(
              v_at_(iblock, 0, lcona, j, i) - v_to_(lcona), 
              v_at_(iblock, 1, lcona, j, i) - v_so_(lcona), 
                  lcona);

          if (rholo[lcona - 1] > rhoup[lcona]) {
            --lcona;
            dztsum += v_dzp_(lcona);
            for (int n = 0; n < 2; ++n) {
              trasum[n] += v_at_(iblock, n, lcona, j, i) * v_dzp_(lcona);
              tramix = trasum[n] / dztsum;
              for (int lmix = lcona; lmix <= lconb; ++lmix) {
                v_at_(iblock, n, lmix, j, i) = tramix;
              }
            }//end for
            continue ;
          }//end if
        }//end if
        break;
      }//end conv_2
      lctot += lconb - lcona + 1;
#ifdef LOWRES
      if (lcona == 0) {
        lcven = lconb - lcona + 1;
      }
#endif
      if (lconb == kcon - 1) {
#ifdef LOWRES
        v_icmon_(iblock, 0, j, i) += lctot;
        v_icmon_(iblock, 1, j, i) += lcven;
#endif
        return;
      }//end if
      lcon = lconb;
      while(true) {
        lcon =lcon + 1;
        if (lcon == kcon - 1) {
#ifdef LOWRES   
          v_icmon_(iblock, 0, j, i) += lctot;
          v_icmon_(iblock, 0, j, i) += lctot; 
          v_icmon_(iblock, 1, j, i) += lcven;
#endif
          return;
        }
        if (lcon == (KM-1)) {
          break ;
        }
        if (rholo[lcon] <= rhoup[lcon + 1]) {
          continue ;
        }
        break;
      }//end conv_3
    }//end conv_1
    return ;
  }
  KOKKOS_INLINE_FUNCTION double
      dens (const double &tq, const double &sq, const int &kk) const {
    double dens;
    dens = (v_c_(0, kk) + (v_c_(3, kk) + v_c_(6, kk) * sq) * sq +
           (v_c_(2, kk) +  v_c_(7, kk) * sq + v_c_(5, kk) * tq) * tq) * tq +
           (v_c_(1, kk) + (v_c_(4, kk) + v_c_(8, kk) * sq) * sq) * sq;
    return dens;
  }
 private:
  const ViewDouble1D v_so_   = *KokkosPconstMod::p_v_so;
  const ViewDouble1D v_to_   = *KokkosPconstMod::p_v_to;
  const ViewDouble1D v_dzp_  = *KokkosPconstMod::p_v_dzp;
  const ViewDouble2D v_c_    = *KokkosPconstMod::p_v_c;
  const ViewInt3D    v_kmt_  = *KokkosGrid::p_v_kmt;
#ifdef LOWRES               
  const ViewFloat4D v_icmon_ = *KokkosOutputMod::p_v_icmon;
#endif                      
  const ViewDouble5D v_at_   = *KokkosTracerMod::p_v_at;
};

class FunctorConvadj2 {
 public:
  FunctorConvadj2 (const double &c2dtts) : c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void operator () (
      const int &n, const int &k, const int &j, const int &i) const {
    const int iblock = 0;

    v_dt_conv_(iblock, n, k, j, i) = (v_at_(iblock, n, k, j, i) 
        - v_atb_(iblock, n, k+1, j, i)) 
            / c2dtts_ * v_vit_(iblock, k, j, i);
    
    v_tend_(iblock, n, k, j, i) += v_dt_conv_(iblock, n, k, j, i);
    v_atb_(iblock, n, k+1, j, i) = v_at_(iblock, n, k, j, i);

    return ;
  }
 private:
  const double c2dtts_;
  const ViewDouble4D v_vit_     = *KokkosPconstMod::p_v_vit;
  const ViewDouble5D v_at_      = *KokkosTracerMod::p_v_at;
  const ViewDouble5D v_atb_     = *KokkosTracerMod::p_v_atb;
  const ViewDouble5D v_tend_    = *KokkosTracerMod::p_v_tend;
  const ViewDouble5D v_dt_conv_ = *KokkosTracerMod::p_v_dt_conv;
};
class FunctorConvadj3 {
 public:
  FunctorConvadj3 (const double &c2dtts) : c2dtts_(c2dtts) {}

  KOKKOS_INLINE_FUNCTION void operator () (const int &j, const int &i) const {
    const int iblock = 0;
    for (int n = 0; n < NTRA; ++n) {
      for (int k = 0; k < KM; ++k) {
        v_dt_conv_(iblock, n, k, j, i) = (v_at_(iblock, n, k, j, i) 
            - v_atb_(iblock, n, k+1, j, i)) 
                / c2dtts_ * v_vit_(iblock, k, j, i);
     
        v_tend_(iblock, n, k, j, i) += v_dt_conv_(iblock, n, k, j, i);
      }
    }
    return ;
  }
 private:
  const double c2dtts_;
  const ViewDouble4D v_vit_     = *KokkosPconstMod::p_v_vit;
  const ViewDouble5D v_at_      = *KokkosTracerMod::p_v_at;
  const ViewDouble5D v_atb_     = *KokkosTracerMod::p_v_atb;
  const ViewDouble5D v_tend_    = *KokkosTracerMod::p_v_tend;
  const ViewDouble5D v_dt_conv_ = *KokkosTracerMod::p_v_dt_conv;
};


KOKKOS_REGISTER_FOR_2D(FunctorConvadj1, FunctorConvadj1)
KOKKOS_REGISTER_FOR_4D(FunctorConvadj2, FunctorConvadj2)
KOKKOS_REGISTER_FOR_2D(FunctorConvadj3, FunctorConvadj3)

#endif // LICOM3_KOKKKOS_SRC_KOKKOS_IMPL_KOKKOS_CONVADJ_HPP_
