#include "../head/def-undef.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_tracer_mod.h"

#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_tracer_mod.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"
#endif // LICOM_ENABLE_KOKKOS

void daily_update_h2d() {

#ifdef LICOM_ENABLE_KOKKOS
  using CppParamMod::MAX_BLOCKS_CLINIC;
  using CppParamMod::KM;
  using CppParamMod::JMT;
  using CppParamMod::IMT;
  using CppParamMod::NX_BLOCK;
  using CppParamMod::NY_BLOCK;
  using CppParamMod::NTRA;

  using KokkosDynMod::p_v_u;
  using KokkosDynMod::p_v_v;
  using KokkosForcMod::p_v_su;
  using KokkosForcMod::p_v_sv;
  using KokkosForcMod::p_v_swv;
  using KokkosForcMod::p_v_sss;
  using KokkosForcMod::p_v_nswv;
  using KokkosTracerMod::p_v_at;
#ifdef FRC_CORE
  using KokkosForcMod::p_v_fresh;
  using KokkosForcMod::p_v_seaice;
#endif // FRC_CORE

  auto dev = Kokkos::DefaultExecutionSpace();

/*
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_su(&(CppForcMod::su[0][0][0]), 
              IMT, JMT, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_sv(&(CppForcMod::sv[0][0][0]), 
              IMT, JMT, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_swv(&(CppForcMod::swv[0][0][0]), 
              IMT, JMT, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_nswv(&(CppForcMod::nswv[0][0][0]), 
              IMT, JMT, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_u(&(CppDynMod::u[0][0][0][0]), 
              IMT, JMT, KM, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_v(&(CppDynMod::v[0][0][0][0]), 
              IMT, JMT, KM, MAX_BLOCKS_CLINIC); 

  static Kokkos::View<double *****, Layout,
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_at(&(CppTracerMod::at[0][0][0][0][0]), 
              IMT, JMT, KM, NTRA, MAX_BLOCKS_CLINIC); 

#ifdef FRC_CORE
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_fresh(&(CppForcMod::fresh[0][0][0]), 
              NX_BLOCK, NY_BLOCK, MAX_BLOCKS_CLINIC); 
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_seaice(&(CppForcMod::seaice[0][0][0]), 
              NX_BLOCK, NY_BLOCK, MAX_BLOCKS_CLINIC); 
#endif // FRC_CORE
  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_sss(&(CppForcMod::sss[0][0][0]), 
              NX_BLOCK, NY_BLOCK, MAX_BLOCKS_CLINIC); 
*/
  static ViewDouble3D::HostMirror h_v_su     = create_mirror_view(*p_v_su);
  static ViewDouble3D::HostMirror h_v_sv     = create_mirror_view(*p_v_sv);
  static ViewDouble3D::HostMirror h_v_sss    = create_mirror_view(*p_v_sss);
  static ViewDouble3D::HostMirror h_v_swv    = create_mirror_view(*p_v_swv);
  static ViewDouble3D::HostMirror h_v_nswv   = create_mirror_view(*p_v_nswv);
#ifdef FRC_CORE
  static ViewDouble3D::HostMirror h_v_fresh  = create_mirror_view(*p_v_fresh);
  static ViewDouble3D::HostMirror h_v_seaice = create_mirror_view(*p_v_seaice);
#endif // FRC_CORE
  static ViewDouble4D::HostMirror h_v_u      = create_mirror_view(*p_v_u);
  static ViewDouble4D::HostMirror h_v_v      = create_mirror_view(*p_v_v);

  static ViewDouble5D::HostMirror h_v_at     = create_mirror_view(*p_v_at);

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_su     (iblock, j, i) = CppForcMod::su[iblock][j][i];
        h_v_sv     (iblock, j, i) = CppForcMod::sv[iblock][j][i];
        h_v_sss    (iblock, j, i) = CppForcMod::sss[iblock][j][i];
        h_v_swv    (iblock, j, i) = CppForcMod::swv[iblock][j][i];
        h_v_nswv   (iblock, j, i) = CppForcMod::nswv[iblock][j][i];
#ifdef FRC_CORE
        h_v_fresh  (iblock, j, i) = CppForcMod::fresh[iblock][j][i];
        h_v_seaice (iblock, j, i) = CppForcMod::seaice[iblock][j][i];
#endif // FRC_CORE
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          h_v_u(iblock, k, j, i) = CppDynMod::u[iblock][k][j][i];
          h_v_v(iblock, k, j, i) = CppDynMod::v[iblock][k][j][i];
        }
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int n = 0; n < NTRA; ++n) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            h_v_at(iblock, n, k, j, i) = CppTracerMod::at[iblock][n][k][j][i];
          }
        }
      }
    }
  }

  Kokkos::deep_copy(dev, *p_v_u,    h_v_u);
  Kokkos::deep_copy(dev, *p_v_v,    h_v_v);
  Kokkos::deep_copy(dev, *p_v_at,   h_v_at);
  Kokkos::deep_copy(dev, *p_v_su,   h_v_su);
  Kokkos::deep_copy(dev, *p_v_sv,   h_v_sv);
  Kokkos::deep_copy(dev, *p_v_swv,  h_v_swv);
  Kokkos::deep_copy(dev, *p_v_sss,  h_v_sss);
  Kokkos::deep_copy(dev, *p_v_nswv, h_v_nswv);
#ifdef FRC_CORE
  Kokkos::deep_copy(dev, *p_v_fresh,  h_v_fresh);
  Kokkos::deep_copy(dev, *p_v_seaice, h_v_seaice);
#endif // FRC_CORE
#endif // LICOM_ENABLE_KOKKOS

  return ;
}

// copy back to do addps, nextstep
// at, h0, u, v, ws
void daily_update_d2h() {

#ifdef LICOM_ENABLE_KOKKOS
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  using CppParamMod::MAX_BLOCKS_CLINIC;
  using CppParamMod::KM;
  using CppParamMod::KMP1;
  using CppParamMod::JMT;
  using CppParamMod::IMT;
  using CppParamMod::NX_BLOCK;
  using CppParamMod::NY_BLOCK;
  using CppParamMod::NTRA;

  using KokkosDynMod   ::p_v_u;
  using KokkosDynMod   ::p_v_v;
  using KokkosDynMod   ::p_v_h0;
  using KokkosDynMod   ::p_v_ws;
  using KokkosForcMod  ::p_v_su;
  using KokkosForcMod  ::p_v_sv;
  using KokkosForcMod  ::p_v_swv;
  using KokkosForcMod  ::p_v_sshf;
  using KokkosForcMod  ::p_v_lthf;
  using KokkosForcMod  ::p_v_fresh;
  using KokkosTracerMod::p_v_at;
  using KokkosTracerMod::p_v_atb;

  auto dev = Kokkos::DefaultExecutionSpace();
  // h0 u v at[0] at[1] ws su sv swv sshf lthf fresh

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_h0(&(CppDynMod::h0[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_h0, *p_v_h0);

  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_u(&(CppDynMod::u[0][0][0][0]), 
              MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_u, *p_v_u);

  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_v(&(CppDynMod::v[0][0][0][0]), 
              MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_v, *p_v_v);

  /*
  static Kokkos::View<double *****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_at_sub(&(CppTracerMod::at[0][0][0][0][0]), 
              MAX_BLOCKS_CLINIC, 2, KM, JMT, IMT); 
  static auto v_at_sub = Kokkos::subview(*p_v_at,
      Kokkos::ALL, Kokkos::make_pair(0, 2), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  Kokkos::deep_copy(dev, h_v_at_sub, v_at_sub);
  */

  static Kokkos::View<double *****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_at(&(CppTracerMod::at[0][0][0][0][0]), 
              MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_at, *p_v_at);

  /*
  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_ws_sub(&(CppDynMod::ws[0][0][0][0]), 
              MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  static auto v_ws_sub = Kokkos::subview(*p_v_ws,
      Kokkos::ALL, Kokkos::make_pair(0, KM), Kokkos::ALL, Kokkos::ALL);
  Kokkos::deep_copy(dev, h_v_ws_sub, v_ws_sub);
  */
  static Kokkos::View<double ****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_ws(&(CppDynMod::ws[0][0][0][0]), 
              MAX_BLOCKS_CLINIC, KMP1, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_ws, *p_v_ws);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_su(&(CppForcMod::su[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_su, *p_v_su);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_sv(&(CppForcMod::sv[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_sv, *p_v_sv);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_swv(&(CppForcMod::swv[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_swv, *p_v_swv);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_sshf(&(CppForcMod::sshf[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_sshf, *p_v_sshf);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_lthf(&(CppForcMod::lthf[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_lthf, *p_v_lthf);

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_fresh(&(CppForcMod::fresh[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_fresh, *p_v_fresh);

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
#endif // LICOM_ENABLE_KOKKOS

  return ;
}