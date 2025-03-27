#include "../head/def-undef.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include "../head/licom_test_time.hpp"

#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_tracer_mod.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"
#endif // LICOM_ENABLE_KOKKOS

#include <iostream>

void daily_update_h2d(const bool &jra_from_device, 
                      const bool &next_step_from_device);

void daily_force_input (const bool &split_jra,
                        const bool &jra_from_device, 
                        const bool &next_step_from_device) {

  using MyTest::my_time;

  if (split_jra && jra_from_device) {
    if (CppParamMod::mytid == 0) {
      std::cout<< "err in daily force input, split_jra: " << split_jra <<
                  " jra_from_device: " << jra_from_device << std::endl;
    }
    fortran_mpi_barrier_();
    fortran_mpi_finalize_();
    exit (0);
  }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_start("jra daily");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.start_once();

  // JRA
  if (split_jra) {
    jra_daily_split_();
  } else {
#if (defined LICOM_ENABLE_FORTRAN) || defined(LICOM_ENABLE_CPP)
    // jra_daily_();
#if defined(LOWRES)
    cpp_jra_daily_low(CppPconstMod::iday);
#elif defined(HIGHRES) || defined(SUPHIGH)
    cpp_jra_daily_high(CppPconstMod::iday);
#endif // defined(HIGHRES) || defined(SUPHIGH)
#elif defined (LICOM_ENABLE_KOKKOS)
    if (jra_from_device) {
#if defined(LOWRES)
      kokkos_jra_daily_low(CppPconstMod::iday);
#elif defined(HIGHRES) || defined(SUPHIGH)
      kokkos_jra_daily_high(CppPconstMod::iday);
#endif // defined(HIGHRES) || defined(SUPHIGH)
    } else {
#if defined(LOWRES)
      cpp_jra_daily_low(CppPconstMod::iday);
#elif defined(HIGHRES) || defined(SUPHIGH)
      cpp_jra_daily_high(CppPconstMod::iday);
#endif // defined(HIGHRES) || defined(SUPHIGH)
    }  // jra from device
#endif // (defined LICOM_ENABLE_VERSION)
  }

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_stop("jra daily");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.end_once();

  if (CppParamMod::mytid == 0) {
    printf("jra_daily time: %.3f s\n", my_time.t_once);
  }
  // =======================
  // time interplate CHLOROPHYLL
  cpp_time_interplate_chlorophyll ();
  // =======================
#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_start("daily_h2d");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.start_daily_h2d();

  // daily memcpy
  daily_update_h2d (jra_from_device, next_step_from_device);

#ifdef LICOM_ENABLE_TEST_TIME
  my_time.testTime_stop("daily_h2d");
#endif // LICOM_ENABLE_TEST_TIME
  my_time.end_daily_h2d();

  return ;
}

void daily_update_h2d(const bool &jra_from_device, 
                      const bool &next_step_from_device) {

#ifdef LICOM_ENABLE_KOKKOS
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  using CppParamMod::MAX_BLOCKS_CLINIC;
  using CppParamMod::KM;
  using CppParamMod::JMT;
  using CppParamMod::IMT;
  using CppParamMod::NTRA;
  using CppParamMod::NX_BLOCK;
  using CppParamMod::NY_BLOCK;

  static auto dev = Kokkos::DefaultExecutionSpace();
  // time_interplate_chlorophyll
  static UnManagedViewDouble3D h_v_sss (&CppForcMod::sss[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy (dev, *KokkosForcMod::p_v_sss, h_v_sss);

  // jra daily
  if (!jra_from_device) {
    static UnManagedViewDouble3D h_v_su (&CppForcMod::su[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_su, h_v_su);
  
    static UnManagedViewDouble3D h_v_sv (&CppForcMod::sv[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_sv, h_v_sv);
  
    static UnManagedViewDouble3D h_v_fresh (&CppForcMod::fresh[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_fresh, h_v_fresh);
  
    static UnManagedViewDouble3D h_v_nswv (&CppForcMod::nswv[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_nswv, h_v_nswv);
  
    static UnManagedViewDouble3D h_v_swv (&CppForcMod::swv[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_swv, h_v_swv);
  
    static UnManagedViewDouble3D h_v_seaice (&CppForcMod::seaice[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosForcMod::p_v_seaice, h_v_seaice);
  }

  // addps
  static UnManagedViewDouble3D h_v_h0 (&CppDynMod::h0[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy (dev, *KokkosDynMod::p_v_h0, h_v_h0);

  // next step
  if (!next_step_from_device) {
    static UnManagedViewDouble4D h_v_up (&CppDynMod::up[0][0][0][0],
        MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_up, h_v_up);

    static UnManagedViewDouble4D h_v_vp (&CppDynMod::vp[0][0][0][0],
        MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_vp, h_v_vp);

    static UnManagedViewDouble4D h_v_utf (&CppDynMod::utf[0][0][0][0],
        MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_utf, h_v_utf);

    static UnManagedViewDouble4D h_v_vtf (&CppDynMod::vtf[0][0][0][0],
        MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_vtf, h_v_vtf);

    static UnManagedViewDouble5D h_v_atb (&(CppTracerMod::atb[0][0][0][0][0]), 
        MAX_BLOCKS_CLINIC, NTRA, KM+1, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosTracerMod::p_v_atb, h_v_atb);

    static UnManagedViewDouble3D h_v_ub (&CppDynMod::ub[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_ub, h_v_ub);

    static UnManagedViewDouble3D h_v_vb (&CppDynMod::vb[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_vb, h_v_vb);

    static UnManagedViewDouble3D h_v_h0p (&CppDynMod::h0p[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_h0p, h_v_h0p);

    static UnManagedViewDouble3D h_v_ubp (&CppDynMod::ubp[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_ubp, h_v_ubp);

    static UnManagedViewDouble3D h_v_vbp (&CppDynMod::vbp[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_vbp, h_v_vbp);

    static UnManagedViewDouble3D h_v_h0f (&CppDynMod::h0f[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_h0f, h_v_h0f);

    static UnManagedViewDouble3D h_v_h0bf (&CppDynMod::h0bf[0][0][0],
        MAX_BLOCKS_CLINIC, JMT, IMT); 
    Kokkos::deep_copy (dev, *KokkosDynMod::p_v_h0bf, h_v_h0bf);
  }

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
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

// copy back to energy
// u v h0 at
void energy_d2h() {

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
  using KokkosTracerMod::p_v_at;

  auto dev = Kokkos::DefaultExecutionSpace();
  // u v h0 at

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

  static Kokkos::View<double ***, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_h0(&(CppDynMod::h0[0][0][0]), 
              MAX_BLOCKS_CLINIC, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_h0, *p_v_h0);

  static Kokkos::View<double *****, Layout, 
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_at(&(CppTracerMod::at[0][0][0][0][0]), 
              MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT); 
  Kokkos::deep_copy(dev, h_v_at, *p_v_at);

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
#endif // LICOM_ENABLE_KOKKOS

  return ;
}