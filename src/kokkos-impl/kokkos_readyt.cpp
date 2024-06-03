#include "../head/def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_readyt.hpp"
extern "C" void readyt_debug_(
    double*wk1,double*wk2,double*wk3,double*wk4,double*wk5,double*wk6,double*wk7,double*wk8,double*wk9,double*wk10,double*wk11,double*wk12,double*wk13,double*wk14
    );
// extern "C" void readyt_debug_(
//     );

#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_tracer_mod.h"

void kokkos_readyt() {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

/*
  static auto dev = Kokkos::DefaultExecutionSpace();
  Kokkos::deep_copy (dev, *p_v_akt,      0.0);
  Kokkos::deep_copy (dev, *p_v_rit,      0.0);
  Kokkos::deep_copy (dev, *p_v_ric,      0.0);
  Kokkos::deep_copy (dev, *p_v_rict,     0.0);
  Kokkos::deep_copy (dev, *p_v_ricdt,    0.0);
  Kokkos::deep_copy (dev, *p_v_ricdttms, 0.0);
*/
using Kokkos::deep_copy;
using Kokkos::create_mirror_view;

//   ViewDouble4D::HostMirror h_v_utl = create_mirror_view(*p_v_utl);
//   ViewDouble4D::HostMirror h_v_vtl = create_mirror_view(*p_v_vtl);
//   ViewDouble5D::HostMirror h_v_akt = create_mirror_view(*p_v_akt);

//   ViewDouble4D::HostMirror h_v_dlu = create_mirror_view(*p_v_dlu);
//   ViewDouble4D::HostMirror h_v_dlv = create_mirror_view(*p_v_dlv);
//   ViewDouble4D::HostMirror h_v_pdensity = create_mirror_view(*p_v_pdensity);
//   ViewDouble4D::HostMirror h_v_rict = create_mirror_view(*p_v_rict);
//   ViewDouble4D::HostMirror h_v_gg = create_mirror_view(*p_v_gg);
//   ViewDouble4D::HostMirror h_v_pp = create_mirror_view(*p_v_pp);
//   ViewDouble4D::HostMirror h_v_alpha = create_mirror_view(*p_v_alpha);
//   ViewDouble4D::HostMirror h_v_beta = create_mirror_view(*p_v_beta);
//   ViewDouble4D::HostMirror h_v_ricdt = create_mirror_view(*p_v_ricdt);

//   UnManagedViewDouble3D h_v_buoytur (&CppForcMod::buoytur[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_buoysol (&CppForcMod::buoysol[0][0][0],
//       1, JMT, IMT); 

//   UnManagedViewDouble3D h_v_h0 (&CppDynMod::h0[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_h0f (&CppDynMod::h0f[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_h0l (&CppDynMod::h0l[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_psa (&CppForcMod::psa[0][0][0],
//       1, JMT, IMT); 
//   UnManagedViewDouble4D h_v_utf (&CppDynMod::utf[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_vtf (&CppDynMod::vtf[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_u (&CppDynMod::u[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble4D h_v_v (&CppDynMod::v[0][0][0][0],
//       1, KM, JMT, IMT); 
//   UnManagedViewDouble5D h_v_at (&CppTracerMod::at[0][0][0][0][0],
//       1, NTRA, KM, JMT, IMT); 
//   UnManagedViewDouble5D h_v_atb (&CppTracerMod::atb[0][0][0][0][0],
//       1, NTRA, KM+1, JMT, IMT); 
//   UnManagedViewDouble3D h_v_pxb (&CppWorkMod::pxb[0][0][0],
//       1, JMT, IMT);
//   UnManagedViewDouble3D h_v_pyb (&CppWorkMod::pyb[0][0][0],
//       1, JMT, IMT);

//   UnManagedViewDouble3D h_v_whx (&CppWorkMod::whx[0][0][0], 
//       1, JMT, IMT);
//   UnManagedViewDouble3D h_v_why (&CppWorkMod::why[0][0][0], 
//       1, JMT, IMT);

//   UnManagedViewDouble3D h_v_wgp (&CppWorkMod::wgp[0][0][0], 
//       1, JMT, IMT);

//   UnManagedViewDouble4D h_v_wka (&CppWorkMod::wka[0][0][0][0], 
//       1, KM, JMT, IMT);

//   UnManagedViewDouble3D h_v_work (&CppWorkMod::work[0][0][0], 
//       1, JMT, IMT);
//   UnManagedViewDouble3D h_v_pay (&CppWorkMod::pay[0][0][0], 
//       1, JMT, IMT);
//   UnManagedViewDouble3D h_v_pax (&CppWorkMod::pax[0][0][0], 
//       1, JMT, IMT);

//   deep_copy(*p_v_utf, h_v_utf);
//   deep_copy(*p_v_vtf, h_v_vtf);
//   deep_copy(*p_v_u, h_v_u);
//   deep_copy(*p_v_v, h_v_v);
//   deep_copy(*p_v_at, h_v_at);
//   deep_copy(*p_v_atb, h_v_atb);
//   deep_copy(*p_v_psa, h_v_psa);
//   deep_copy(*p_v_h0f, h_v_h0f);
//   deep_copy(*p_v_h0, h_v_h0);

  parallel_for ("readyt_1", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyt1());

  parallel_for ("readyt_2", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyt2());

  parallel_for ("readyt_3", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyt3());

  parallel_for ("readyt_4", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyt4());

  parallel_for ("readyt_5", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyt5());

  parallel_for ("readyt_6", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KMM1, JMT, IMT}, tile3D), FunctorReadyt6());

  parallel_for ("readyt_7", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorReadyt7());

  parallel_for ("readyt_8", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyt8());

  parallel_for ("readyt_9", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyt9());

  parallel_for ("readyt_10", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyt10());

  {
    parallel_for ("readyt_11", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorReadyt11());
  }

  // {
  //   int team_size   = 128;
  //   int league_size = (IMT * JMT + team_size - 1) / team_size;
  //   parallel_for ("readyt_tgrid_to_ugrid", Kokkos::TeamPolicy<>(
  //       league_size, team_size), FunctorReadyt111());
  // }
  // {
  //   const int tile_x = 16;
  //   const int tile_y = 8;
  //   const int calc_x = 15;
  //   const int calc_y = 7;
  //   const int team_size   = tile_x * tile_y;
  //   const int league_size = ((IMT + calc_x - 1) / calc_x) * ((JMT + calc_y - 1) / calc_y);
 
  //   using ScratchView = Kokkos::View<double*, Kokkos::DefaultExecutionSpace::scratch_memory_space>;
  //   size_t shmem_size = ScratchView::shmem_size (team_size);
 
  //   parallel_for ("readyt_tgrid_to_ugrid", Kokkos::TeamPolicy<>(
  //       league_size, team_size).set_scratch_size(0, Kokkos::PerTeam(shmem_size)), FunctorReadyt112());
  // }

  parallel_for ("readyt_12", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}), FunctorReadyt12());

  parallel_for ("readyt_13", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}), FunctorReadyt13());

//   deep_copy(h_v_utl, *p_v_utl);
//   deep_copy(h_v_vtl, *p_v_vtl);
//   deep_copy(h_v_utf, *p_v_utf);
//   deep_copy(h_v_vtf, *p_v_vtf);
//   deep_copy(h_v_akt, *p_v_akt);
//   deep_copy(h_v_rict, *p_v_rict);
//   deep_copy(h_v_pdensity, *p_v_pdensity);
//   deep_copy(h_v_gg, *p_v_gg);
//   deep_copy(h_v_pp, *p_v_pp);
//   deep_copy(h_v_alpha, *p_v_alpha);
//   deep_copy(h_v_beta, *p_v_beta);
//   deep_copy(h_v_buoytur, *p_v_buoytur);
//   deep_copy(h_v_buoysol, *p_v_buoysol);
//   deep_copy(h_v_ricdt, *p_v_ricdt);
//   deep_copy(h_v_dlu, *p_v_dlu);
//   deep_copy(h_v_dlv, *p_v_dlv);
//   deep_copy(h_v_h0l, *p_v_h0l);
//   deep_copy(h_v_h0f, *p_v_h0f);
//   deep_copy(h_v_pxb, *p_v_pxb);
//   deep_copy(h_v_pyb, *p_v_pyb);
//   deep_copy(h_v_wgp, *p_v_wgp);
//   deep_copy(h_v_work, *p_v_work);
//   deep_copy(h_v_whx, *p_v_whx);
//   deep_copy(h_v_why, *p_v_why);
//   deep_copy(h_v_pay, *p_v_pay);
//   deep_copy(h_v_pax, *p_v_pax);

//   readyt_debug_(
//     h_v_utl.data(),
//     h_v_vtl.data(),
//     h_v_utf.data(),
//     h_v_vtf.data(),
//     h_v_akt.data(),
//     h_v_rict.data(),
//     h_v_pdensity.data(),
//     h_v_gg.data(),
//     h_v_pp.data(),
//     h_v_alpha.data(),
//     h_v_beta.data(),
//     h_v_ricdt.data(),
//     h_v_dlu.data(),
//     h_v_dlv.data()
//     );
  return ;
}
#endif // LICOM_ENABLE_KOKKOS
