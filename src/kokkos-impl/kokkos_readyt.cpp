#include "../head/def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_readyt.hpp"

extern "C" void readyt_debug_(
    double* wk1,double* wk2,double* wk3,double* wk4,double* wk5,
    double* wk6,
    double* wk7
);

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

  readyt_debug_((*p_v_rict).data(),
  (*p_v_gg).data(),
  (*p_v_pp).data(),
  (*p_v_alpha).data(),
  (*p_v_beta).data(),
  (*p_v_dlu).data(),
  (*p_v_dlv).data()
  );
  return ;
}
#endif // LICOM_ENABLE_KOKKOS
