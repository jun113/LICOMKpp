#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS


#include "kokkos_icesnow.hpp"
//---------------------------------------------
//    ICESNOW
//---------------------------------------------
void kokkos_icesnow() {

  Kokkos::parallel_for("icesnow_1", Kokkos::MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorIcesnow1());

  Kokkos::parallel_for("icesnow_2", Kokkos::MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorIcesnow2());

  return ;
}
//--------------------------------------
// End icesnow
//--------------------------------------

#endif // LICOM_ENABLE_KOKKOS
