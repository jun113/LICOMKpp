#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS


#include "kokkos_convadj.hpp"
//---------------------------------------------
//    CONVADJ
//---------------------------------------------
void kokkos_convadj() {

  using CppParamMod::mytid;
  using CppPconstMod::dts;
  using CppPconstMod::ist;
  using CppPconstMod::adv_tracer;

  double c2dtts;
  std::string str_adv_tracer(adv_tracer);
  int flag_adv_tracer;
  if (str_adv_tracer.find("centered") != 
      str_adv_tracer.npos) {
    flag_adv_tracer = 0;
    if (ist >= 1) {
      c2dtts = dts * 2.0;
    } else {
      c2dtts = dts;
    }
  } else if (str_adv_tracer.find("tspas") !=
      str_adv_tracer.npos) {
    flag_adv_tracer = 1;
    c2dtts = dts;
  } else {
    if (mytid == 0) {
      printf("error in convadj\n");
      printf("The false advection option for tracer in convadj\n");
    }
    exit(0);
  }

#ifdef LOWRES
  Kokkos::parallel_for("convadj_1", Kokkos::MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorConvadj1 ());
#endif // LOWRES

  // adv tracer == tspas
  if (flag_adv_tracer == 1) {
  Kokkos::parallel_for("convadj_2", Kokkos::MDRangePolicy<Kokkos::Rank<4>> (
      koArr4D{0, 0, 0, 0}, koArr4D{NTRA, KM, JMT, IMT}), FunctorConvadj2(c2dtts));
  } else {
  Kokkos::parallel_for("convadj_3", Kokkos::MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorConvadj3(c2dtts));
  }

  return ;
}
//--------------------------------------
// End icesnow
//--------------------------------------
#endif // LICOM_ENABLE_KOKKOS
