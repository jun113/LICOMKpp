#include "../head/def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_nextstep.hpp"

#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_pconst_mod.h"

void kokkos_nextstep() {
  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;
  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::KM;
  using CppParamMod::MAX_BLOCKS_CLINIC;

  // update h0 in addps.
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  auto dev = Kokkos::DefaultExecutionSpace();
  static Kokkos::View<double ***, Layout,
      Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          h_v_h0(&(CppDynMod::h0[0][0][0]),
              MAX_BLOCKS_CLINIC, JMT, IMT);
  Kokkos::deep_copy(dev, *p_v_h0, h_v_h0);
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  parallel_for ("next_step_1", MDRangePolicy<Kokkos::Rank<3>> (
      koArr3D{0, 0, 0}, koArr3D{KM, JMT, IMT}, tile3D), FunctorNextStep1());

  parallel_for ("next_step_2", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorNextStep2());

  CppPconstMod::isb = 0;
  CppPconstMod::isc = 0;
  CppPconstMod::ist = 0;

  return;
}
#endif // LICOM_ENABLE_KOKKOS
