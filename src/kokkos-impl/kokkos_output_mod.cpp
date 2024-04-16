#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosOutputMod {

#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES

ViewFloat4D *p_v_icmon = nullptr;

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
#endif // ISO_TYPE_BF

#ifdef ISO
#ifdef ISOOUT
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
#endif // SMAG_OUT

#endif // LOWRES


} // namespace KokkosOutputMod
#endif // LICOM_ENABLE_KOKKOS
