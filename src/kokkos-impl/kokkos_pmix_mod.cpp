#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosPmixMod {

ViewDouble4D *p_v_ric          = nullptr;

ViewDouble4D *p_v_rict         = nullptr;
ViewDouble4D *p_v_rict_replace = nullptr;

ViewDouble3D *p_v_rict_ref     = nullptr;
ViewDouble4D *p_v_rit          = nullptr;
ViewDouble4D *p_v_riu          = nullptr;

ViewDouble4D *p_v_ricdt    = nullptr;
ViewDouble4D *p_v_ricdttms = nullptr;

ViewDouble4D *p_v_ridt     = nullptr;

ViewDouble4D *p_v_s2u = nullptr;
ViewDouble4D *p_v_s2t = nullptr;

#ifdef SOLAR
ViewDouble1D *p_v_pen = nullptr;
#endif // SOLAR

#ifdef SOLARCHLORO
ViewDouble4D *p_v_pen_chl = nullptr;
#endif // SOLARCHLORO
} // namespace KokkosPmixMod
#endif // LICOM_ENABLE_KOKKOS
