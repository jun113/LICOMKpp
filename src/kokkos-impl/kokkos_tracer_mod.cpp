#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosTracerMod {

ViewDouble5D *p_v_atb = nullptr;
ViewDouble4D *p_v_net = nullptr;
ViewDouble5D *p_v_at  = nullptr;

//extern double (*(&restore_at))[NTRA][JMT][IMT];

ViewDouble4D *p_v_pdensity = nullptr;

ViewDouble3D *p_v_amld = nullptr;

ViewDouble5D *p_v_tend = nullptr;

ViewDouble5D *p_v_ax = nullptr;
ViewDouble5D *p_v_ay = nullptr;
ViewDouble5D *p_v_az = nullptr;

ViewDouble5D *p_v_dx = nullptr;
ViewDouble5D *p_v_dy = nullptr;
ViewDouble5D *p_v_dz = nullptr;

ViewDouble4D *p_v_penetrate = nullptr;

ViewDouble5D *p_v_dt_diff = nullptr;

ViewDouble5D *p_v_ddy = nullptr;

ViewDouble5D *p_v_dt_conv = nullptr;
ViewDouble5D *p_v_dt_restore = nullptr;

#ifdef ISO
ViewDouble5D *p_v_aay_iso = nullptr;
ViewDouble5D *p_v_ddy_iso = nullptr;

ViewDouble5D *p_v_ax_iso = nullptr;
ViewDouble5D *p_v_ay_iso = nullptr;
ViewDouble5D *p_v_az_iso = nullptr;

ViewDouble5D *p_v_dx_iso = nullptr;
ViewDouble5D *p_v_dy_iso = nullptr;
ViewDouble5D *p_v_dz_iso = nullptr;
#endif // ISO

ViewDouble3D *p_v_licomqice  = nullptr;
} // namespace KokkosTracerMod 
#endif // LICOM_ENABLE_KOKKOS
