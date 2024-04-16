#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosDynMod {

ViewDouble3D *p_v_ub = nullptr;
ViewDouble3D *p_v_vb = nullptr;

ViewDouble3D *p_v_ubp = nullptr;
ViewDouble3D *p_v_vbp = nullptr;

ViewDouble3D *p_v_h0p = nullptr;

ViewDouble4D *p_v_up = nullptr;
ViewDouble4D *p_v_vp = nullptr;

ViewDouble4D *p_v_ws = nullptr;

ViewDouble3D *p_v_h0l = nullptr;
ViewDouble3D *p_v_h0f = nullptr;
ViewDouble3D *p_v_h0bl = nullptr;
ViewDouble3D *p_v_h0bf = nullptr;

ViewDouble4D *p_v_utl = nullptr;
ViewDouble4D *p_v_utf = nullptr;
ViewDouble4D *p_v_vtl = nullptr;
ViewDouble4D *p_v_vtf = nullptr;

ViewDouble3D *p_v_sbcx = nullptr;
ViewDouble3D *p_v_bbcx = nullptr;
ViewDouble3D *p_v_sbcy = nullptr;
ViewDouble3D *p_v_bbcy = nullptr;

ViewDouble2D *p_v_buffer = nullptr;

ViewDouble3D *p_v_h0 = nullptr;

ViewDouble4D *p_v_u = nullptr;
ViewDouble4D *p_v_v = nullptr;

ViewDouble4D *p_v_gg = nullptr;

ViewDouble4D *p_v_dlu = nullptr;
ViewDouble4D *p_v_dlv = nullptr;

ViewDouble3D *p_v_dlub = nullptr;
ViewDouble3D *p_v_dlvb = nullptr;

} // namespace KokkosDynMod 
#endif // LICOM_ENABLE_KOKKOS
