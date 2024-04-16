#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosHmixDel2 {

ViewDouble3D *p_v_dtn = nullptr;
ViewDouble3D *p_v_dts = nullptr;
ViewDouble3D *p_v_dte = nullptr;
ViewDouble3D *p_v_dtw = nullptr;
ViewDouble3D *p_v_duc = nullptr;
ViewDouble3D *p_v_dun = nullptr;
ViewDouble3D *p_v_dus = nullptr;
ViewDouble3D *p_v_due = nullptr;
ViewDouble3D *p_v_duw = nullptr;
ViewDouble3D *p_v_dmc = nullptr;
ViewDouble3D *p_v_dmn = nullptr;
ViewDouble3D *p_v_dms = nullptr;
ViewDouble3D *p_v_dme = nullptr;
ViewDouble3D *p_v_dmw = nullptr;
ViewDouble3D *p_v_dum = nullptr;
ViewDouble3D *p_v_ahf = nullptr;
ViewDouble3D *p_v_amf = nullptr;
} // namespace KokkosHmixDel2
#endif // LICOM_ENABLE_KOKKOS
