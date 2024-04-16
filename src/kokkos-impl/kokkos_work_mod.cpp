#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosWorkMod {

ViewDouble3D *p_v_pxb = nullptr;
ViewDouble3D *p_v_pyb = nullptr;

ViewDouble3D *p_v_pax = nullptr;
ViewDouble3D *p_v_pay = nullptr;

ViewDouble3D *p_v_whx = nullptr;
ViewDouble3D *p_v_why = nullptr;

ViewDouble3D *p_v_wgp = nullptr;

ViewDouble4D *p_v_wka = nullptr;

ViewDouble4D *p_v_work_1 = nullptr;
ViewDouble4D *p_v_work_2 = nullptr;
ViewDouble4D *p_v_work_3 = nullptr;

ViewDouble4D *p_v_temp = nullptr;

ViewDouble4D *p_v_uk = nullptr;
ViewDouble4D *p_v_vk = nullptr;

ViewDouble3D *p_v_work = nullptr;

ViewDouble1D *p_v_wkk = nullptr;

ViewDouble4D *p_v_wkb = nullptr;
ViewDouble4D *p_v_wkc = nullptr;
ViewDouble4D *p_v_wkd = nullptr;

ViewDouble4D *p_v_tf  = nullptr;
ViewDouble3D *p_v_stf = nullptr;

ViewFloat2D *p_v_buffer_real4 = nullptr;

ViewDouble2D *p_v_work1_g = nullptr;
ViewDouble2D *p_v_work2_g = nullptr;
ViewDouble2D *p_v_work3_g = nullptr;
} // namespace KokkosWorkMod
#endif // LICOM_ENABLE_KOKKOS
