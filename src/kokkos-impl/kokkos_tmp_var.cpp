#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"
namespace KokkosTmpVar {

// READYT
ViewDouble3D *p_v_work1 = nullptr;
ViewDouble3D *p_v_work2 = nullptr;

ViewDouble4D *p_v_pp  = nullptr;
ViewDouble4D *p_v_ppa = nullptr;
ViewDouble4D *p_v_ppb = nullptr;
ViewDouble4D *p_v_ppc = nullptr;

ViewDouble4D *p_v_alpha = nullptr;
ViewDouble4D *p_v_beta  = nullptr;

// READYC
ViewDouble4D *p_v_wp12 = nullptr;
ViewDouble4D *p_v_wp13 = nullptr;

ViewDouble4D *p_v_riv1 = nullptr;
ViewDouble4D *p_v_riv2 = nullptr;

ViewDouble3D *p_v_wk1 = nullptr;
ViewDouble3D *p_v_wk2 = nullptr;
ViewDouble3D *p_v_wk3 = nullptr;
ViewDouble3D *p_v_wp3 = nullptr;

#ifdef BCKMEX
ViewDouble3D *p_v_diff_back = nullptr;
ViewDouble3D *p_v_diff_back_sh = nullptr;
ViewDouble3D *p_v_diff_back_nn = nullptr;
#endif // BCKMEX
ViewDouble4D *p_v_uv_ws_face = nullptr;

#ifndef SMAG
ViewDouble2D *p_v_hduk = nullptr;
ViewDouble2D *p_v_hdvk = nullptr;

ViewDouble2D *p_v_div_out = nullptr;

#ifdef BIHAR
ViewDouble2D *p_v_curl = nullptr;

ViewDouble2D *p_v_cc = nullptr;
ViewDouble3D *p_v_d2uk = nullptr;
ViewDouble3D *p_v_d2vk = nullptr;
#endif // BIHAR
#endif // SMAG

// TRACER
ViewDouble4D *p_v_vtl_ori = nullptr;

ViewDouble3D *p_v_adv_tt = nullptr;

ViewDouble4D *p_v_at_00_max_min = nullptr;

#ifdef BIHAR
ViewDouble3D *p_v_dt2k = nullptr;
#endif // BIHAR

ViewInt1D    *p_v_nn = nullptr;
ViewDouble1D *p_v_xs = nullptr;

ViewDouble4D *p_v_c_cnsew = nullptr;
#ifdef ISO
#ifdef LDD97
ViewDouble4D *p_v_f1 = nullptr;
ViewDouble4D *p_v_f2 = nullptr;
#endif // LDD97
#endif // ISO

// BAROTR
ViewDouble2D *p_v_gradx = nullptr;
ViewDouble2D *p_v_grady = nullptr;

// POP Halo Update
ViewDouble1D *p_v_halo_buffer = nullptr;
} // namespace KokkosTmpVar

#endif // LICOM_ENABLE_KOKKOS
