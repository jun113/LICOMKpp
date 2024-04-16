#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_TRACER_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_TRACER_MOD_H_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosTracerMod {

extern ViewDouble5D *p_v_atb;
extern ViewDouble4D *p_v_net;
extern ViewDouble5D *p_v_at;

//extern double (*(&restore_at))[NTRA][JMT][IMT];

extern ViewDouble4D *p_v_pdensity;

extern ViewDouble3D *p_v_amld;

extern ViewDouble5D *p_v_tend;

extern ViewDouble5D *p_v_ax;
extern ViewDouble5D *p_v_ay;
extern ViewDouble5D *p_v_az;

extern ViewDouble5D *p_v_dx;
extern ViewDouble5D *p_v_dy;
extern ViewDouble5D *p_v_dz;

extern ViewDouble4D *p_v_penetrate;

extern ViewDouble5D *p_v_dt_diff;

extern ViewDouble5D *p_v_ddy;

extern ViewDouble5D *p_v_dt_conv;
extern ViewDouble5D *p_v_dt_restore;

#ifdef ISO
extern ViewDouble5D *p_v_aay_iso;
extern ViewDouble5D *p_v_ddy_iso;

extern ViewDouble5D *p_v_ax_iso;
extern ViewDouble5D *p_v_ay_iso;
extern ViewDouble5D *p_v_az_iso;

extern ViewDouble5D *p_v_dx_iso;
extern ViewDouble5D *p_v_dy_iso;
extern ViewDouble5D *p_v_dz_iso;
#endif // ISO

extern ViewDouble3D *p_v_licomqice;
} // namespace KokkosTracerMod 

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_TRACER_MOD_H_
