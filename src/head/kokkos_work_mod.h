#ifndef LICOM3_KOKKOS_SRC_KOKKOS_WORK_MOD_H__
#define LICOM3_KOKKOS_SRC_KOKKOS_WORK_MOD_H__
#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosWorkMod {

extern ViewDouble3D *p_v_pxb;
extern ViewDouble3D *p_v_pyb;

extern ViewDouble3D *p_v_pax;
extern ViewDouble3D *p_v_pay;

extern ViewDouble3D *p_v_whx;
extern ViewDouble3D *p_v_why;

extern ViewDouble3D *p_v_wgp;

extern ViewDouble4D *p_v_wka;

extern ViewDouble4D *p_v_work_1;
extern ViewDouble4D *p_v_work_2;
extern ViewDouble4D *p_v_work_3;

extern ViewDouble4D *p_v_temp;

extern ViewDouble4D *p_v_uk;
extern ViewDouble4D *p_v_vk;

extern ViewDouble3D *p_v_work;

extern ViewDouble1D *p_v_wkk;

extern ViewDouble4D *p_v_wkb;
extern ViewDouble4D *p_v_wkc;
extern ViewDouble4D *p_v_wkd;

extern ViewDouble4D *p_v_tf;
extern ViewDouble3D *p_v_stf;

extern ViewFloat2D *p_v_buffer_real4;

extern ViewDouble2D *p_v_work1_g;
extern ViewDouble2D *p_v_work2_g;
extern ViewDouble2D *p_v_work3_g;
} // namespace KokkosWorkMod

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_KOKKOS_WORK_MOD_H__
