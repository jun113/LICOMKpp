#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_DYN_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_DYN_MOD_H_
#include "def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosDynMod {

extern ViewDouble3D *p_v_ub;
extern ViewDouble3D *p_v_vb;

extern ViewDouble3D *p_v_ubp;
extern ViewDouble3D *p_v_vbp;

extern ViewDouble3D *p_v_h0p;

extern ViewDouble4D *p_v_up;
extern ViewDouble4D *p_v_vp;

extern ViewDouble4D *p_v_ws;

extern ViewDouble3D *p_v_h0l;
extern ViewDouble3D *p_v_h0f;
extern ViewDouble3D *p_v_h0bl;
extern ViewDouble3D *p_v_h0bf;

extern ViewDouble4D *p_v_utl;
extern ViewDouble4D *p_v_utf;
extern ViewDouble4D *p_v_vtl;
extern ViewDouble4D *p_v_vtf;

extern ViewDouble3D *p_v_sbcx;
extern ViewDouble3D *p_v_bbcx;
extern ViewDouble3D *p_v_sbcy;
extern ViewDouble3D *p_v_bbcy;

extern ViewDouble2D *p_v_buffer;

extern ViewDouble3D *p_v_h0;

extern ViewDouble4D *p_v_u;
extern ViewDouble4D *p_v_v;

extern ViewDouble4D *p_v_gg;

extern ViewDouble4D *p_v_dlu;
extern ViewDouble4D *p_v_dlv;

extern ViewDouble3D *p_v_dlub;
extern ViewDouble3D *p_v_dlvb;
} // namespace KokkosDynMod 

#endif // LICOM_ENABLE_KOKKOS

#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_DYN_MOD_H_
