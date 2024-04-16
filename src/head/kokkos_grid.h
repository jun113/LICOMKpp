#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_GRID_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_GRID_H_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosGrid {

extern ViewDouble1D *p_v_area_t_k;
                                                         
extern ViewDouble1D *p_v_volume_t_k;
extern ViewDouble1D *p_v_volume_t_marg_k;

// extern ViewInt2D *p_v_kmt_g;
// extern ViewInt2D *p_v_basin_g;

extern ViewDouble3D *p_v_dxu;
extern ViewDouble3D *p_v_dyu;
extern ViewDouble3D *p_v_dxt;
extern ViewDouble3D *p_v_dyt;
extern ViewDouble3D *p_v_dxur;
extern ViewDouble3D *p_v_dyur;
extern ViewDouble4D *p_v_dxyur;
extern ViewDouble3D *p_v_dxtr;
extern ViewDouble3D *p_v_dytr;
extern ViewDouble3D *p_v_hts;
extern ViewDouble3D *p_v_htw;
extern ViewDouble3D *p_v_hun;
extern ViewDouble3D *p_v_hue;
extern ViewDouble5D *p_v_h_tu_swen;
extern ViewDouble3D *p_v_ulat;
extern ViewDouble3D *p_v_ulon;
extern ViewDouble3D *p_v_tlat;
extern ViewDouble3D *p_v_tlon;
extern ViewDouble3D *p_v_angle;
extern ViewDouble3D *p_v_anglet;
extern ViewDouble3D *p_v_fcor;
extern ViewDouble3D *p_v_fcort;
extern ViewDouble3D *p_v_uarea;
extern ViewDouble3D *p_v_tarea;
extern ViewDouble3D *p_v_uarea_r;
extern ViewDouble3D *p_v_tarea_r;
extern ViewDouble3D *p_v_ht;
extern ViewDouble3D *p_v_hu;
extern ViewDouble3D *p_v_hur;

extern ViewInt3D *p_v_basin;
                                                         
extern ViewInt3D *p_v_kmt;
extern ViewInt3D *p_v_kmu;
extern ViewInt3D *p_v_kmtold;
                                                         
extern ViewDouble3D *p_v_rcalct;
extern ViewDouble3D *p_v_rcalcu;
                                                         
extern ViewInt3D *p_v_kmtn;
extern ViewInt3D *p_v_kmts;
extern ViewInt3D *p_v_kmte;
extern ViewInt3D *p_v_kmtw;
extern ViewInt3D *p_v_kmun;
extern ViewInt3D *p_v_kmus;
extern ViewInt3D *p_v_kmue;
extern ViewInt3D *p_v_kmuw;

// kmtn kmts kmte kmtw
// kmun kmus kmue kmuw

extern ViewInt4D *p_v_kmt_nsew;

extern ViewInt3D *p_v_kmtee;
extern ViewInt3D *p_v_kmtnn;
                                                         
// extern ViewInt3D *p_v_region_mask;

// extern ViewDouble2D *p_v_tlat_g;
// extern ViewDouble2D *p_v_tlon_g;
                                                         
extern ViewDouble3D *p_v_at0;
extern ViewDouble3D *p_v_atn;
extern ViewDouble3D *p_v_ate;
extern ViewDouble3D *p_v_atne;
extern ViewDouble3D *p_v_au0;
extern ViewDouble3D *p_v_aus;
extern ViewDouble3D *p_v_auw;
extern ViewDouble3D *p_v_ausw;

} // namespace KokkosGrid

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_GRID_H_
