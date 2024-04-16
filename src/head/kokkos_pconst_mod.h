#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PCONST_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PCONST_MOD_H_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosPconstMod {

extern ViewInt1D *p_v_j_global;
extern ViewInt1D *p_v_i_global;

extern ViewDouble4D *p_v_vit;
extern ViewDouble4D *p_v_viv;

extern ViewDouble1D *p_v_ahv_back;

extern ViewInt3D *p_v_na;

#if (defined NETCDF) || (defined ALL)
extern ViewFloat1D *p_v_lon;
extern ViewFloat1D *p_v_lat;

extern ViewFloat2D *p_v_lon_o;
extern ViewFloat2D *p_v_lat_o;
extern ViewFloat2D *p_v_ulon_o;
extern ViewFloat2D *p_v_ulat_o;

extern ViewFloat1D *p_v_lev;
extern ViewFloat1D *p_v_lev1;
#endif // NETCDF || ALL

extern ViewDouble2D *p_v_s_lon;
extern ViewDouble2D *p_v_s_lat;

extern ViewDouble1D *p_v_zkt;
extern ViewDouble1D *p_v_dzp;
extern ViewDouble1D *p_v_odzp;
extern ViewDouble1D *p_v_odzt;

extern ViewDouble2D *p_v_odz_pt;

extern ViewDouble1D *p_v_zkp;

extern ViewDouble3D *p_v_ebea;
extern ViewDouble3D *p_v_ebeb;
extern ViewDouble3D *p_v_ebla;
extern ViewDouble3D *p_v_eblb;

extern ViewDouble3D *p_v_epea;
extern ViewDouble3D *p_v_epeb;
extern ViewDouble3D *p_v_epla;
extern ViewDouble3D *p_v_eplb;

extern ViewDouble3D *p_v_rrd1;
extern ViewDouble3D *p_v_rrd2;

extern ViewDouble3D *p_v_ohbt;
extern ViewDouble3D *p_v_ohbu;
extern ViewDouble3D *p_v_dzph;

extern ViewDouble3D *p_v_hbx;
extern ViewDouble3D *p_v_hby;

extern ViewDouble3D *p_v_snlat;

extern ViewDouble1D *p_v_to;
extern ViewDouble1D *p_v_so;
extern ViewDouble2D *p_v_c;
extern ViewDouble1D *p_v_po;

extern ViewDouble4D *p_v_akmu;
extern ViewDouble4D *p_v_akmt;

extern ViewDouble5D *p_v_akt;

extern ViewDouble1D *p_v_am;
extern ViewDouble1D *p_v_ah;

extern ViewDouble4D *p_v_am3;
extern ViewDouble4D *p_v_ah3;

extern ViewDouble4D *p_v_amx;
extern ViewDouble4D *p_v_amy;

extern ViewInt1D *p_v_nmonth ;
extern ViewInt1D *p_v_nnmonth;

#ifdef TIDEMIX
extern ViewDouble4D *p_v_fz_tide;
extern ViewDouble4D *p_v_ak_tide;

extern ViewDouble4D *p_v_fztidal;
extern ViewDouble4D *p_v_richardson;
extern ViewDouble4D *p_v_wp3_tidal;
extern ViewDouble4D *p_v_ak_tide1;
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
extern ViewDouble4D *p_v_wp1_canuto;
extern ViewDouble4D *p_v_wp2_canuto;
extern ViewDouble4D *p_v_wp3_canuto;
extern ViewDouble4D *p_v_wp4_canuto;
extern ViewDouble4D *p_v_wp5_canuto;
extern ViewDouble4D *p_v_wp6_canuto;
extern ViewDouble4D *p_v_wp7_canuto;
extern ViewDouble4D *p_v_wp8_canuto;
extern ViewDouble4D *p_v_wp10_canuto;
extern ViewDouble4D *p_v_wp11_canuto;
extern ViewDouble4D *p_v_wp12_canuto;
extern ViewDouble4D *p_v_wp13_canuto;

extern ViewDouble4D *p_v_wk1_canuto;
extern ViewDouble4D *p_v_wk2_canuto;
extern ViewDouble4D *p_v_wk3_canuto;
extern ViewDouble4D *p_v_wk4_canuto;

extern ViewDouble4D *p_v_fcor_canuto;
extern ViewDouble4D *p_v_fcort_canuto;
extern ViewDouble4D *p_v_alpha_canuto;
extern ViewDouble4D *p_v_beta_canuto;
#endif // CANUTOMIXOUT
} // namespace KokkosPconstMod

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PCONST_MOD_H_
