#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosPconstMod {

ViewInt1D *p_v_j_global = nullptr;
ViewInt1D *p_v_i_global = nullptr;

ViewDouble4D *p_v_vit = nullptr;
ViewDouble4D *p_v_viv = nullptr;

ViewDouble1D *p_v_ahv_back = nullptr;

ViewInt3D *p_v_na = nullptr;

#if (defined NETCDF) || (defined ALL)
ViewFloat1D *p_v_lon = nullptr;
ViewFloat1D *p_v_lat = nullptr;

ViewFloat2D *p_v_lon_o  = nullptr;
ViewFloat2D *p_v_lat_o  = nullptr;
ViewFloat2D *p_v_ulon_o = nullptr;
ViewFloat2D *p_v_ulat_o = nullptr;

ViewFloat1D *p_v_lev  = nullptr;
ViewFloat1D *p_v_lev1 = nullptr;
#endif // NETCDF || ALL

ViewDouble2D *p_v_s_lon = nullptr;
ViewDouble2D *p_v_s_lat = nullptr;

ViewDouble1D *p_v_zkt  = nullptr;
ViewDouble1D *p_v_dzp  = nullptr;
ViewDouble1D *p_v_odzp = nullptr;
ViewDouble1D *p_v_odzt = nullptr;

ViewDouble2D *p_v_odz_pt = nullptr;

ViewDouble1D *p_v_zkp  = nullptr;

ViewDouble3D *p_v_ebea = nullptr;
ViewDouble3D *p_v_ebeb = nullptr;
ViewDouble3D *p_v_ebla = nullptr;
ViewDouble3D *p_v_eblb = nullptr;

ViewDouble3D *p_v_epea = nullptr;
ViewDouble3D *p_v_epeb = nullptr;
ViewDouble3D *p_v_epla = nullptr;
ViewDouble3D *p_v_eplb = nullptr;

ViewDouble3D *p_v_rrd1 = nullptr;
ViewDouble3D *p_v_rrd2 = nullptr;

ViewDouble3D *p_v_ohbt = nullptr;
ViewDouble3D *p_v_ohbu = nullptr;
ViewDouble3D *p_v_dzph = nullptr;

ViewDouble3D *p_v_hbx  = nullptr;
ViewDouble3D *p_v_hby  = nullptr;

ViewDouble3D *p_v_snlat = nullptr;

ViewDouble1D *p_v_to = nullptr;
ViewDouble1D *p_v_so = nullptr;
ViewDouble2D *p_v_c  = nullptr;
ViewDouble1D *p_v_po = nullptr;

ViewDouble4D *p_v_akmu = nullptr;
ViewDouble4D *p_v_akmt = nullptr;

ViewDouble5D *p_v_akt = nullptr;

ViewDouble1D *p_v_am = nullptr;
ViewDouble1D *p_v_ah = nullptr;

ViewDouble4D *p_v_am3 = nullptr;
ViewDouble4D *p_v_ah3 = nullptr;

ViewDouble4D *p_v_amx = nullptr;
ViewDouble4D *p_v_amy = nullptr;

ViewInt1D *p_v_nmonth  = nullptr;
ViewInt1D *p_v_nnmonth = nullptr;

#ifdef TIDEMIX
ViewDouble4D *p_v_fz_tide = nullptr;
ViewDouble4D *p_v_ak_tide = nullptr;

ViewDouble4D *p_v_fztidal    = nullptr;
ViewDouble4D *p_v_richardson = nullptr;
ViewDouble4D *p_v_wp3_tidal  = nullptr;
ViewDouble4D *p_v_ak_tide1   = nullptr;
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
ViewDouble4D *p_v_wp1_canuto  = nullptr;
ViewDouble4D *p_v_wp2_canuto  = nullptr;
ViewDouble4D *p_v_wp3_canuto  = nullptr;
ViewDouble4D *p_v_wp4_canuto  = nullptr;
ViewDouble4D *p_v_wp5_canuto  = nullptr;
ViewDouble4D *p_v_wp6_canuto  = nullptr;
ViewDouble4D *p_v_wp7_canuto  = nullptr;
ViewDouble4D *p_v_wp8_canuto  = nullptr;
ViewDouble4D *p_v_wp10_canuto = nullptr;
ViewDouble4D *p_v_wp11_canuto = nullptr;
ViewDouble4D *p_v_wp12_canuto = nullptr;
ViewDouble4D *p_v_wp13_canuto = nullptr;

ViewDouble4D *p_v_wk1_canuto  = nullptr;
ViewDouble4D *p_v_wk2_canuto  = nullptr;
ViewDouble4D *p_v_wk3_canuto  = nullptr;
ViewDouble4D *p_v_wk4_canuto  = nullptr;

ViewDouble4D *p_v_fcor_canuto  = nullptr;
ViewDouble4D *p_v_fcort_canuto = nullptr;
ViewDouble4D *p_v_alpha_canuto = nullptr;
ViewDouble4D *p_v_beta_canuto  = nullptr;
#endif // CANUTOMIXOUT
} // namespace KokkosPconstMod
#endif // LICOM_ENABLE_KOKKOS
