#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_ISOPYC_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_ISOPYC_MOD_H_

#include "def-undef.h"

#if (defined(ISO) && defined(LICOM_ENABLE_KOKKOS))
#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosIsopycMod {

extern ViewDouble1D *p_v_dptlim;

extern ViewDouble1D *p_v_fzisop;

extern ViewDouble3D *p_v_ahisop;

extern ViewDouble3D *p_v_athkdf;

extern ViewDouble5D *p_v_e;

extern ViewDouble5D *p_v_rhoi;

extern ViewDouble5D *p_v_k1;
extern ViewDouble5D *p_v_k2;
extern ViewDouble5D *p_v_k3;

extern ViewDouble4D *p_v_adv_vetiso;
extern ViewDouble4D *p_v_adv_vbtiso;
extern ViewDouble4D *p_v_adv_vntiso;

#ifdef isopycmixspatialvar
extern ViewFloat4D *p_v_dciso1;
extern ViewFloat4D *p_v_dciso2;
#endif // isopycmixspatialvar

extern ViewDouble1D *p_v_kisrpl;

extern ViewInt1D *p_v_krplin;

extern ViewDouble1D *p_v_zt;

extern ViewDouble1D *p_v_dzw;
extern ViewDouble1D *p_v_dzwr;

extern ViewDouble1D *p_v_dzr;

extern ViewDouble4D *p_v_tmask;
extern ViewDouble3D *p_v_f3;

} // namespace KokkosIsopycMod

#endif // (defined(ISO) && defined(LICOM_ENABLE_KOKKOS))
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_ISOPYC_MOD_H_
