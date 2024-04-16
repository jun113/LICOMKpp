#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#ifdef ISO

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"
namespace KokkosIsopycMod {

ViewDouble1D *p_v_dptlim = nullptr;

ViewDouble1D *p_v_fzisop = nullptr;

ViewDouble3D *p_v_ahisop = nullptr;

ViewDouble3D *p_v_athkdf = nullptr;

ViewDouble5D *p_v_e = nullptr;

ViewDouble5D *p_v_rhoi = nullptr;

ViewDouble5D *p_v_k1 = nullptr;
ViewDouble5D *p_v_k2 = nullptr;
ViewDouble5D *p_v_k3 = nullptr;

ViewDouble4D *p_v_adv_vetiso = nullptr;
ViewDouble4D *p_v_adv_vbtiso = nullptr;
ViewDouble4D *p_v_adv_vntiso = nullptr;

#ifdef isopycmixspatialvar
ViewFloat4D *p_v_dciso1 = nullptr;
ViewFloat4D *p_v_dciso2 = nullptr;
#endif // isopycmixspatialvar

ViewDouble1D *p_v_kisrpl = nullptr;

ViewInt1D *p_v_krplin = nullptr;

ViewDouble1D *p_v_zt = nullptr;

ViewDouble1D *p_v_dzw = nullptr;
ViewDouble1D *p_v_dzwr = nullptr;

ViewDouble1D *p_v_dzr = nullptr;

ViewDouble4D *p_v_tmask = nullptr;
ViewDouble3D *p_v_f3 = nullptr;
} // namespace KokkosIsopycMod

#endif // ISO
#endif // LICOM_ENABLE_KOKKOS
