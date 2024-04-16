#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosGrid {

ViewDouble1D *p_v_area_t_k        = nullptr;
                                                         
ViewDouble1D *p_v_volume_t_k      = nullptr;
ViewDouble1D *p_v_volume_t_marg_k = nullptr;

// ViewInt2D *p_v_kmt_g   = nullptr;
// ViewInt2D *p_v_basin_g = nullptr;

ViewDouble3D *p_v_dxu     = nullptr;
ViewDouble3D *p_v_dyu     = nullptr;
ViewDouble3D *p_v_dxt     = nullptr;
ViewDouble3D *p_v_dyt     = nullptr;
ViewDouble3D *p_v_dxur    = nullptr;
ViewDouble3D *p_v_dyur    = nullptr;
ViewDouble4D *p_v_dxyur   = nullptr;
ViewDouble3D *p_v_dxtr    = nullptr;
ViewDouble3D *p_v_dytr    = nullptr;
ViewDouble3D *p_v_hts     = nullptr;
ViewDouble3D *p_v_htw     = nullptr;
ViewDouble3D *p_v_hun     = nullptr;
ViewDouble3D *p_v_hue     = nullptr;
ViewDouble5D *p_v_h_tu_swen = nullptr;
ViewDouble3D *p_v_ulat    = nullptr;
ViewDouble3D *p_v_ulon    = nullptr;
ViewDouble3D *p_v_tlat    = nullptr;
ViewDouble3D *p_v_tlon    = nullptr;
ViewDouble3D *p_v_angle   = nullptr;
ViewDouble3D *p_v_anglet  = nullptr;
ViewDouble3D *p_v_fcor    = nullptr;
ViewDouble3D *p_v_fcort   = nullptr;
ViewDouble3D *p_v_uarea   = nullptr;
ViewDouble3D *p_v_tarea   = nullptr;
ViewDouble3D *p_v_uarea_r = nullptr;
ViewDouble3D *p_v_tarea_r = nullptr;
ViewDouble3D *p_v_ht      = nullptr;
ViewDouble3D *p_v_hu      = nullptr;
ViewDouble3D *p_v_hur     = nullptr;

ViewInt3D *p_v_basin  = nullptr;
                                                         
ViewInt3D *p_v_kmt    = nullptr;
ViewInt3D *p_v_kmu    = nullptr;
ViewInt3D *p_v_kmtold = nullptr;
                                                         
ViewDouble3D *p_v_rcalct = nullptr;
ViewDouble3D *p_v_rcalcu = nullptr;
                                                         
ViewInt3D *p_v_kmtn = nullptr;
ViewInt3D *p_v_kmts = nullptr;
ViewInt3D *p_v_kmte = nullptr;
ViewInt3D *p_v_kmtw = nullptr;
ViewInt3D *p_v_kmun = nullptr;
ViewInt3D *p_v_kmus = nullptr;
ViewInt3D *p_v_kmue = nullptr;
ViewInt3D *p_v_kmuw = nullptr;

// kmtn kmts kmte kmtw
// kmun kmus kmue kmuw
ViewInt4D *p_v_kmt_nsew = nullptr;
                                                         
ViewInt3D *p_v_kmtee = nullptr;
ViewInt3D *p_v_kmtnn = nullptr;
                                                         
// ViewInt3D *p_v_region_mask = nullptr;

// ViewDouble2D *p_v_tlat_g = nullptr;
// ViewDouble2D *p_v_tlon_g = nullptr;
                                                         
ViewDouble3D *p_v_at0    = nullptr;
ViewDouble3D *p_v_atn    = nullptr;
ViewDouble3D *p_v_ate    = nullptr;
ViewDouble3D *p_v_atne   = nullptr;
ViewDouble3D *p_v_au0    = nullptr;
ViewDouble3D *p_v_aus    = nullptr;
ViewDouble3D *p_v_auw    = nullptr;
ViewDouble3D *p_v_ausw   = nullptr;

} // namespace KokkosGrid
#endif // LICOM_ENABLE_KOKKOS
