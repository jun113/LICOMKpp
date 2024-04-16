#include "../head/def-undef.h"
#if (defined(LICOM_ENABLE_KOKKOS) && defined(CANUTO))
#include "../head/fortran_canuto_mod.h"
#include "../head/kokkos_config.hpp"
namespace KokkosCanutoMod {
                                     
ViewDouble1D *p_v_and2on2a1 = nullptr;
ViewDouble1D *p_v_amtaun2a1 = nullptr;
                                     
ViewDouble1D *p_v_sma1 = nullptr;
ViewDouble1D *p_v_sha1 = nullptr;
ViewDouble1D *p_v_ssa1 = nullptr;
                                     
ViewDouble1D *p_v_rib  = nullptr;
ViewDouble1D *p_v_ridb = nullptr;
                                     
ViewDouble2D *p_v_slq2b = nullptr;

ViewDouble2D *p_v_smb = nullptr;
ViewDouble2D *p_v_shb = nullptr;
ViewDouble2D *p_v_ssb = nullptr;
                                     
ViewInt1D *p_v_irimax = nullptr;
                                     
ViewDouble1D *p_v_sisamax = nullptr;
ViewDouble1D *p_v_ra_rmax = nullptr;
ViewDouble1D *p_v_c_y_r0  = nullptr;

ViewDouble1D *p_v_back_ra_r = nullptr;

ViewDouble1D *p_v_sm_r1 = nullptr;
ViewDouble1D *p_v_sh_r1 = nullptr;
ViewDouble1D *p_v_ss_r1 = nullptr;

ViewDouble1D *p_v_slq2_r1 = nullptr;

ViewDouble1D *p_v_ria   = nullptr;
ViewDouble1D *p_v_slq2a = nullptr;

ViewDouble1D *p_v_sma = nullptr;
ViewDouble1D *p_v_sha = nullptr;

} // namespace KokkosCanutoMod
#endif // (defined(LICOM_ENABLE_KOKKOS) && defined(CANUTO))
