#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CANUTO_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CANUTO_MOD_H_

#include "def-undef.h"
#if (defined(LICOM_ENABLE_KOKKOS) && defined(CANUTO))
#include "fortran_canuto_mod.h"
#include "kokkos_config.hpp"
namespace KokkosCanutoMod {
                                     
extern ViewDouble1D *p_v_and2on2a1;
extern ViewDouble1D *p_v_amtaun2a1;
                                     
extern ViewDouble1D *p_v_sma1;
extern ViewDouble1D *p_v_sha1;
extern ViewDouble1D *p_v_ssa1;
                                     
extern ViewDouble1D *p_v_rib;
extern ViewDouble1D *p_v_ridb;
                                     
extern ViewDouble2D *p_v_slq2b;

extern ViewDouble2D *p_v_smb;
extern ViewDouble2D *p_v_shb;
extern ViewDouble2D *p_v_ssb;
                                     
extern ViewInt1D *p_v_irimax;
                                     
extern ViewDouble1D *p_v_sisamax;
extern ViewDouble1D *p_v_ra_rmax;
extern ViewDouble1D *p_v_c_y_r0;
extern ViewDouble1D *p_v_back_ra_r;

extern ViewDouble1D *p_v_sm_r1;
extern ViewDouble1D *p_v_sh_r1;
extern ViewDouble1D *p_v_ss_r1;

extern ViewDouble1D *p_v_slq2_r1;

extern ViewDouble1D *p_v_ria;
extern ViewDouble1D *p_v_slq2a;

extern ViewDouble1D *p_v_sma;
extern ViewDouble1D *p_v_sha;

} // namespace KokkosCanutoMod
#endif // (defined(LICOM_ENABLE_KOKKOS) && defined(CANUTO))
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CANUTO_MOD_H_
