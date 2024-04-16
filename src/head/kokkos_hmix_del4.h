#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_hmix_del4_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_hmix_del4_H_
#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS


#include "Kokkos_Core.hpp"
namespace KokkosHmixDel4 {

extern ViewDouble3D *p_v_dtn;
extern ViewDouble3D *p_v_dts;
extern ViewDouble3D *p_v_dte;
extern ViewDouble3D *p_v_dtw;
extern ViewDouble3D *p_v_duc;
extern ViewDouble3D *p_v_dun;
extern ViewDouble3D *p_v_dus;
extern ViewDouble3D *p_v_due;
extern ViewDouble3D *p_v_duw;
extern ViewDouble3D *p_v_dmc;
extern ViewDouble3D *p_v_dmn;
extern ViewDouble3D *p_v_dms;
extern ViewDouble3D *p_v_dme;
extern ViewDouble3D *p_v_dmw;
extern ViewDouble3D *p_v_dum;
extern ViewDouble3D *p_v_ahf;
extern ViewDouble3D *p_v_amf;
extern ViewDouble3D *p_v_ratio_dxy;

extern ViewDouble4D *p_v_dt_nsew;
extern ViewDouble4D *p_v_du_cnsewm;
extern ViewDouble4D *p_v_dm_cnsew;
} // namespace KokkosHmixDel4

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_hmix_del4_H_
