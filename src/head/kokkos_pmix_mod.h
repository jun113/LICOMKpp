#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PMIX_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PMIX_MOD_H_

#include "def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "Kokkos_Core.hpp"

namespace KokkosPmixMod {

extern ViewDouble4D *p_v_ric;

extern ViewDouble4D *p_v_rict;
extern ViewDouble4D *p_v_rict_replace;

extern ViewDouble3D *p_v_rict_ref;
extern ViewDouble4D *p_v_rit;
extern ViewDouble4D *p_v_riu;

extern ViewDouble4D *p_v_ricdt;
extern ViewDouble4D *p_v_ricdttms;

extern ViewDouble4D *p_v_ridt;

extern ViewDouble4D *p_v_s2u;
extern ViewDouble4D *p_v_s2t;

#ifdef SOLAR
extern ViewDouble1D *p_v_pen;
#endif // SOLAR

#ifdef SOLARCHLORO
extern ViewDouble4D *p_v_pen_chl;
#endif // SOLARCHLORO
} // namespace KokkosPmixMod

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_PMIX_MOD_H_
