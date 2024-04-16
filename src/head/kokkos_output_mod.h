#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_OUTPUT_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_OUTPUT_MOD_H_

#include "def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

namespace KokkosOutputMod {

#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
extern ViewFloat4D *p_v_icmon;

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
#endif // ISO_TYPE_BF

#ifdef ISO
#ifdef ISOOUT
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
#endif // SMAG_OUT

#endif // LOWRES

} // namespace KokkosOutputMod

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_OUTPUT_MOD_H_
