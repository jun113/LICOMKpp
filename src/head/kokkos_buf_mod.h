#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_BUF_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_BUF_MOD_H_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosBufMod {

#ifdef USE_OCN_CARBON

#endif // USE_OCN_CARBON

extern ViewDouble3D *p_v_ifrac;

} // namespace KokkosBufMod

#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_BUF_MOD_H_
