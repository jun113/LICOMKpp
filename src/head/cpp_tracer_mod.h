#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_TRACER_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_TRACER_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"
namespace CppTracerMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

extern double (&atb)[MAX_BLOCKS_CLINIC][NTRA][KM+1][JMT][IMT];
extern double (&net)[MAX_BLOCKS_CLINIC][NTRA][JMT][IMT];

extern double (&at)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
// extern double (*&restore_at)[NTRA][JMT][IMT];

extern double (&pdensity)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&amld)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&tend)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&ax)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&ay)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&az)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&dx)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&dy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&dz)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&penetrate)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&dt_diff)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&ddy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&dt_conv)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&dt_restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

#ifdef ISO
extern double (&aay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&ddy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&ax_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&ay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&az_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&dx_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&dy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&dz_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
#endif // ISO

extern double (&licomqice)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double &fw_norm2;

} // namespace CppTracerMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_TRACER_MOD_H_
