#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_OUTPUT_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_OUTPUT_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

namespace CppOutputMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

#ifdef DAILYACC
extern float (&z0daily)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern float (&mlddaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&lthfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&sshfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&lwvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&swvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&precdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&evapdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&roffdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&sudaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&svdaily)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&tsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&ssdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&usdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&vsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&wsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
extern float (&z0mon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&himon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&hdmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&qicemon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&lthfmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&sshfmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&lwvmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&swvmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&sumon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&svmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&runoffmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&freshmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&wsmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&tsmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&ssmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&usmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&vsmon)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float (&icmon)[MAX_BLOCKS_CLINIC][2][JMT][IMT];

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
extern float (&athkdfmon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISO_TYPE_BF

#ifdef ISO
extern float (&axmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&aymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&azmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&dxmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&dymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&dzmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float (&wsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&vsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#ifdef ISOOUT
extern float (&azmon_vetisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&azmon_vntisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float (&azmon_vbtisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
extern float (&azmon_a3mon)[KM][JMT][IMT];
#endif // SMAG_OUT

extern double &err_norm2mon;
#endif // LOWRES

constexpr float SPVAL = 1.0e+35;
} // namespace CppOutputMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_OUTPUT_MOD_H_
