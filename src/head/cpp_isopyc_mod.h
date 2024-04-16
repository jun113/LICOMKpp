#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_ISOPYC_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_ISOPYC_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

namespace CppIsopycMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
#ifdef ISO
constexpr int NRPL = 5; 

extern double (&dptlim)[NRPL+1];

extern double (&fzisop)[KM];

extern double &slmxr;

extern double (&ahisop)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&athkdf)[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern double (*(&e))[3][JMT][KMP1][IMT];
// extern double (*(&rhoi))[NRPL][JMT][KM+1][IMT];

// extern double (*(&k1))[1][JMT][KM+1][IMT];
// extern double (*(&k2))[1][JMT][KM+1][IMT];
// extern double (*(&k3))[3][JMT][KM+1][IMT];

// extern double (*(&adv_vetiso))[JMT][KM][IMT];

extern double (&adv_vbtiso)[MAX_BLOCKS_CLINIC][JMT][KM+1][IMT];
extern double (&adv_vntiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT];

#ifdef isopycmixspatialvar
extern float (&dciso1)[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern float (&dciso2)[MAX_BLOCKS_CLINIC][JMT][KM][IMT];

extern float &dslope;
extern float &slopec;
#endif // isopycmixspatialvar

extern double (&kisrpl)[KM];

extern int (&krplin)[NRPL];

extern double (&zt)[KM];

extern double (&dzw)[KM+1];
extern double (&dzwr)[KM+1];

extern double (&dzr)[KM];

extern double (&tmask)[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern double (&f3)[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // ISO
} // namespace CppIsopycMod
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_ISOPYC_MOD_H_
