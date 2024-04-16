#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_PMIX_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_PMIX_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

namespace CppPmixMod {

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMM1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

extern int &rtst;
extern int &rtend;
extern int &rust;
extern int &ruend;

extern double &wndmix;
extern double &fricmx;

extern double &diff_cbt_back;
extern double &diff_cbt_limit;

extern double &visc_cbu_back;
extern double &visc_cbu_limit;

// extern double (*&ric)[KMM1][JMT][IMT];
// extern double (*&rict)[KMM1][JMT][IMT];
// extern double (*&rict_replace)[KMM1][JMT][IMT];

extern double ric[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double rict[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double rict_replace[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

// extern double (*&rict_ref)[JMT][IMT];

// extern double (*&rit)[KMM1][JMT][IMT];
extern double rit[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
// extern double (*&riu)[KM+1][JMT][IMT];

extern double (&ricdt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double (&ricdttms)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double (&ridt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double (&s2u)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double (&s2t)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

#ifdef SOLAR
extern double (&pen)[KMM1];
#endif // SOLAR

#ifdef SOLARCHLORO
extern double (&pen_chl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // SOLARCHLORO
} // namespace CppPmixMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_PMIX_MOD_H_
