#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_DYN_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_DYN_MOD_H_

#include "cpp_param_mod.h"

namespace CppDynMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;

extern double (&ub)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&vb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ubp)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&vbp)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&h0p)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&up)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&vp)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&ws)[MAX_BLOCKS_CLINIC][KMP1][JMT][IMT];

extern double (&h0l)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&h0f)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&h0bl)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&h0bf)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&utl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&utf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&vtl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&vtf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&sbcx)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&bbcx)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&sbcy)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&bbcy)[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern double (*&buffer)[JMT_GLOBAL][IMT_GLOBAL];

extern double (&h0)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&u)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&v)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*&gg)[KM][JMT][IMT];
// extern double (*&dlu)[KM][JMT][IMT];
// extern double (*&dlv)[KM][JMT][IMT];

// extern double (*&dlub)[JMT][IMT];
// extern double (*&dlvb)[JMT][IMT];

extern double gg[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dlu[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dlv[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double dlub[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dlvb[MAX_BLOCKS_CLINIC][JMT][IMT];

} // namespace CppDynMod
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_DYN_MOD_H_
