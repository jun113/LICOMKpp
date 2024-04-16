#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL4_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL4_H_
#include "cpp_param_mod.h"
namespace CppHmixDel4 {

using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;

extern double* dt_nsew;
extern double* du_cnsewm;
extern double* dm_cnsew;

extern double (&dtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&duc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&due)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&duw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dmc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dmn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dms)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dme)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dmw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dum)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ahf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&amf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ratio_dxy)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern double &ah;
extern double &am;

} // CppHmixDel4
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL4_H_
