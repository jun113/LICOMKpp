#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL2_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL2_H_
#include "cpp_param_mod.h"
namespace CppHmixDel2 {

using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

extern double (*&dtn)[NY_BLOCK][NX_BLOCK];
extern double (*&dts)[NY_BLOCK][NX_BLOCK];
extern double (*&dte)[NY_BLOCK][NX_BLOCK];
extern double (*&dtw)[NY_BLOCK][NX_BLOCK];
extern double (*&duc)[NY_BLOCK][NX_BLOCK];
extern double (*&dun)[NY_BLOCK][NX_BLOCK];
extern double (*&dus)[NY_BLOCK][NX_BLOCK];
extern double (*&due)[NY_BLOCK][NX_BLOCK];
extern double (*&duw)[NY_BLOCK][NX_BLOCK];
extern double (*&dmc)[NY_BLOCK][NX_BLOCK];
extern double (*&dmn)[NY_BLOCK][NX_BLOCK];
extern double (*&dms)[NY_BLOCK][NX_BLOCK];
extern double (*&dme)[NY_BLOCK][NX_BLOCK];
extern double (*&dmw)[NY_BLOCK][NX_BLOCK];
extern double (*&dum)[NY_BLOCK][NX_BLOCK];
extern double (*&ahf)[NY_BLOCK][NX_BLOCK];
extern double (*&amf)[NY_BLOCK][NX_BLOCK];

extern double &ah;
extern double &am;

} // CppHmixDel2
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_HMIX_DEL2_H_
