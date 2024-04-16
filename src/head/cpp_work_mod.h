#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_WORK_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_WORK_MOD_H_

#include "cpp_param_mod.h"
namespace CppWorkMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;

extern double (&pxb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&pyb)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&pax)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&pay)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&whx)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&why)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&wgp)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&wka)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*&work_1)[KM][JMT][IMT];
// extern double (*&work_2)[KM][JMT][IMT];
// extern double (*&work_3)[KM+1][JMT][IMT];
// extern double (*&temp)[KM][JMT][IMT];

// extern double (*&uk)[KM][JMT][IMT];
// extern double (*&vk)[KM][JMT][IMT];

extern double (&work)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&wkk)[KMP1];

// extern double (*&wkb)[KM][JMT][IMT];
// extern double (*&wkc)[KM][JMT][IMT];
// extern double (*&wkd)[KM][JMT][IMT];

// extern double (*&tf)[KM][JMT][IMT];
// extern double (*&stf)[JMT][IMT];

// extern float  (*&buffer_real4)[IMT_GLOBAL];

// extern double (*&work1_g)[IMT_GLOBAL];
// extern double (*&work2_g)[IMT_GLOBAL];
// extern double (*&work3_g)[IMT_GLOBAL];

extern double* work_1;
extern double* work_2;
extern double* work_3;
extern double* temp;
extern double* uk;
extern double* vk;
extern double* wkb;
extern double* wkc;
extern double* wkd;
extern double* tf;
extern double* stf;
extern double* work1_g;
extern double* work2_g;
} // namespace CppWorkMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_WORK_MOD_H_
