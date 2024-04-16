#ifndef LICOM3_KOKKOS_SRC_HAED_CPP_PARAM_MOD_H_
#define LICOM3_KOKKOS_SRC_HAED_CPP_PARAM_MOD_H_

#include "def-undef.h"
namespace CppParamMod {

constexpr int LICOM_BLOCKSIZEX  = BLCKX;
constexpr int LICOM_BLOCKSIZEY  = BLCKY;

constexpr int MAX_BLOCKS_CLINIC = MXBLCKS;
constexpr int MAX_BLOCKS_TROPIC = MXBLCKS;

constexpr int NGHOST            = 2;
constexpr int NX_BLOCK          = LICOM_BLOCKSIZEX + 2 * NGHOST;
constexpr int NY_BLOCK          = LICOM_BLOCKSIZEY + 2 * NGHOST;

constexpr int JMT_GLOBAL        = NJMT;
constexpr int JMM_GLOBAL        = JMT_GLOBAL - 1;
constexpr int IMT_GLOBAL        = NIMT;
constexpr int KM                = NKM;

constexpr int NUM_OVERLAP       = 2;

constexpr int JST               = 1;
constexpr int JSM               = JST + 1;

constexpr int JET               = NY_BLOCK;
constexpr int JEM               = JET - 1;

constexpr int JMT               = NY_BLOCK;
constexpr int IMT               = NX_BLOCK;

constexpr int IMM_GLOBAL        = IMT_GLOBAL - 1;
constexpr int IMM               = IMT - 1;
constexpr int JMM               = JMT - 1;
constexpr int KMP1              = KM + 1;
constexpr int KMM1              = KM - 1;

constexpr int NTRA              = 2;

extern int &ierr;
extern int &mytid;

extern int &jj_start;
extern int &jj_end;

extern int &my_task;
extern int &master_task;

extern int &am_hor;
extern int &an_hor;

extern int &am_bihar;
extern int &an_bihar;

extern int &num_cpl;

constexpr int S_IMT             = 640;
constexpr int S_JMT             = 320;
} // namespace CppParamMod

#endif // LICOM3_KOKKOS_SRC_HAED_CPP_PARAM_MOD_H_
