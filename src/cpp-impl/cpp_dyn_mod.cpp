#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_dyn_mod.h"

namespace CppDynMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::IMT_GLOBAL;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&ub)[MAX_BLOCKS_CLINIC][JMT][IMT]       = dyn_mod_mp_ub_;
double (&vb)[MAX_BLOCKS_CLINIC][JMT][IMT]       = dyn_mod_mp_vb_;
double (&ubp)[MAX_BLOCKS_CLINIC][JMT][IMT]      = dyn_mod_mp_ubp_;
double (&vbp)[MAX_BLOCKS_CLINIC][JMT][IMT]      = dyn_mod_mp_vbp_;
double (&h0p)[MAX_BLOCKS_CLINIC][JMT][IMT]      = dyn_mod_mp_h0p_;

double (&up)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = dyn_mod_mp_up_;
double (&vp)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = dyn_mod_mp_vp_;

double (&ws)[MAX_BLOCKS_CLINIC][KMP1][JMT][IMT] = dyn_mod_mp_ws_;

double (&h0l)[MAX_BLOCKS_CLINIC][JMT][IMT]      = dyn_mod_mp_h0l_;
double (&h0f)[MAX_BLOCKS_CLINIC][JMT][IMT]      = dyn_mod_mp_h0f_;
double (&h0bl)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_h0bl_;
double (&h0bf)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_h0bf_;

double (&utl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = dyn_mod_mp_utl_;
double (&utf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = dyn_mod_mp_utf_;
double (&vtl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = dyn_mod_mp_vtl_;
double (&vtf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = dyn_mod_mp_vtf_;

double (&sbcx)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_sbcx_;
double (&bbcx)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_bbcx_;
double (&sbcy)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_sbcy_;
double (&bbcy)[MAX_BLOCKS_CLINIC][JMT][IMT]     = dyn_mod_mp_bbcy_;

double (*(&buffer))[IMT_GLOBAL]                 = dyn_mod_mp_buffer_;

double (&h0)[MAX_BLOCKS_CLINIC][JMT][IMT]       = dyn_mod_mp_h0_;

double (&u)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = dyn_mod_mp_u_;
double (&v)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = dyn_mod_mp_v_;

// double (*(&gg))[KM][JMT][IMT]                   = dyn_mod_mp_gg_;
// double (*(&dlu))[KM][JMT][IMT]                  = dyn_mod_mp_dlu_;
// double (*(&dlv))[KM][JMT][IMT]                  = dyn_mod_mp_dlv_;

// double (*(&dlub))[JMT][IMT]                     = dyn_mod_mp_dlub_;
// double (*(&dlvb))[JMT][IMT]                     = dyn_mod_mp_dlvb_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&ub)[MAX_BLOCKS_CLINIC][JMT][IMT]       = __dyn_mod_MOD_ub;
double (&vb)[MAX_BLOCKS_CLINIC][JMT][IMT]       = __dyn_mod_MOD_vb;
double (&ubp)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __dyn_mod_MOD_ubp;
double (&vbp)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __dyn_mod_MOD_vbp;
double (&h0p)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __dyn_mod_MOD_h0p;

double (&up)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __dyn_mod_MOD_up;
double (&vp)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __dyn_mod_MOD_vp;

double (&ws)[MAX_BLOCKS_CLINIC][KMP1][JMT][IMT] = __dyn_mod_MOD_ws;

double (&h0l)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __dyn_mod_MOD_h0l;
double (&h0f)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __dyn_mod_MOD_h0f;
double (&h0bl)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_h0bl;
double (&h0bf)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_h0bf;

double (&utl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __dyn_mod_MOD_utl;
double (&utf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __dyn_mod_MOD_utf;
double (&vtl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __dyn_mod_MOD_vtl;
double (&vtf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __dyn_mod_MOD_vtf;

double (&sbcx)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_sbcx;
double (&bbcx)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_bbcx;
double (&sbcy)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_sbcy;
double (&bbcy)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __dyn_mod_MOD_bbcy;

double (*(&buffer))[IMT_GLOBAL]                 = __dyn_mod_MOD_buffer;

double (&h0)[MAX_BLOCKS_CLINIC][JMT][IMT]       = __dyn_mod_MOD_h0;

double (&u)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = __dyn_mod_MOD_u;
double (&v)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = __dyn_mod_MOD_v;

// double (*(&gg))[KM][JMT][IMT]                   = __dyn_mod_MOD_gg;
// double (*(&dlu))[KM][JMT][IMT]                  = __dyn_mod_MOD_dlu;
// double (*(&dlv))[KM][JMT][IMT]                  = __dyn_mod_MOD_dlv;

// double (*(&dlub))[JMT][IMT]                     = __dyn_mod_MOD_dlub;
// double (*(&dlvb))[JMT][IMT]                     = __dyn_mod_MOD_dlvb;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

double gg[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
double dlu[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
double dlv[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
double dlub[MAX_BLOCKS_CLINIC][JMT][IMT];
double dlvb[MAX_BLOCKS_CLINIC][JMT][IMT];
} // namespace CppDynMod
#endif // LICOM_ENABLE_FORTRAN
