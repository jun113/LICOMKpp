#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_TMP_VAR_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_TMP_VAR_H_

#include "def-undef.h"
namespace CppTmpVar {

// READYT
extern double *work1;
extern double *work2;

extern double *pp ;
extern double *ppa;
extern double *ppb;
extern double *ppc;

extern double *alpha;
extern double *beta ;

// READYC
extern double *wp12;
extern double *wp13;

extern double *riv1;
extern double *riv2;

#ifdef BCKMEX
extern double *diff_back;
extern double *diff_back_sh;
extern double *diff_back_nn;
#endif // BCKMEX

extern double *u_wface;
extern double *v_sface;
extern double *uv_ws_face;

#ifndef SMAG
extern double *hduk;
extern double *hdvk;

extern double *div_out;

extern double *am_factor;

#ifdef BIHAR
extern double *curl;
extern double *gradx1;
extern double *grady1;
extern double *gradx2;
extern double *grady2;

extern double *grad_xy_12;

extern double *am_factor;

extern double *cc;
extern double *d2uk;
extern double *d2vk;
#endif // BIHAR
#endif // SMAG

// TRACER
extern double *vtl_ori;

extern double *adv_tt;

extern double *at00;
extern double *atmax;
extern double *atmin;

extern double *at_00_max_min;

extern double *hdtk;
#ifdef BIHAR
extern double *dt2k;
#endif // BIHAR

extern int    *nn;
extern double *xs;

extern double *cn;
extern double *cs;
extern double *ce;
extern double *cw;

extern double *c_cnsew;
#ifdef ISO
#ifdef LDD97
extern double *f1;
extern double *f2;
#endif // LDD97
#endif // ISO

// BAROTR
extern double *gradx;
extern double *grady;

// POP Halo Update
extern double *halo_buffer;
} // namespace CppTmpVar

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_TMP_VAR_H_
