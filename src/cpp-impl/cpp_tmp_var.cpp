#include "../head/def-undef.h"
namespace CppTmpVar {
// READYT
double *work1 = nullptr;
double *work2 = nullptr;

double *pp  = nullptr;
double *ppa = nullptr;
double *ppb = nullptr;
double *ppc = nullptr;

double *alpha = nullptr;
double *beta  = nullptr;

// READYC
double *wp12 = nullptr;
double *wp13 = nullptr;

double *riv1 = nullptr;
double *riv2 = nullptr;

#ifdef BCKMEX
double *diff_back = nullptr;
double *diff_back_sh = nullptr;
double *diff_back_nn = nullptr;
#endif // BCKMEX

double *u_wface = nullptr;
double *v_sface = nullptr;
double *uv_ws_face = nullptr;

#ifndef SMAG
double *hduk = nullptr;
double *hdvk = nullptr;

double *div_out = nullptr;

#ifdef BIHAR
double *curl = nullptr;
double *gradx1 = nullptr;
double *grady1 = nullptr;
double *gradx2 = nullptr;
double *grady2 = nullptr;

double *grad_xy_12 = nullptr;

double *am_factor = nullptr;

double *cc = nullptr;
double *d2uk = nullptr;
double *d2vk = nullptr;
#endif // BIHAR
#endif // SMAG

// TRACER
double *vtl_ori = nullptr;

double *adv_tt = nullptr;

double *at00 = nullptr;
double *atmax = nullptr;
double *atmin = nullptr;

double *at_00_max_min = nullptr;

double *hdtk = nullptr;
#ifdef BIHAR
double *dt2k = nullptr;
#endif // BIHAR

int    *nn = nullptr;
double *xs = nullptr;

double *cn = nullptr;
double *cs = nullptr;
double *ce = nullptr;
double *cw = nullptr;

double *c_cnsew = nullptr;
#ifdef ISO
#ifdef LDD97
double *f1 = nullptr;
double *f2 = nullptr;
#endif // LDD97
#endif // ISO

// BAROTR
double *gradx = nullptr;
double *grady = nullptr;

// POP Halo Update
double *halo_buffer = nullptr;
} // namespace CppTmpVar