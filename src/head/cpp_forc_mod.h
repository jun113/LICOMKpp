#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_FORC_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_FORC_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"
namespace CppForcMod {

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::S_JMT;
using CppParamMod::S_IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;

extern float *t10;
extern float *u10;
extern float *v10;
extern float *slp;
extern float *q10;
extern float *swhf;
extern float *lwhf;
extern float *precr;
extern float *precs;
extern float *rf;
extern float *si;

extern double* buffer;

extern double* tsa3;
extern double* wspdu3;
extern double* wspdv3;
extern double* psa3;
extern double* qar3;
extern double* swv3;
extern double* lwv3;
extern double* rain3;
extern double* snow3;
extern double* runoff3;
extern double* seaice3;

extern double* wspd3;

extern double* s_wx;
extern double* s_wy;
extern double* s_work;

extern double* windx;
extern double* windy;
extern double* model_sst;
extern double* zz;
extern double* qs;
extern double* theta;
extern double* core_sensible;
extern double* core_latent;
extern double* core_tau;

// extern double (*&su3)[12][JMT][IMT];
// extern double (*&sv3)[12][JMT][IMT];
// extern double (*&psa3)[12][JMT][IMT];
// extern double (&tsa3)[MAX_BLOCKS_CLINIC][JMT][IMT];
// extern double (*&qar3)[12][JMT][IMT];
// extern double (*&uva3)[12][JMT][IMT];
// extern double (*&swv3)[12][JMT][IMT];
// extern double (*&cld3)[12][JMT][IMT];
// extern double (*&sss3)[12][JMT][IMT];
// extern double (*&sst3)[12][JMT][IMT];
// extern double (*&nswv3)[12][JMT][IMT];
// extern double (*&dqdt3)[12][JMT][IMT];
// extern double (*&chloro3)[12][JMT][IMT];
// extern double (*&wspd3)[12][JMT][IMT];
// extern double (*&wspdu3)[12][JMT][IMT];
// extern double (*&wspdv3)[12][JMT][IMT];
// extern double (*&lwv3)[12][JMT][IMT];
// extern double (*&seaice3)[12][JMT][IMT];
// extern double (*&rain3)[12][JMT][IMT];
// extern double (*&snow3)[12][JMT][IMT];

extern double (&su)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&sv)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&psa)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&tsa)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&sss)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&swv)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&uva)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&qar)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&cld)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ddd)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&qqq)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&sst)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&nswv)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&dqdt)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&chloro)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&lwv)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&seaice)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&rain)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&snow)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&fresh)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&runoff)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&lthf)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&sshf)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (&ustar)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&buoytur)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&buoysol)[MAX_BLOCKS_CLINIC][JMT][IMT];

/*
extern double *(&su3_io)[];
extern double *(&sv3_io)[];
extern double *(&psa3_io)[];
extern double *(&tsa3_io)[];
extern double *(&qar3_io)[];
extern double *(&uva3_io)[];
extern double *(&swv3_io)[];
extern double *(&cld3_io)[];
extern double *(&sss3_io)[];
extern double *(&sst3_io)[];
extern double *(&nswv3_io)[];
extern double *(&dqdt3_io)[];
extern double *(&chloro3_io)[];
extern double *(&wspd3_io)[];
extern double *(&wspdu3_io)[];
extern double *(&wspdv3_io)[];
extern double *(&lwv3_io)[];
extern double *(&seaice_io)[];
extern double *(&rain_io)[];
extern double *(&snow_io)[];
extern double *(&fresh_io)[];
extern double *(&runoff_io)[];
*/

extern double (&restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double (&tsf)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ssf)[MAX_BLOCKS_CLINIC][JMT][IMT];

#ifdef TIDEMIX
extern double (&wave_dis)[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // TIDEMIX

/*
extern double (*(&t10)[S_JMT][S_IMT];
extern double (*(&u10)[S_JMT][S_IMT];
extern double (*(&v10)[S_JMT][S_IMT];
extern double (*(&slp)[S_JMT][S_IMT];
extern double (*(&q10)[S_JMT][S_IMT];
extern double (*(&swhf)[S_JMT][S_IMT];
extern double (*(&lwhf)[S_JMT][S_IMT];

extern double (*(&precr)[S_JMT][S_IMT];
extern double (*(&precs)[S_JMT][S_IMT];

extern double (*(&rf)[S_JMT][S_IMT];
extern double (*(&si)[S_JMT][S_IMT];
extern double (*(&buffer3d)[KM][JMT][IMT];
extern double (*(&w3d)[JMT_GLOBAL][IMT_GLOBAL];
*/
} // namespace CppForcMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_FORC_MOD_H_
