#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_forc_mod.h"
namespace CppForcMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::S_IMT;
using CppParamMod::S_JMT;
using CppParamMod::NTRA;

float *t10   = nullptr;
float *u10   = nullptr;
float *v10   = nullptr;
float *slp   = nullptr;
float *q10   = nullptr;
float *swhf  = nullptr;
float *lwhf  = nullptr;
float *precr = nullptr;
float *precs = nullptr;
float *rf    = nullptr;
float *si    = nullptr;

double* buffer = nullptr;

double* tsa3    = nullptr;
double* wspdu3  = nullptr;
double* wspdv3  = nullptr;
double* psa3    = nullptr;
double* qar3    = nullptr;
double* swv3    = nullptr;
double* lwv3    = nullptr;
double* rain3   = nullptr;
double* snow3   = nullptr;
double* runoff3 = nullptr;
double* seaice3 = nullptr;

double* wspd3   = nullptr;

double* s_wx   = nullptr;
double* s_wy   = nullptr;
double* s_work = nullptr;

double* windx         = nullptr;
double* windy         = nullptr;
double* model_sst     = nullptr;
double* zz            = nullptr;
double* qs            = nullptr;
double* theta         = nullptr;
double* core_sensible = nullptr;
double* core_latent   = nullptr;
double* core_tau      = nullptr;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
// double (*(&su3))[12][JMT][IMT]                       = forc_mod_mp_su3_;
// double (*(&sv3))[12][JMT][IMT]                       = forc_mod_mp_sv3_;
// double (*(&psa3))[12][JMT][IMT]                      = forc_mod_mp_psa3_;
// double (&tsa3)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_tsa3_;
// double (*(&qar3))[12][JMT][IMT]                      = forc_mod_mp_qar3_;
// double (*(&uva3))[12][JMT][IMT]                      = forc_mod_mp_uva3_;
// double (*(&swv3))[12][JMT][IMT]                      = forc_mod_mp_swv3_;
// double (*(&cld3))[12][JMT][IMT]                      = forc_mod_mp_cld3_;
// double (*(&sss3))[12][JMT][IMT]                      = forc_mod_mp_sss3_;
// double (*(&sst3))[12][JMT][IMT]                      = forc_mod_mp_sst3_;
// double (*(&nswv3))[12][JMT][IMT]                     = forc_mod_mp_nswv3_;
// double (*(&dqdt3))[12][JMT][IMT]                     = forc_mod_mp_dqdt3_;
// double (*(&chloro3))[12][JMT][IMT]                   = forc_mod_mp_chloro3_;
// double (*(&wspd3))[12][JMT][IMT]                     = forc_mod_mp_wspd3_;
// double (*(&wspdu3))[12][JMT][IMT]                    = forc_mod_mp_wspdu3_;
// double (*(&wspdv3))[12][JMT][IMT]                    = forc_mod_mp_wspdv3_;
// double (*(&lwv3))[12][JMT][IMT]                      = forc_mod_mp_lwv3_;
// double (*(&seaice3))[12][JMT][IMT]                   = forc_mod_mp_seaice3_;
// double (*(&rain3))[12][JMT][IMT]                     = forc_mod_mp_rain3_;
// double (*(&snow3))[12][JMT][IMT]                     = forc_mod_mp_snow3_;

double (&su)[MAX_BLOCKS_CLINIC][JMT][IMT]            = forc_mod_mp_su_;
double (&sv)[MAX_BLOCKS_CLINIC][JMT][IMT]            = forc_mod_mp_sv_;
double (&psa)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_psa_;
double (&tsa)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_tsa_;
double (&sss)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_sss_;
double (&swv)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_swv_;
double (&uva)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_uva_;
double (&qar)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_qar_;
double (&cld)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_cld_;
double (&ddd)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_ddd_;
double (&qqq)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_qqq_;
double (&sst)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_sst_;
double (&nswv)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_nswv_;
double (&dqdt)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_dqdt_;
double (&chloro)[MAX_BLOCKS_CLINIC][JMT][IMT]        = forc_mod_mp_chloro_;
double (&lwv)[MAX_BLOCKS_CLINIC][JMT][IMT]           = forc_mod_mp_lwv_;
double (&seaice)[MAX_BLOCKS_CLINIC][JMT][IMT]        = forc_mod_mp_seaice_;
double (&rain)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_rain_;
double (&snow)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_snow_;
double (&fresh)[MAX_BLOCKS_CLINIC][JMT][IMT]         = forc_mod_mp_fresh_;
double (&runoff)[MAX_BLOCKS_CLINIC][JMT][IMT]        = forc_mod_mp_runoff_;
double (&lthf)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_lthf_;
double (&sshf)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_sshf_;

double (&ustar)[MAX_BLOCKS_CLINIC][JMT][IMT]         = forc_mod_mp_ustar_;
double (&buoytur)[MAX_BLOCKS_CLINIC][JMT][IMT]       = forc_mod_mp_buoytur_;
double (&buoysol)[MAX_BLOCKS_CLINIC][JMT][IMT]       = forc_mod_mp_buoysol_;

/*
double *(&su3_io)[]                                  = forc_mod_mp_su3_io_;
double *(&sv3_io)[]                                  = forc_mod_mp_sv3_io_;
double *(&psa3_io)[]                                 = forc_mod_mp_psa3_io_;
double *(&tsa3_io)[]                                 = forc_mod_mp_tsa3_io_;
double *(&qar3_io)[]                                 = forc_mod_mp_qar3_io_;
double *(&uva3_io)[]                                 = forc_mod_mp_uva3_io_;
double *(&swv3_io)[]                                 = forc_mod_mp_swv3_io_;
double *(&cld3_io)[]                                 = forc_mod_mp_cld3_io_;
double *(&sss3_io)[]                                 = forc_mod_mp_sss3_io_;
double *(&sst3_io)[]                                 = forc_mod_mp_sst3_io_;
double *(&nswv3_io)[]                                = forc_mod_mp_nswv3_io_;
double *(&dqdt3_io)[]                                = forc_mod_mp_dqdt3_io_;
double *(&chloro3_io)[]                              = forc_mod_mp_chloro3_io_;
double *(&wspd3_io)[]                                = forc_mod_mp_wspd3_io_;
double *(&wspdu3_io)[]                               = forc_mod_mp_wspdu_io_;
double *(&wspdv3_io)[]                               = forc_mod_mp_wspdv_io_;
double *(&lwv3_io)[]                                 = forc_mod_mp_lwv3_io_;
double *(&seaice_io)[]                               = forc_mod_mp_seaice_io_;
double *(&rain_io)[]                                 = forc_mod_mp_rain_io_;
double *(&snow_io)[]                                 = forc_mod_mp_snow_io_;
double *(&fresh_io)[]                                = forc_mod_mp_fresh_io_;
double *(&runoff_io)[]                               = forc_mod_mp_runoff_io_;
*/

double (&restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = forc_mod_mp_restore_;
double (&tsf)[MAX_BLOCKS_CLINIC][JMT][IMT]               = forc_mod_mp_tsf_;
double (&ssf)[MAX_BLOCKS_CLINIC][JMT][IMT]               = forc_mod_mp_ssf_;
                                                    
#ifdef TIDEMIX                                      
double (&wave_dis)[MAX_BLOCKS_CLINIC][JMT][IMT]          = forc_mod_mp_wave_dis_;
#endif // TIDEMIX                                   
                                                    
/*
double (*(&t10))[S_JMT][S_IMT]                           = forc_mod_mp_t10_;
double (*(&u10))[S_JMT][S_IMT]                           = forc_mod_mp_u10_;
double (*(&v10))[S_JMT][S_IMT]                           = forc_mod_mp_v10_;
double (*(&slp))[S_JMT][S_IMT]                           = forc_mod_mp_slp_;
double (*(&q10))[S_JMT][S_IMT]                           = forc_mod_mp_q10_;
double (*(&swhf))[S_JMT][S_IMT]                          = forc_mod_mp_swhf_;
double (*(&lwhf))[S_JMT][S_IMT]                          = forc_mod_mp_lwhf_;
                                                    
double (*(&precr))[S_JMT][S_IMT]                         = forc_mod_mp_precr_;
double (*(&precs))[S_JMT][S_IMT]                         = forc_mod_mp_precs_;
                                                    
double (*(&rf))[S_JMT][S_IMT]                            = forc_mod_mp_rf_;
double (*(&si))[S_JMT][S_IMT]                            = forc_mod_mp_si_;
double (*(&buffer3d))[KM][JMT][IMT]                      = forc_mod_mp_buffer3d_;
double (*(&w3d))[JMT_GLOBAL][IMT_GLOBAL]                 = forc_mod_mp_w3d_;
*/
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
// double (*(&su3))[12][JMT][IMT]                       = __forc_mod_MOD_su3;
// double (*(&sv3))[12][JMT][IMT]                       = __forc_mod_MOD_sv3;
// double (*(&psa3))[12][JMT][IMT]                      = __forc_mod_MOD_psa3;
// double (&tsa3)[MAX_BLOCKS_CLINIC][JMT][IMT]                      = __forc_mod_MOD_tsa3;
// double (*(&qar3))[12][JMT][IMT]                      = __forc_mod_MOD_qar3;
// double (*(&uva3))[12][JMT][IMT]                      = __forc_mod_MOD_uva3;
// double (*(&swv3))[12][JMT][IMT]                      = __forc_mod_MOD_swv3;
// double (*(&cld3))[12][JMT][IMT]                      = __forc_mod_MOD_cld3;
// double (*(&sss3))[12][JMT][IMT]                      = __forc_mod_MOD_sss3;
// double (*(&sst3))[12][JMT][IMT]                      = __forc_mod_MOD_sst3;
// double (*(&nswv3))[12][JMT][IMT]                     = __forc_mod_MOD_nswv3;
// double (*(&dqdt3))[12][JMT][IMT]                     = __forc_mod_MOD_dqdt3;
// double (*(&chloro3))[12][JMT][IMT]                   = __forc_mod_MOD_chloro3;
// double (*(&wspd3))[12][JMT][IMT]                     = __forc_mod_MOD_wspd3;
// double (*(&wspdu3))[12][JMT][IMT]                    = __forc_mod_MOD_wspdu3;
// double (*(&wspdv3))[12][JMT][IMT]                    = __forc_mod_MOD_wspdv3;
// double (*(&lwv3))[12][JMT][IMT]                      = __forc_mod_MOD_lwv3;
// double (*(&seaice3))[12][JMT][IMT]                   = __forc_mod_MOD_seaice3;
// double (*(&rain3))[12][JMT][IMT]                     = __forc_mod_MOD_rain3;
// double (*(&snow3))[12][JMT][IMT]                     = __forc_mod_MOD_snow3;

double (&su)[MAX_BLOCKS_CLINIC][JMT][IMT]            = __forc_mod_MOD_su;
double (&sv)[MAX_BLOCKS_CLINIC][JMT][IMT]            = __forc_mod_MOD_sv;
double (&psa)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_psa;
double (&tsa)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_tsa;
double (&sss)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_sss;
double (&swv)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_swv;
double (&uva)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_uva;
double (&qar)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_qar;
double (&cld)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_cld;
double (&ddd)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_ddd;
double (&qqq)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_qqq;
double (&sst)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_sst;
double (&nswv)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_nswv;
double (&dqdt)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_dqdt;
double (&chloro)[MAX_BLOCKS_CLINIC][JMT][IMT]        = __forc_mod_MOD_chloro;
double (&lwv)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __forc_mod_MOD_lwv;
double (&seaice)[MAX_BLOCKS_CLINIC][JMT][IMT]        = __forc_mod_MOD_seaice;
double (&rain)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_rain;
double (&snow)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_snow;
double (&fresh)[MAX_BLOCKS_CLINIC][JMT][IMT]         = __forc_mod_MOD_fresh;
double (&runoff)[MAX_BLOCKS_CLINIC][JMT][IMT]        = __forc_mod_MOD_runoff;
double (&lthf)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_lthf;
double (&sshf)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_sshf;

double (&ustar)[MAX_BLOCKS_CLINIC][JMT][IMT]         = __forc_mod_MOD_ustar;
double (&buoytur)[MAX_BLOCKS_CLINIC][JMT][IMT]       = __forc_mod_MOD_buoytur;
double (&buoysol)[MAX_BLOCKS_CLINIC][JMT][IMT]       = __forc_mod_MOD_buoysol;

/*
double *(&su3_io)[]                                  = __forc_mod_MOD_su3_io;
double *(&sv3_io)[]                                  = __forc_mod_MOD_sv3_io;
double *(&psa3_io)[]                                 = __forc_mod_MOD_psa3_io;
double *(&tsa3_io)[]                                 = __forc_mod_MOD_tsa3_io;
double *(&qar3_io)[]                                 = __forc_mod_MOD_qar3_io;
double *(&uva3_io)[]                                 = __forc_mod_MOD_uva3_io;
double *(&swv3_io)[]                                 = __forc_mod_MOD_swv3_io;
double *(&cld3_io)[]                                 = __forc_mod_MOD_cld3_io;
double *(&sss3_io)[]                                 = __forc_mod_MOD_sss3_io;
double *(&sst3_io)[]                                 = __forc_mod_MOD_sst3_io;
double *(&nswv3_io)[]                                = __forc_mod_MOD_nswv3_io;
double *(&dqdt3_io)[]                                = __forc_mod_MOD_dqdt3_io;
double *(&chloro3_io)[]                              = __forc_mod_MOD_chloro3_io;
double *(&wspd3_io)[]                                = __forc_mod_MOD_wspd3_io;
double *(&wspdu3_io)[]                               = __forc_mod_MOD_wspdu_io;
double *(&wspdv3_io)[]                               = __forc_mod_MOD_wspdv_io;
double *(&lwv3_io)[]                                 = __forc_mod_MOD_lwv3_io;
double *(&seaice_io)[]                               = __forc_mod_MOD_seaice_io;
double *(&rain_io)[]                                 = __forc_mod_MOD_rain_io;
double *(&snow_io)[]                                 = __forc_mod_MOD_snow_io;
double *(&fresh_io)[]                                = __forc_mod_MOD_fresh_io;
double *(&runoff_io)[]                               = __forc_mod_MOD_runoff_io;
*/

double (&restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __forc_mod_MOD_restore;
double (&tsf)[MAX_BLOCKS_CLINIC][JMT][IMT]               = __forc_mod_MOD_tsf;
double (&ssf)[MAX_BLOCKS_CLINIC][JMT][IMT]               = __forc_mod_MOD_ssf;
                                                    
#ifdef TIDEMIX                                      
double (&wave_dis)[MAX_BLOCKS_CLINIC][JMT][IMT]          = __forc_mod_MOD_wave_dis;
#endif // TIDEMIX                                   
                                                    
/*
double (*(&t10))[S_JMT][S_IMT]                           = __forc_mod_MOD_t10;
double (*(&u10))[S_JMT][S_IMT]                           = __forc_mod_MOD_u10;
double (*(&v10))[S_JMT][S_IMT]                           = __forc_mod_MOD_v10;
double (*(&slp))[S_JMT][S_IMT]                           = __forc_mod_MOD_slp;
double (*(&q10))[S_JMT][S_IMT]                           = __forc_mod_MOD_q10;
double (*(&swhf))[S_JMT][S_IMT]                          = __forc_mod_MOD_swhf;
double (*(&lwhf))[S_JMT][S_IMT]                          = __forc_mod_MOD_lwhf;
                                                    
double (*(&precr))[S_JMT][S_IMT]                         = __forc_mod_MOD_precr;
double (*(&precs))[S_JMT][S_IMT]                         = __forc_mod_MOD_precs;
                                                    
double (*(&rf))[S_JMT][S_IMT]                            = __forc_mod_MOD_rf;
double (*(&si))[S_JMT][S_IMT]                            = __forc_mod_MOD_si;
double (*(&buffer3d))[KM][JMT][IMT]                      = __forc_mod_MOD_buffer3d;
double (*(&w3d))[JMT_GLOBAL][IMT_GLOBAL]                 = __forc_mod_MOD_w3d;
*/
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppForcMod
#endif // LICOM_ENABLE_FORTRAN
