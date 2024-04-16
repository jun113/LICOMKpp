#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_FORC_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_FORC_MOD_H_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosForcMod {
/*
double *(&su3)[12][JMT][IMT]                         = forc_mod_mp_su3_;
double *(&sv3)[12][JMT][IMT]                         = forc_mod_mp_sv3_;
double *(&psa3)[12][JMT][IMT]                        = forc_mod_mp_psa3_;
double *(&tsa3)[12][JMT][IMT]                        = forc_mod_mp_tsa3_;
double *(&qar3)[12][JMT][IMT]                        = forc_mod_mp_qar3_;
double *(&uva3)[12][JMT][IMT]                        = forc_mod_mp_uva3_;
double *(&swv3)[12][JMT][IMT]                        = forc_mod_mp_swv3_;
double *(&cld3)[12][JMT][IMT]                        = forc_mod_mp_cld3_;
double *(&sss3)[12][JMT][IMT]                        = forc_mod_mp_sss3_;
double *(&sst3)[12][JMT][IMT]                        = forc_mod_mp_sst3_;
double *(&nswv3)[12][JMT][IMT]                       = forc_mod_mp_nswv3_;
double *(&dqdt3)[12][JMT][IMT]                       = forc_mod_mp_dqdt3_;
double *(&chloro3)[12][JMT][IMT]                     = forc_mod_mp_chloro3_;
double *(&wspd3)[12][JMT][IMT]                       = forc_mod_mp_wspd3_;
double *(&wspdu3)[12][JMT][IMT]                      = forc_mod_mp_wspdu3_;
double *(&wspdv3)[12][JMT][IMT]                      = forc_mod_mp_wspdv3_;
double *(&lwv3)[12][JMT][IMT]                        = forc_mod_mp_lwv3_;
double *(&seaice3)[12][JMT][IMT]                     = forc_mod_mp_seaice3_;
double *(&rain3)[12][JMT][IMT]                       = forc_mod_mp_rain3_;
double *(&snow3)[12][JMT][IMT]                       = forc_mod_mp_snow3_;
*/

extern ViewDouble3D *p_v_su;
extern ViewDouble3D *p_v_sv;

extern ViewDouble3D *p_v_psa;
extern ViewDouble3D *p_v_tsa;
extern ViewDouble3D *p_v_sss;
extern ViewDouble3D *p_v_swv;

extern ViewDouble3D *p_v_uva;
extern ViewDouble3D *p_v_qar;
extern ViewDouble3D *p_v_cld;
extern ViewDouble3D *p_v_ddd;
extern ViewDouble3D *p_v_qqq;
extern ViewDouble3D *p_v_sst;
extern ViewDouble3D *p_v_nswv;
extern ViewDouble3D *p_v_dqdt;
extern ViewDouble3D *p_v_chloro;
extern ViewDouble3D *p_v_lwv;
extern ViewDouble3D *p_v_seaice;
extern ViewDouble3D *p_v_rain;
extern ViewDouble3D *p_v_snow;
extern ViewDouble3D *p_v_fresh;
extern ViewDouble3D *p_v_runoff;
extern ViewDouble3D *p_v_lthf;
extern ViewDouble3D *p_v_sshf;

extern ViewDouble3D *p_v_ustar;
extern ViewDouble3D *p_v_buoytur;
extern ViewDouble3D *p_v_buoysol;

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

extern ViewDouble5D *p_v_restore;

extern ViewDouble3D *p_v_tsf;
extern ViewDouble3D *p_v_ssf;

#ifdef TIDEMIX
extern ViewDouble3D *p_v_wave_dis;
#endif // TIDEMIX

extern ViewFloat3D *p_v_t10;
extern ViewFloat3D *p_v_u10;
extern ViewFloat3D *p_v_v10;
extern ViewFloat3D *p_v_slp;
extern ViewFloat3D *p_v_q10;
extern ViewFloat3D *p_v_swhf;
extern ViewFloat3D *p_v_lwhf;
extern ViewFloat3D *p_v_precr;
extern ViewFloat3D *p_v_precs;
extern ViewFloat3D *p_v_rf;
extern ViewFloat3D *p_v_si;

extern ViewDouble1D *p_v_buffer;

extern ViewDouble2D *p_v_s_wx;
extern ViewDouble2D *p_v_s_wy;
extern ViewDouble2D *p_v_s_work;

extern ViewDouble2D *p_v_tsa3;
extern ViewDouble2D *p_v_wspd3;
extern ViewDouble2D *p_v_wspdu3;
extern ViewDouble2D *p_v_wspdv3;
extern ViewDouble2D *p_v_psa3;
extern ViewDouble2D *p_v_qar3;
extern ViewDouble2D *p_v_swv3;
extern ViewDouble2D *p_v_lwv3;
extern ViewDouble2D *p_v_rain3;
extern ViewDouble2D *p_v_snow3;
extern ViewDouble2D *p_v_runoff3;
extern ViewDouble2D *p_v_seaice3;

extern ViewDouble2D *p_v_windx;
extern ViewDouble2D *p_v_windy;
extern ViewDouble2D *p_v_model_sst;
extern ViewDouble2D *p_v_zz;
extern ViewDouble2D *p_v_qs;
extern ViewDouble2D *p_v_theta;
extern ViewDouble2D *p_v_core_sensible;
extern ViewDouble2D *p_v_core_latent;
extern ViewDouble2D *p_v_core_tau;

} // namespace KokkosForcMod
#endif // LICOM_ENABLE_KOKKOS
#endif // LICOM3_KOKKOS_SRC_HEAD_KOKKOS_FORC_MOD_H_
