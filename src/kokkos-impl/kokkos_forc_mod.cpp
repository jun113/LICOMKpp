#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

namespace KokkosForcMod {

ViewDouble3D *p_v_su = nullptr;
ViewDouble3D *p_v_sv = nullptr;

ViewDouble3D *p_v_psa = nullptr;
ViewDouble3D *p_v_tsa = nullptr;
ViewDouble3D *p_v_sss = nullptr;
ViewDouble3D *p_v_swv = nullptr;

ViewDouble3D *p_v_uva    = nullptr;
ViewDouble3D *p_v_qar    = nullptr;
ViewDouble3D *p_v_cld    = nullptr;
ViewDouble3D *p_v_ddd    = nullptr;
ViewDouble3D *p_v_qqq    = nullptr;
ViewDouble3D *p_v_sst    = nullptr;
ViewDouble3D *p_v_nswv   = nullptr;
ViewDouble3D *p_v_dqdt   = nullptr;
ViewDouble3D *p_v_chloro = nullptr;
ViewDouble3D *p_v_lwv    = nullptr;
ViewDouble3D *p_v_seaice = nullptr;
ViewDouble3D *p_v_rain   = nullptr;
ViewDouble3D *p_v_snow   = nullptr;
ViewDouble3D *p_v_fresh  = nullptr;
ViewDouble3D *p_v_runoff = nullptr;
ViewDouble3D *p_v_lthf   = nullptr;
ViewDouble3D *p_v_sshf   = nullptr;

ViewDouble3D *p_v_ustar   = nullptr;
ViewDouble3D *p_v_buoytur = nullptr;
ViewDouble3D *p_v_buoysol = nullptr;

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

ViewDouble5D *p_v_restore = nullptr;

ViewDouble3D *p_v_tsf = nullptr;
ViewDouble3D *p_v_ssf = nullptr;

#ifdef TIDEMIX
ViewDouble3D *p_v_wave_dis = nullptr;
#endif // TIDEMIX

ViewFloat3D	*p_v_t10   = nullptr;
ViewFloat3D	*p_v_u10   = nullptr;
ViewFloat3D *p_v_v10   = nullptr;
ViewFloat3D *p_v_slp   = nullptr;
ViewFloat3D *p_v_q10   = nullptr;
ViewFloat3D *p_v_swhf  = nullptr;
ViewFloat3D *p_v_lwhf  = nullptr;
ViewFloat3D *p_v_precr = nullptr;
ViewFloat3D *p_v_precs = nullptr;
ViewFloat3D *p_v_rf    = nullptr;
ViewFloat3D *p_v_si    = nullptr;

ViewDouble1D *p_v_buffer    = nullptr;

ViewDouble2D *p_v_s_wx   = nullptr;
ViewDouble2D *p_v_s_wy   = nullptr;
ViewDouble2D *p_v_s_work = nullptr;

ViewDouble2D *p_v_tsa3   = nullptr;
ViewDouble2D *p_v_wspd3  = nullptr;
ViewDouble2D *p_v_wspdu3 = nullptr;
ViewDouble2D *p_v_wspdv3 = nullptr;
ViewDouble2D *p_v_psa3 = nullptr;
ViewDouble2D *p_v_qar3 = nullptr;
ViewDouble2D *p_v_swv3 = nullptr;
ViewDouble2D *p_v_lwv3 = nullptr;
ViewDouble2D *p_v_rain3 = nullptr;
ViewDouble2D *p_v_snow3 = nullptr;
ViewDouble2D *p_v_runoff3 = nullptr;
ViewDouble2D *p_v_seaice3 = nullptr;

ViewDouble2D *p_v_windx = nullptr;
ViewDouble2D *p_v_windy = nullptr;
ViewDouble2D *p_v_model_sst = nullptr;
ViewDouble2D *p_v_zz = nullptr;
ViewDouble2D *p_v_qs = nullptr;
ViewDouble2D *p_v_theta = nullptr;
ViewDouble2D *p_v_core_sensible = nullptr;
ViewDouble2D *p_v_core_latent = nullptr;
ViewDouble2D *p_v_core_tau = nullptr;
} // namespace KokkosForcMod
#endif // LICOM_ENABLE_KOKKOS
