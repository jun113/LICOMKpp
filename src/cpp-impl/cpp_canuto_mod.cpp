#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/fortran_canuto_mod.h"
namespace CppCanutoMod {

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double &eps_bot__under                 = canuto_mod_mp_eps_bot__under_;
                                     
double (&and2on2a1)[2*MT+1]            = canuto_mod_mp_and2on2a1_;
double (&amtaun2a1)[2*MT+1]            = canuto_mod_mp_amtaun2a1_;
                                     
double &dand2on2                       = canuto_mod_mp_dand2on2_;
                                     
double (&sma1)[2*MT+1]                 = canuto_mod_mp_sma1_;
double (&sha1)[2*MT+1]                 = canuto_mod_mp_sha1_;
double (&ssa1)[2*MT+1]                 = canuto_mod_mp_ssa1_;
                                     
double &rri                            = canuto_mod_mp_rri_;
double &rnd2on2                        = canuto_mod_mp_rnd2on2_;
double &dri                            = canuto_mod_mp_dri_;
double &deltheta_r                     = canuto_mod_mp_deltheta_r_;
double &b1                             = canuto_mod_mp_b1_;
double &theta_rcrp                     = canuto_mod_mp_theta_rcrp_;
double &theta_rcrn                     = canuto_mod_mp_theta_rcrn_;
double &ako                            = canuto_mod_mp_ako_;
double &back_l_0                       = canuto_mod_mp_back_l_0_;
                                     
double (&rib)[2*MT+1]                  = canuto_mod_mp_rib_;
double (&ridb)[2*MT+1]                 = canuto_mod_mp_ridb_;
                                     
double (&slq2b)[2*MT+1][2*MT+1]        = canuto_mod_mp_slq2b_;
double (&smb)[2*MT+1][2*MT+1]          = canuto_mod_mp_smb_;
double (&shb)[2*MT+1][2*MT+1]          = canuto_mod_mp_shb_;
double (&ssb)[2*MT+1][2*MT+1]          = canuto_mod_mp_ssb_;
                                     
int    (&irimax)[2*MT+1]               = canuto_mod_mp_irimax_;
                                     
double (&sisamax)[4*N_THETA_R_OCT+1]   = canuto_mod_mp_sisamax_;
double (&ra_rmax)[4*N_THETA_R_OCT+1]   = canuto_mod_mp_ra_rmax_;
double (&c_y_r0)[4*N_THETA_R_OCT+1]    = canuto_mod_mp_c_y_r0_;
double (&back_ra_r)[4*N_THETA_R_OCT+1] = canuto_mod_mp_back_ra_r_;
double (&sm_r1)[4*N_THETA_R_OCT+1]     = canuto_mod_mp_sm_r1_;
double (&sh_r1)[4*N_THETA_R_OCT+1]     = canuto_mod_mp_sh_r1_;
double (&ss_r1)[4*N_THETA_R_OCT+1]     = canuto_mod_mp_ss_r1_;
double (&slq2_r1)[4*N_THETA_R_OCT+1]   = canuto_mod_mp_slq2_r1_;

double (&ria)[NTBL]                    = canuto_mod_mp_ria_;
double (&slq2a)[NTBL]                  = canuto_mod_mp_slq2a_;
double (&sma)[NTBL]                    = canuto_mod_mp_sma_;
double (&sha)[NTBL]                    = canuto_mod_mp_sha_;

int &ifirst                            = canuto_mod_mp_ifirst_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double &eps_bot__under                 = __canuto_mod_MOD_eps_bot__under;
                                     
double (&and2on2a1)[2*MT+1]            = __canuto_mod_MOD_and2on2a1;
double (&amtaun2a1)[2*MT+1]            = __canuto_mod_MOD_amtaun2a1;
                                     
double &dand2on2                       = __canuto_mod_MOD_dand2on2;
                                     
double (&sma1)[2*MT+1]                 = __canuto_mod_MOD_sma1;
double (&sha1)[2*MT+1]                 = __canuto_mod_MOD_sha1;
double (&ssa1)[2*MT+1]                 = __canuto_mod_MOD_ssa1;
                                     
double &rri                            = __canuto_mod_MOD_rri;
double &rnd2on2                        = __canuto_mod_MOD_rnd2on2;
double &dri                            = __canuto_mod_MOD_dri;
double &deltheta_r                     = __canuto_mod_MOD_deltheta_r;
double &b1                             = __canuto_mod_MOD_b1;
double &theta_rcrp                     = __canuto_mod_MOD_theta_rcrp;
double &theta_rcrn                     = __canuto_mod_MOD_theta_rcrn;
double &ako                            = __canuto_mod_MOD_ako;
double &back_l_0                       = __canuto_mod_MOD_back_l_0;
                                     
double (&rib)[2*MT+1]                  = __canuto_mod_MOD_rib;
double (&ridb)[2*MT+1]                 = __canuto_mod_MOD_ridb;
                                     
double (&slq2b)[2*MT+1][2*MT+1]        = __canuto_mod_MOD_slq2b;
double (&smb)[2*MT+1][2*MT+1]          = __canuto_mod_MOD_smb;
double (&shb)[2*MT+1][2*MT+1]          = __canuto_mod_MOD_shb;
double (&ssb)[2*MT+1][2*MT+1]          = __canuto_mod_MOD_ssb;
                                     
int    (&irimax)[2*MT+1]               = __canuto_mod_MOD_irimax;
                                     
double (&sisamax)[4*N_THETA_R_OCT+1]   = __canuto_mod_MOD_sisamax;
double (&ra_rmax)[4*N_THETA_R_OCT+1]   = __canuto_mod_MOD_ra_rmax;
double (&c_y_r0)[4*N_THETA_R_OCT+1]    = __canuto_mod_MOD_c_y_r0;
double (&back_ra_r)[4*N_THETA_R_OCT+1] = __canuto_mod_MOD_back_ra_r;
double (&sm_r1)[4*N_THETA_R_OCT+1]     = __canuto_mod_MOD_sm_r1;
double (&sh_r1)[4*N_THETA_R_OCT+1]     = __canuto_mod_MOD_sh_r1;
double (&ss_r1)[4*N_THETA_R_OCT+1]     = __canuto_mod_MOD_ss_r1;
double (&slq2_r1)[4*N_THETA_R_OCT+1]   = __canuto_mod_MOD_slq2_r1;

double (&ria)[NTBL]                    = __canuto_mod_MOD_ria;
double (&slq2a)[NTBL]                  = __canuto_mod_MOD_slq2a;
double (&sma)[NTBL]                    = __canuto_mod_MOD_sma;
double (&sha)[NTBL]                    = __canuto_mod_MOD_sha;

int &ifirst                            = __canuto_mod_MOD_ifirst;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppCanutoMod
#endif // LICOM_ENABLE_FORTRAN
