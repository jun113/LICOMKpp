#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_CANUTO_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_CANUTO_MOD_H_

constexpr int      NMODEL           =   1;
constexpr int      NTBL             =   501;
constexpr double   RI0              = - 4.e0;
constexpr double   E                =   2.71828182845904509e0;
constexpr double   PI_CONST         =   3.14159265358979312e0;
constexpr int      IFBACK           =   5;
constexpr int      IFSALI           =   1;
constexpr bool     IFSHEARMIN       =   true;
constexpr double   S2MIN            =   1.e-14;
constexpr bool     IFZEROSHEAR      =   true;
constexpr int      ICONDEAR         =   0;
constexpr int      IFEPSON2         =   2;
constexpr double   EPSON2__         =   0.288e0;
constexpr int      IFBOTENHANCE     =   0;
constexpr double   EPS_BOT0         =   2.e-5;
constexpr double   SCALE_BOT        =   5.e4;
constexpr int      IFDEEPLAT        =   1;
constexpr double   EPLATIDEPENDMIN  =   7.e-2;
constexpr int      IFCHECKBOTTOMEPS =   0;
constexpr int      IFRAFGMAX        =   1;
constexpr int      IFSALBACK        =   5;

constexpr int      NEXTRTBL0        =   62;
constexpr int      IFEXPABSTABLE    =   1;
constexpr int      IFAST            =   1;
constexpr int      IFASTEXPABS      =   IFAST * IFEXPABSTABLE;
constexpr int      NEXTRTBL1        =   1000;
constexpr int      NEXTRTBL         =   NEXTRTBL0 + IFEXPABSTABLE * NEXTRTBL1;
constexpr int      NPOSAPPROX       =   101;
constexpr int      MT0              =   NTBL - NPOSAPPROX;
constexpr int      MT               =   MT0 + NEXTRTBL;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern    double   canuto_mod_mp_eps_bot__under_;
extern    double   canuto_mod_mp_and2on2a1_[2*MT+1];
extern    double   canuto_mod_mp_amtaun2a1_[2*MT+1];
extern    double   canuto_mod_mp_dand2on2_;
extern    double   canuto_mod_mp_sma1_[2*MT+1];
extern    double   canuto_mod_mp_sha1_[2*MT+1];
extern    double   canuto_mod_mp_ssa1_[2*MT+1];
extern    double   canuto_mod_mp_rri_;
extern    double   canuto_mod_mp_rnd2on2_;
extern    double   canuto_mod_mp_dri_;
extern    double   canuto_mod_mp_deltheta_r_;
extern    double   canuto_mod_mp_b1_;
extern    double   canuto_mod_mp_theta_rcrp_;
extern    double   canuto_mod_mp_theta_rcrn_;
extern    double   canuto_mod_mp_ako_;
extern    double   canuto_mod_mp_back_l_0_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern    double   __canuto_mod_MOD_eps_bot__under;
extern    double   __canuto_mod_MOD_and2on2a1[2*MT+1];
extern    double   __canuto_mod_MOD_amtaun2a1[2*MT+1];
extern    double   __canuto_mod_MOD_dand2on2;
extern    double   __canuto_mod_MOD_sma1[2*MT+1];
extern    double   __canuto_mod_MOD_sha1[2*MT+1];
extern    double   __canuto_mod_MOD_ssa1[2*MT+1];
extern    double   __canuto_mod_MOD_rri;
extern    double   __canuto_mod_MOD_rnd2on2;
extern    double   __canuto_mod_MOD_dri;
extern    double   __canuto_mod_MOD_deltheta_r;
extern    double   __canuto_mod_MOD_b1;
extern    double   __canuto_mod_MOD_theta_rcrp;
extern    double   __canuto_mod_MOD_theta_rcrn;
extern    double   __canuto_mod_MOD_ako;
extern    double   __canuto_mod_MOD_back_l_0;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

constexpr int      IFCHENGCON = 1;
constexpr int      IDEFMLD    = 0;
constexpr double   DELTEMLD   = 0.5e0;
constexpr double   DELRHMLD   = 0.125e-3;
constexpr int      ILOMEGA    = 0;
constexpr double   AMLDMINLOM = 5.e4;
constexpr int      NBIG       = 100;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern    double   canuto_mod_mp_rib_[2*MT+1];
extern    double   canuto_mod_mp_ridb_[2*MT+1];
extern    double   canuto_mod_mp_slq2b_[2*MT+1][2*MT+1];
extern    double   canuto_mod_mp_smb_[2*MT+1][2*MT+1];
extern    double   canuto_mod_mp_shb_[2*MT+1][2*MT+1];
extern    double   canuto_mod_mp_ssb_[2*MT+1][2*MT+1];
                   
extern    int      canuto_mod_mp_irimax_[2*MT+1];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern    double   __canuto_mod_MOD_rib[2*MT+1];
extern    double   __canuto_mod_MOD_ridb[2*MT+1];
extern    double   __canuto_mod_MOD_slq2b[2*MT+1][2*MT+1];
extern    double   __canuto_mod_MOD_smb[2*MT+1][2*MT+1];
extern    double   __canuto_mod_MOD_shb[2*MT+1][2*MT+1];
extern    double   __canuto_mod_MOD_ssb[2*MT+1][2*MT+1];
                   
extern    int      __canuto_mod_MOD_irimax[2*MT+1];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

constexpr int      MT_RA_R        = NPOSAPPROX - 1;
constexpr int      N_THETA_R_OCT  = 
    static_cast<int>((((PI_CONST / 4.e0) * MT_RA_R) / 15.e0)) * 15;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern    double   canuto_mod_mp_sisamax_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_ra_rmax_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_c_y_r0_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_back_ra_r_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_sm_r1_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_sh_r1_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_ss_r1_[4*N_THETA_R_OCT+1];
extern    double   canuto_mod_mp_slq2_r1_[4*N_THETA_R_OCT+1];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern    double   __canuto_mod_MOD_sisamax[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_ra_rmax[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_c_y_r0[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_back_ra_r[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_sm_r1[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_sh_r1[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_ss_r1[4*N_THETA_R_OCT+1];
extern    double   __canuto_mod_MOD_slq2_r1[4*N_THETA_R_OCT+1];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

constexpr int      IFPOLARTABLEWRITE = 0;
constexpr int      IFBG_THETA_INTERP = 1;

constexpr double   BACK_PH_0         =   
    (6.e-5) * ((1.e2) / (2.e0 * PI_CONST));

constexpr double   ADJUST_GARGETT    = 1.e0;

constexpr double   BACK_K_0          =   
    (0.1e0) * (2.e0) * PI_CONST * (1.e-2) * ADJUST_GARGETT;

constexpr double   BACK_DEL_0        = PI_CONST / BACK_K_0;
constexpr double   BACK_S2           = BACK_PH_0 * BACK_K_0;
constexpr double   BACK_SM2          = 1.e0 / BACK_S2;
constexpr double   V_BACK0           = 0.01;
constexpr double   T_BACK0           = 0.01;
constexpr double   S_BACK0           = 0.01;
constexpr double   RI_INTERNAL       = 1.e0;
constexpr double   BACKFRAC          = 85.e-2;
constexpr double   BACKFACT          = 1 / E;



#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern    double   canuto_mod_mp_ria_[NTBL];
extern    double   canuto_mod_mp_slq2a_[NTBL];
extern    double   canuto_mod_mp_sma_[NTBL];
extern    double   canuto_mod_mp_sha_[NTBL];

extern int canuto_mod_mp_ifirst_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern    double   __canuto_mod_MOD_ria[NTBL];
extern    double   __canuto_mod_MOD_slq2a[NTBL];
extern    double   __canuto_mod_MOD_sma[NTBL];
extern    double   __canuto_mod_MOD_sha[NTBL];

extern int __canuto_mod_MOD_ifirst;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_CANUTO_MOD_H_
