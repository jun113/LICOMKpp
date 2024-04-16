#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_CANUTO_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_CANUTO_MOD_H_

#include "fortran_canuto_mod.h"
namespace CppCanutoMod {

extern double  &eps_bot__under;   
extern double (&and2on2a1)[2 * MT + 1];
extern double (&amtaun2a1)[2 * MT + 1];

extern double  &dand2on2;

extern double (&sma1)[2 * MT + 1];
extern double (&sha1)[2 * MT + 1];
extern double (&ssa1)[2 * MT + 1];

extern double  &rri;
extern double  &rnd2on2;
extern double  &dri;
extern double  &deltheta_r;
extern double  &b1;
extern double  &theta_rcrp;
extern double  &theta_rcrn;
extern double  &ako;
extern double  &back_l_0;

extern double (&rib)[2 * MT + 1]; 
extern double (&ridb)[2 * MT + 1]; 
extern double (&slq2b)[2 * MT + 1][2 * MT + 1]; 
extern double (&smb)[2 * MT + 1][2 * MT + 1]; 
extern double (&shb)[2 * MT + 1][2 * MT + 1]; 
extern double (&ssb)[2 * MT + 1][2 * MT + 1];
extern int    (&irimax)[2 * MT + 1];

extern double (&sisamax)[4 * N_THETA_R_OCT + 1]; 
extern double (&ra_rmax)[4 * N_THETA_R_OCT + 1]; 
extern double (&c_y_r0)[4 * N_THETA_R_OCT + 1]; 
extern double (&back_ra_r)[4 * N_THETA_R_OCT + 1];
extern double (&sm_r1)[4 * N_THETA_R_OCT + 1]; 
extern double (&sh_r1)[4 * N_THETA_R_OCT + 1]; 
extern double (&ss_r1)[4 * N_THETA_R_OCT + 1]; 
extern double (&slq2_r1)[4 * N_THETA_R_OCT + 1];

extern double (&ria)[NTBL]; 
extern double (&slq2a)[NTBL]; 
extern double (&sma)[NTBL]; 
extern double (&sha)[NTBL];
extern int    &ifirst;
} // namespace CppCanutoMod



#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_CANUTO_MOD_H_
