#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PARAM_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PARAM_MOD_H_
// i j k m n
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int param_mod_mp_ierr_;
extern int param_mod_mp_mytid_;
extern int param_mod_mp_jj_start_;
extern int param_mod_mp_jj_end_;

extern int param_mod_mp_my_task_;
extern int param_mod_mp_master_task_;

extern double param_mod_mp_am_hor_;
extern double param_mod_mp_ah_hor_;
extern double param_mod_mp_am_bihar_;
extern double param_mod_mp_ah_bihar_;

extern int param_mod_mp_num_cpl_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __param_mod_MOD_ierr;
extern int __param_mod_MOD_mytid;
extern int __param_mod_MOD_jj_start;
extern int __param_mod_MOD_jj_end;

extern int __param_mod_MOD_my_task;
extern int __param_mod_MOD_master_task;

extern double __param_mod_MOD_am_hor;
extern double __param_mod_MOD_ah_hor;
extern double __param_mod_MOD_am_bihar;
extern double __param_mod_MOD_ah_bihar;

extern int __param_mod_MOD_num_cpl;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PARAM_MOD_H_
