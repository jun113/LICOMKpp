#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_EXTERN_FUNCTION_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_EXTERN_FUNCTION_H_

#include "def-undef.h"

extern "C" void readyt_();
extern "C" void readyc_();
extern "C" void barotr_();
extern "C" void bclinc_();
extern "C" void tracer_();
extern "C" void icesnow_();
extern "C" void convadj_();
extern "C" void energy_();
extern "C" void accumm_();

extern "C" void allocate_readyt_();
extern "C" void allocate_readyc_();
extern "C" void allocate_tracer_();

extern "C" void get_block_info_();
extern "C" void get_f_pop_distrb_info_(int *, int *, int *);

extern "C" void pop_halo_create_from_fortran_get_(int*, int*, int*, int*,
    double*, double*, double*, double*, double*, double*, double*, double*);

extern "C" void pop_halo_create_from_fortran_get_1_(int*, int*, int*, int*);
extern "C" void pop_halo_create_from_fortran_get_2_(int*, int*, int*, 
    int*, int*, int*, int*, int*);
// extern "C" void pop_halo_create_from_fortran_get_(int*, int*, int*, int*);
// extern "C" void pop_halo_create_from_fortran_free_buf_();

extern "C" void scatter_global_jra_r4_(float*, float*, int*);
extern "C" void scatter_global_jra_r8_(double*, double*, int*);

extern "C" void jra_daily_();
extern "C" void addps_();
extern "C" void accumm_();
extern "C" void ssaveins_();
extern "C" void nextstep_();

extern "C" void turb_2_(
    double *, // z
    double *, // t
    double *, // s
    double *, // rh
    double *, // ri
    double *, // rid
    double *, // s2
    double *, // fricmx 
    double *, // wndmix 
    double *, // v_back 
    double *, // t_back 
    double *, // s_back 
    double *, // an2 
    double *, // ustar_ 
    double *, // buoytur
    double *, // buoysol
    double *, // coriol
    double *, // amld
    double *, // akm
    double *, // akh
    double *, // aks
    int    *, // n
    int    *, // na
    int    *, // nmax
    int    *, // isurfuse
    int    *, // ifextermld
    int    *, // ifoutput
    int    *, // ii
    int    *, // jj
    int    *  /* iblock*/);

extern "C" void fortran_mpi_barrier_();

extern "C" void pop_haloupdate_readyc_(int *);

extern "C" void pop_haloupdate_barotr1_(int *);
extern "C" void pop_haloupdate_barotr2_(int *, int *);

extern "C" void pop_haloupdate_bclinc2_(int *);
extern "C" void pop_haloupdate_bclinc3_(int *);
extern "C" void pop_haloupdate_bclinc33_(int *);

extern "C" void pop_haloupdate_tracer1_(int *);
extern "C" void pop_haloupdate_tracer2_(int *);

extern "C" void pop_haloupdate_smts_(int *);

extern "C" void mpi_bcast_tracer_(double *);
extern "C" void mpi_reduce_tracer_(double *, double *);

extern "C" void free_buf_mod_();
#ifdef CANUTO
extern "C" void free_canuto_mod_();
#endif // CANUTO
extern "C" void free_dyn_mod_();
extern "C" void free_forc_mod_();
extern "C" void free_grid_();
extern "C" void free_forc_mod_();
#ifndef BIHAR
extern "C" void free_hmix_del2_();
#else // BIHAR
extern "C" void free_hmix_del4_();
#endif // BIHAR
#ifdef ISO
extern "C" void free_isopyc_mod_();
#endif // ISO

extern "C" void free_output_mod_();
extern "C" void free_pconst_mod_();
extern "C" void free_pmix_mod_();
extern "C" void free_work_mod_();
extern "C" void free_tracer_mod_();
extern "C" void free_tmp_var_();

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_EXTERN_FUNCTION_H_
