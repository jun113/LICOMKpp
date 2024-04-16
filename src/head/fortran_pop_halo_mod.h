#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_HALO_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_HALO_MOD_H_

#include "def-undef.h"

#include "cpp_pop_blocks_mod.h"

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int     pop_halomod_mp_bufsizesend_;
extern int     pop_halomod_mp_bufsizerecv_;
              
extern int*    pop_halomod_mp_bufsendi4_;
extern int*    pop_halomod_mp_bufrecvi4_;

extern float*  pop_halomod_mp_bufsendr4_;
extern float*  pop_halomod_mp_bufrecvr4_;

extern double* pop_halomod_mp_bufsendr8_;
extern double* pop_halomod_mp_bufrecvr8_;

extern int*    pop_halomod_mp_buftripolei4_;
extern float*  pop_halomod_mp_buftripoler4_;
extern double* pop_halomod_mp_buftripoler8_;

extern int*    pop_halomod_mp_pophalo_bufrecvtask_;
extern int*    pop_halomod_mp_pophalo_bufsendtask_;
extern int*    pop_halomod_mp_pophalo_bufsizesend_;
extern int*    pop_halomod_mp_pophalo_bufsizerecv_;
extern int*    pop_halomod_mp_pophalo_bufsrclocaladdr_;
extern int*    pop_halomod_mp_pophalo_bufdstlocaladdr_;
extern int*    pop_halomod_mp_pophalo_bufsendaddr_;
extern int*    pop_halomod_mp_pophalo_bufrecvaddr_;

// extern int*    pop_halomod_mp_pop_halomsgcreate_iglobal_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int     __pop_halomod_MOD_bufsizesend;
extern int     __pop_halomod_MOD_bufsizerecv;
              
extern int*    __pop_halomod_MOD_bufsendi4;
extern int*    __pop_halomod_MOD_bufrecvi4;

extern float*  __pop_halomod_MOD_bufsendr4;
extern float*  __pop_halomod_MOD_bufrecvr4;

extern double* __pop_halomod_MOD_bufsendr8;
extern double* __pop_halomod_MOD_bufrecvr8;

extern int*    __pop_halomod_MOD_buftripolei4;
extern float*  __pop_halomod_MOD_buftripoler4;
extern double* __pop_halomod_MOD_buftripoler8;

extern int*    __pop_halomod_MOD_pophalo_bufrecvtask;
extern int*    __pop_halomod_MOD_pophalo_bufsendtask;
extern int*    __pop_halomod_MOD_pophalo_bufsizesend;
extern int*    __pop_halomod_MOD_pophalo_bufsizerecv;
extern int*    __pop_halomod_MOD_pophalo_bufsrclocaladdr;
extern int*    __pop_halomod_MOD_pophalo_bufdstlocaladdr;
extern int*    __pop_halomod_MOD_pophalo_bufsendaddr;
extern int*    __pop_halomod_MOD_pophalo_bufrecvaddr;

// extern int*    __pop_halomod_MOD_pop_halomsgcreate_iglobal;
#endif // LICOM_ENABLE_FORTRAN_COMODILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_HALO_MOD_H_
