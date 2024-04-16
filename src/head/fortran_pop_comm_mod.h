#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_COMM_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_COMM_MOD_H_
#include "def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int pop_commmod_mp_pop_communicator_;
extern int pop_commmod_mp_pop_mytask_;
extern int pop_commmod_mp_pop_mastertask_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __pop_commmod_MOD_pop_communicator;
extern int __pop_commmod_MOD_pop_mytask;
extern int __pop_commmod_MOD_pop_mastertask;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM_ENABLE_FORTRAN

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_COMM_MOD_H_