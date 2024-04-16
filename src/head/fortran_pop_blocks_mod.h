#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_BLOCKS_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_BLOCKS_MOD_H_

#include "def-undef.h"

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int pop_blocksmod_mp_pop_numblocks_;
extern int pop_blocksmod_mp_pop_numblocksx_;
extern int pop_blocksmod_mp_pop_numblocksy_;

extern int* pop_blocksmod_mp_pop_halomsgcreate_iglobal_;

extern "C" void pop_blocksmod_mp_get_allblocks_ij_(int *, int *, int *);

extern "C" void pop_blocksmod_mp_pop_blocksgetblockinfo1_(int*, int*, int*, int*, int*);
extern "C" void pop_blocksmod_mp_pop_blocksgetblockinfo2_(int*, int*, int*, int*, int*);
extern "C" void pop_blocksmod_mp_pop_blocksgetblockinfo3_(int*, int*, int*);
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __pop_blocksmod_MOD_pop_numblocks;
extern int __pop_blocksmod_MOD_pop_numblocksx;
extern int __pop_blocksmod_MOD_pop_numblocksy;

extern int* __pop_blocksmod_MOD_pop_halomsgcreate_iglobal;

extern "C" void __pop_blocksmod_MOD_get_allblocks_ij(int *, int *, int *);

extern "C" void __pop_blocksmod_MOD_pop_blocksgetblockinfo1(int*, int*, int*, int*, int*);
extern "C" void __pop_blocksmod_MOD_pop_blocksgetblockinfo2(int*, int*, int*, int*, int*);
extern "C" void __pop_blocksmod_MOD_pop_blocksgetblockinfo3(int*, int*, int*);
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_POP_BLOCKS_MOD_H_
