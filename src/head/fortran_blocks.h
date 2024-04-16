#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BLOCKS_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BLOCKS_H_

struct block {
  int block_id;
  int local_id;
  int ib, ie, jb, je;
  int iblock, jblock; 
  //int i_glob[18], j_glob[18];
  int* i_glob, j_glob;
};

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern "C" struct block blocks_mp_get_block_(int*, int*);
extern block *blocks_mp_all_blocks_;

extern int blocks_mp_ib_;
extern int blocks_mp_ie_;
extern int blocks_mp_jb_;
extern int blocks_mp_je_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern "C" struct block __blocks_MOD_get_block(int*, int*);
extern block *__blocks_MOD_all_blocks;

extern int __blocks_MOD_ib;
extern int __blocks_MOD_ie;
extern int __blocks_MOD_jb;
extern int __blocks_MOD_je;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BLOCKS_H_
