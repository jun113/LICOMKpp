#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_BLOCKS_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_BLOCKS_H_

#include "fortran_blocks.h"
namespace CppBlocks {

using func_ptr_get_block = struct block (*) (int*, int*);

extern func_ptr_get_block get_block;

extern block* &all_blocks;

extern int &ib;
extern int &ie;
extern int &jb;
extern int &je;

} // namespace CppBlocks

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_BLOCKS_H_
