#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/fortran_blocks.h"
namespace CppBlocks {

using func_ptr_get_block = struct block (*) (int*, int*);

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
func_ptr_get_block get_block = blocks_mp_get_block_;

block* &all_blocks = blocks_mp_all_blocks_;

int &ib = blocks_mp_ib_;
int &ie = blocks_mp_ie_;
int &jb = blocks_mp_jb_;
int &je = blocks_mp_je_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU

func_ptr_get_block get_block = __blocks_MOD_get_block;

block* &all_blocks = __blocks_MOD_all_blocks;

int &ib = __blocks_MOD_ib;
int &ie = __blocks_MOD_ie;
int &jb = __blocks_MOD_jb;
int &je = __blocks_MOD_je;

#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppBlocks
#endif // LICOM_ENABLE_FORTRAN
