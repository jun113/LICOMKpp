#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_POP_VLOCKS_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_POP_VLOCKS_MOD_H_

#include "def-undef.h"

#ifndef LICOM_ENABLE_FORTRAN
#include "cpp_param_mod.h"
#include "fortran_pop_blocks_mod.h"
namespace CppPOPBlocksMod {

using CppParamMod::LICOM_BLOCKSIZEX;
using CppParamMod::LICOM_BLOCKSIZEY;

constexpr int POP_HALO_WIDTH = 2;

constexpr int POP_NX_BLOCK = LICOM_BLOCKSIZEX + 2 * POP_HALO_WIDTH;
constexpr int POP_NY_BLOCK = LICOM_BLOCKSIZEY + 2 * POP_HALO_WIDTH;

constexpr int POP_BLOCKS_NORTH            = 1;
constexpr int POP_BLOCKS_SOUTH            = 2;
constexpr int POP_BLOCKS_EAST             = 3;
constexpr int POP_BLOCKS_WEST             = 4;
constexpr int POP_BLOCKS_NORTH_EAST       = 5;
constexpr int POP_BLOCKS_NORTH_WEST       = 6;
constexpr int POP_BLOCKS_SOUTH_EAST       = 7;
constexpr int POP_BLOCKS_SOUTH_WEST       = 8;
constexpr int POP_BLOCKS_NORTH2           = 9;
constexpr int POP_BLOCKS_SOUTH2           = 10;
constexpr int POP_BLOCKS_EAST2            = 11;
constexpr int POP_BLOCKS_WEST2            = 12;
constexpr int POP_BLOCKS_NORTH_EAST2      = 13;
constexpr int POP_BLOCKS_NORTH_WEST2      = 14;
constexpr int POP_BLOCKS_SOUTH_EAST2      = 15;
constexpr int POP_BLOCKS_SOUTH_WEST2      = 16;
constexpr int POP_BLOCKS_EAST_NORTH_EAST  = 17;
constexpr int POP_BLOCKS_EAST_SOUTH_EAST  = 18;
constexpr int POP_BLOCKS_WEST_NORTH_WEST  = 19;
constexpr int POP_BLOCKS_WEST_SOUTH_WEST  = 20;
constexpr int POP_BLOCKS_NORTH_NORTH_EAST = 21;
constexpr int POP_BLOCKS_SOUTH_SOUTH_EAST = 22;
constexpr int POP_BLOCKS_NORTH_NORTH_WEST = 23;
constexpr int POP_BLOCKS_SOUTH_SOUTH_WEST = 24;

extern int &POP_numBlocks;
extern int &POP_numBlocksX;
extern int &POP_numBlocksY;

extern int* &POP_HaloMsgCreate_iGlobal;

extern void pop_blocks_get_nbr_id (const int &, const int &, 
		const char *, const char *, int &);

extern void(&pop_blocks_get_block_info1)(int *, int *, int *, int *, int *);
extern void(&pop_blocks_get_block_info2)(int *, int *, int *, int *, int *);
extern void(&pop_blocks_get_block_info3)(int *, int *, int *);
} // namespace CppPOPBlocksMod

#endif // LICOM_ENABLE_FORTRAN
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_POP_VLOCKS_MOD_H_
