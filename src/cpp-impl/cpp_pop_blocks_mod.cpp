#include "def-undef.h"

#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pop_blocks_mod.h"
#include "../head/cpp_extern_functions.h"
#include "../head/fortran_pop_blocks_mod.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
namespace CppPOPBlocksMod {

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &POP_numBlocks  = pop_blocksmod_mp_pop_numblocks_;
int &POP_numBlocksX = pop_blocksmod_mp_pop_numblocksx_;
int &POP_numBlocksY = pop_blocksmod_mp_pop_numblocksy_;

int* &POP_HaloMsgCreate_iGlobal = pop_blocksmod_mp_pop_halomsgcreate_iglobal_;

void(&get_allblocks_ij)(int *, int *, int *) = pop_blocksmod_mp_get_allblocks_ij_;

void(&pop_blocks_get_block_info1)(int *, int *, int *, int *, int *) = pop_blocksmod_mp_pop_blocksgetblockinfo1_;
void(&pop_blocks_get_block_info2)(int *, int *, int *, int *, int *) = pop_blocksmod_mp_pop_blocksgetblockinfo2_;
void(&pop_blocks_get_block_info3)(int *, int *, int *)               = pop_blocksmod_mp_pop_blocksgetblockinfo3_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &POP_numBlocks  = __pop_blocksmod_MOD_pop_numblocks;
int &POP_numBlocksX = __pop_blocksmod_MOD_pop_numblocksx;
int &POP_numBlocksY = __pop_blocksmod_MOD_pop_numblocksy;
int* &POP_HaloMsgCreate_iGlobal = __pop_blocksmod_MOD_pop_halomsgcreate_iglobal;

void(&get_allblocks_ij)(int *, int *, int *) = __pop_blocksmod_MOD_get_allblocks_ij;

void(&pop_blocks_get_block_info1)(int *, int *, int *, int *, int *) = __pop_blocksmod_MOD_pop_blocksgetblockinfo2;
void(&pop_blocks_get_block_info2)(int *, int *, int *, int *, int *) = __pop_blocksmod_MOD_pop_blocksgetblockinfo2;
void(&pop_blocks_get_block_info3)(int *, int *, int *)               = __pop_blocksmod_MOD_pop_blocksgetblockinfo3;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

// TODO
// void pop_blocks_get_block_info (const int &blockId,
// 		int &ib, int &ie, int &jb, int &je) {
// 	if (blockId < 1 || blockId > POP_numBlocks) {
// 		printf ("POP_BlocksGetBlockInfo: invalid blockID\n");
// 		return ;
// 	}
	
// 	return ;
// }
void pop_blocks_get_nbr_id (const int &blockId,
		const int &direction, const char *iBoundary,
				const char *jBoundary, int &nbrID) {
	// DESCRIPTION:
	// This function returns the block id of a neighboring block in a
	// requested direction.  Supported directions currently include:
	//     POP\_blocksNorth             (i  ,j-1)
	//     POP\_blocksSouth             (i  ,j+1)
	//     POP\_blocksEast              (i+1,j  )
	//     POP\_blocksWest              (i-1,j  )
	//     POP\_blocksNorthEast         (i+1,j-1)
	//     POP\_blocksNorthWest         (i-1,j-1)
	//     POP\_blocksSouthEast         (i  ,j+1)
	//     POP\_blocksSouthWest         (i-1,j+1)
	//     POP\_blocksNorth2            (i  ,j-2)
	//     POP\_blocksSouth2            (i  ,j+2)
	//     POP\_blocksEast2             (i+2,j  )
	//     POP\_blocksWest2             (i-2,j  )
	//     POP\_blocksNorthEast2        (i+2,j-2)
	//     POP\_blocksNorthWest2        (i-2,j-2)
	//     POP\_blocksSouthEast2        (i+2,j+2)
	//     POP\_blocksSouthWest2        (i-2,j+2)
	//     POP\_blocksEastNorthEast     (i+2,j-1)
	//     POP\_blocksEastSouthEast     (i+2,j+1)
	//     POP\_blocksWestNorthWest     (i-2,j-1)
	//     POP\_blocksWestSouthWest     (i-2,j+1)
	//     POP\_blocksNorthNorthEast    (i+1,j-2)
	//     POP\_blocksSouthSouthEast    (i+1,j+2)
	//     POP\_blocksNorthNorthWest    (i-1,j-2)
	//     POP\_blocksSouthSouthWest    (i-1,j+2)

	// i,j block location of current block
	int iBlock, jBlock;
	// ------------------------------------------------
	// retrieve info for current block
	// ------------------------------------------------
	int blkId = blockId;
	pop_blocks_get_block_info3(&blkId, &iBlock, &jBlock);

	// ---------------------------------------------------
	// compute i,j block location of neighbor
	// ---------------------------------------------------
	// i,j block location of neighboring block
	int inbr, jnbr;
	if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH) {
		inbr = iBlock;
		jnbr = jBlock - 1;
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = 1;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 1;
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH) {
		inbr = iBlock;
		jnbr = jBlock + 1;
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_EAST) {
		inbr = iBlock + 1;
		jnbr = jBlock;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_WEST) {
		inbr = iBlock - 1;
		jnbr = jBlock;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_EAST) {
		inbr = iBlock + 1;
		jnbr = jBlock - 1;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = 1;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock;
				if (inbr == 0) {
					inbr = POP_numBlocksX;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_WEST) {
		inbr = iBlock - 1;
		jnbr = jBlock - 1;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = 1;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 2;
				if (inbr > POP_numBlocksX) {
					inbr = 1;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_EAST) {
		inbr = iBlock + 1;
		jnbr = jBlock + 1;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_WEST) {
		inbr = iBlock - 1;
		jnbr = jBlock + 1;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH2) {
		inbr = iBlock;
		jnbr = jBlock - 2;
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr - POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 1;
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH2) {
		inbr = iBlock;
		jnbr = jBlock + 2;
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_EAST2) {
		inbr = iBlock + 2;
		jnbr = jBlock;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_WEST2) {
		inbr = iBlock - 2;
		jnbr = jBlock;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_EAST2) {
		inbr = iBlock + 2;
		jnbr = jBlock - 2;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr - POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock - 1;
				if (inbr == 0) {
					inbr = POP_numBlocksX;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_WEST2) {
		inbr = iBlock - 2;
		jnbr = jBlock - 2;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = 1;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 3;
				if (inbr > POP_numBlocksX) {
					inbr = 0;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_EAST2) {
		inbr = iBlock + 2;
		jnbr = jBlock + 2;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_WEST2) {
		inbr = iBlock - 2;
		jnbr = jBlock + 2;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_EAST_NORTH_EAST) {
		inbr = iBlock + 2;
		jnbr = jBlock - 1;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr - POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock - 1;
				if (inbr == 0) {
					inbr = POP_numBlocksX;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_WEST_NORTH_WEST) {
		inbr = iBlock - 2;
		jnbr = jBlock - 1;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr + POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 3;
				if (inbr > POP_numBlocksX) {
					inbr = 0;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_EAST_SOUTH_EAST) {
		inbr = iBlock + 2;
		jnbr = jBlock + 1;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_WEST_SOUTH_WEST) {
		inbr = iBlock - 2;
		jnbr = jBlock + 1;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_NORTH_EAST) {
		inbr = iBlock + 1;
		jnbr = jBlock - 2;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr - POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock;
				if (inbr == 0) {
					inbr = POP_numBlocksX;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_NORTH_NORTH_WEST) {
		inbr = iBlock - 1;
		jnbr = jBlock - 2;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr < 1) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = jnbr - POP_numBlocksY;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
       // return negative j value to flag tripole
       // i index of main northern neighbor across the
       // tripole cut - may also need i+1,i-1 to get
       // other points if there has been padding or
       // if the block size does not divide the domain
       // evenly
				inbr = POP_numBlocksX - iBlock + 2;
				if (inbr > POP_numBlocksX) {
					inbr = 0;
				}
				jnbr = jnbr - 1;
			} else {
				printf ("POP_BlocksGetNbrID: unknown north boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_SOUTH_EAST) {
		inbr = iBlock + 1;
		jnbr = jBlock + 2;
		if (inbr > POP_numBlocksX) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = inbr - POP_numBlocksX;
			} else {
				printf ("POP_BlocksGetNbrID: unknown east boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else if (direction == CppPOPBlocksMod::POP_BLOCKS_SOUTH_SOUTH_WEST) {
		inbr = iBlock - 1;
		jnbr = jBlock + 2;
		if (inbr < 1) {
			if (str_trim_cmp(iBoundary, "closed") == 0) {
				inbr = 0;
			} else if (str_trim_cmp(iBoundary, "cyclic") == 0) {
				inbr = POP_numBlocksX + inbr;
			} else {
				printf ("POP_BlocksGetNbrID: unknown west boundary\n");
				exit(0);
			}
		}
		if (jnbr > POP_numBlocksY) {
			if (str_trim_cmp(jBoundary, "closed") == 0) {
				jnbr = 0;
			} else if (str_trim_cmp(jBoundary, "cyclic") == 0) {
				jnbr = POP_numBlocksY + jnbr;
			} else if (str_trim_cmp(jBoundary, "tripole") == 0) {
				jnbr = 0;
			} else {
				printf ("POP_BlocksGetNbrID: unknown south boundary\n");
				exit(0);
			}
		}
	} else {
		printf ("POP_BlocksGetNbrID: unknown direction\n");
		exit(0);
	}
	// ---------------------------------------------------------
	// now get block id for this neighbor block
	// ---------------------------------------------------------
	if (inbr > 0 && jnbr > 0) {
		get_allblocks_ij(&inbr, &jnbr, &nbrID);
	} else if (inbr > 0 && jnbr < 0) {
		// tripole upper boundary return negative value to flag tripole
		int abs_jnbr = abs(jnbr);
		get_allblocks_ij(&inbr, &abs_jnbr, &nbrID);
		nbrID = - nbrID;
	} else {
		// neighbor outside domain
		nbrID = 0;
	}
  return ;
}

} // namespace CppPOPBlocksMod

#endif // LICOM_ENABLE_FORTRAN