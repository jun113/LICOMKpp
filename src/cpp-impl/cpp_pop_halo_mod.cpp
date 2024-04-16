#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "cpp_pop_reductions_mod.hpp"

#include "../head/cpp_domain.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pop_comm_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_blocks_mod.h"
#include "../head/cpp_pop_comm_mod.h"
#include "../head/fortran_pop_halo_mod.h"
#include "../head/cpp_pop_grid_horz_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <cstdio>

#include <iostream>
using namespace::std;
// memory footprint

namespace CppPOPHaloMod {

double* arrCommPriorK = new double[CppParamMod::KM 
		* CppParamMod::JMT * CppParamMod::IMT];

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &bufSizeSend        = pop_halomod_mp_bufsizesend_;
int &bufSizeRecv        = pop_halomod_mp_bufsizesend_;

// int*    (&bufSendI4)    = pop_halomod_mp_bufsendi4_;
// int*    (&bufRecvI4)    = pop_halomod_mp_bufrecvi4_;
                     
// float*  (&bufSendR4)    = pop_halomod_mp_bufsendr4_;
// float*  (&bufRecvR4)    = pop_halomod_mp_bufrecvr4_;
                     
// double* (&bufSendR8)    = pop_halomod_mp_bufsendr8_;
// double* (&bufRecvR8)    = pop_halomod_mp_bufrecvr8_;

// int*    (&bufTripoleI4) = pop_halomod_mp_buftripolei4_;
// float*  (&bufTripoleR4) = pop_halomod_mp_buftripoler4_;
// double* (&bufTripoleR8) = pop_halomod_mp_buftripoler8_;
// int* (&POPHalo_bufRecvTask)     = pop_halomod_mp_pophalo_bufrecvtask_;
// int* (&POPHalo_bufSendTask)     = pop_halomod_mp_pophalo_bufsendtask_;
// int* (&POPHalo_bufSizeSend)     = pop_halomod_mp_pophalo_bufsizesend_;
// int* (&POPHalo_bufSizeRecv)     = pop_halomod_mp_pophalo_bufsizerecv_;
// int* (&POPHalo_bufSrcLocalAddr) = pop_halomod_mp_pophalo_bufsrclocaladdr_;
// int* (&POPHalo_bufDstLocalAddr) = pop_halomod_mp_pophalo_bufdstlocaladdr_;
// int* (&POPHalo_bufRecvAddr)     = pop_halomod_mp_pophalo_bufrecvaddr_;
// int* (&POPHalo_bufSendAddr)     = pop_halomod_mp_pophalo_bufsendaddr_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int     &bufSizeSend    = __pop_halomod_MOD_bufsizesend;
int     &bufSizeRecv    = __pop_halomod_MOD_bufsizesend;

// int*    &bufSendI4      = __pop_halomod_MOD_bufsendi4;
// int*    &bufRecvI4      = __pop_halomod_MOD_bufrecvi4;
                        
// float*  &bufSendR4      = __pop_halomod_MOD_bufsendr4;
// float*  &bufRecvR4      = __pop_halomod_MOD_bufrecvr4;
                        
// double* &bufSendR8      = __pop_halomod_MOD_bufsendr8;
// double* &bufRecvR8      = __pop_halomod_MOD_bufrecvr8;

// int*    &bufTripoleI4   = __pop_halomod_MOD_buftripolei4;
// float*  &bufTripoleR4   = __pop_halomod_MOD_buftripoler4;
// double* &bufTripoleR8   = __pop_halomod_MOD_buftripoler8;

// int* (&POPHalo_bufRecvTask)     = __pop_halomod_MOD_pophalo_bufrecvtask;
// int* (&POPHalo_bufSendTask)     = __pop_halomod_MOD_pophalo_bufsendtask;
// int* (&POPHalo_bufSizeSend)     = __pop_halomod_MOD_pophalo_bufsizesend;
// int* (&POPHalo_bufSizeRecv)     = __pop_halomod_MOD_pophalo_bufsizerecv;
// int* (&POPHalo_bufSrcLocalAddr) = __pop_halomod_MOD_pophalo_bufsrclocaladdr;
// int* (&POPHalo_bufDstLocalAddr) = __pop_halomod_MOD_pophalo_bufdstlocaladdr;
// int* (&POPHalo_bufRecvAddr)     = __pop_halomod_MOD_pophalo_bufrecvaddr;
// int* (&POPHalo_bufSendAddr)     = __pop_halomod_MOD_pophalo_bufsendaddr;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

// buffer for use to send max in 3D r8 halo updates
// [numMsgSend][bufSizeSend * KM]
void* bufSend = nullptr;
// buffer for use to recv max in 3D r8 halo updates
// [numMsgRecv][bufSizeRecv * KM]
void* bufRecv = nullptr;
// [KM][POP_haloWidth+1][nxGlobal]
void* bufTripole = nullptr;

// MPI request ids
MPI_Request* sndRequest = nullptr;
MPI_Request* rcvRequest = nullptr;
// MPI status flags
MPI_Status* sndStatus = nullptr;
MPI_Status* rcvStatus = nullptr;

// POP_HaloIncrementMsgCount
void pop_halo_increment_msg_count(
		int* const sndCounter, int* const rcvCounter,
		const int &srcProc, const int &dstProc,
		const int &msgSize, const int &numProcs) {
	// DESCRIPTION:
	// This is a utility routine to increment the arrays for counting
	// whether messages are required.  It checks the source and destination
	// task to see whether the current task needs to send, receive or
	// copy messages to fill halo regions (ghost cells).
	// INPUT:
	// srcProc: source processor for communication
	// dstProc: destination processor for communication
	// msgSize: number of words for this message
	// numProcs: total of proces
	// OUTPUT:
	// sndCounter: array for counting messages to be sent
	// rcvCounter: array for counting messages to be sent
	using CppPOPCommMod::POP_myTask;
	if (srcProc < 0 || dstProc < 0 || 
			srcProc > numProcs || dstProc > numProcs) {
		printf ("POP_HaloIncrementMsgCount: invalid processor number\n");
		exit(0);
	}
	// ------------------------------------
	// if destination all land or outside closed boundary (dstProc = 0), 
	// then no send is necessary, so do the rest only for dstProc /= 0
	// ------------------------------------
	if (dstProc == 0) return;
	// ------------------------------------
	// if the current processor is the source, must send data 
	// local copy if dstProc = srcProc
	// ------------------------------------
	if (srcProc == POP_myTask + 1) {
		// TODO
		sndCounter[dstProc - 1] = sndCounter[dstProc - 1] + msgSize; 
	}
	// ------------------------------------
	// if the current processor is the destination, must receive data 
	// local copy if dstProc = srcProc
	// ------------------------------------
	if (dstProc == POP_myTask + 1) {
		if (srcProc > 0) {
      // the source block has ocean points
      // count as a receive from srcProc
			rcvCounter[srcProc - 1] = rcvCounter[srcProc - 1] + msgSize;
		} else {
      // if the source block has been dropped, create
      // a local copy to fill halo with a fill value
			rcvCounter[dstProc - 1] = rcvCounter[dstProc - 1] + msgSize;
		}
	}
	return ;
}

// POP_HaloMsgCreate
void pop_halo_msg_create(CppPOPHaloMod::pop_halo &halo,
		const int &srcBlock, const int &srcProc, const int &srcLocalID,
		const int &dstBlock, const int &dstProc, const int &dstLocalID,
				const int &nxGlobal, const char *direction) {

	using CppPOPCommMod::POP_myTask;
	using CppPOPHaloMod::bufSizeSend;
	using CppPOPHaloMod::bufSizeRecv;
	using CppPOPBlocksMod::POP_HALO_WIDTH;
	
  int msgIndx;               			// message counter and index into msg array
  // int blockIndx;             			// block counter and index into msg array
  int bufSize;               			// size of message buffer
  int ibSrc, ieSrc, jbSrc, jeSrc; // phys domain info for source block
	int ibDst, ieDst, jbDst, jeDst; // phys domain info for dest   block
  int i,j;                  		  // dummy loop index

  // global i index for location in tripole
	int *iGlobal = nullptr;

	// -----------------------------------------------------
	// if destination all land or outside closed boundary (dstProc = 0), 
	// then no send is necessary, so do the rest only for dstProc /= 0
		if (dstProc == 0) return;
	// -----------------------------------------------------

	// get block information if either block is local
	if (srcProc == (POP_myTask + 1) || dstProc == (POP_myTask + 1)) {

		if (srcBlock >= 0 && dstBlock >= 0) {
			int srcBlk = srcBlock;
			CppPOPBlocksMod::pop_blocks_get_block_info1(
					&srcBlk, &ibSrc, &ieSrc, &jbSrc, &jeSrc);
		} else { 
			// tripole - need iGlobal info
			int abs_srcBlock = abs(srcBlock);
			CppPOPBlocksMod::pop_blocks_get_block_info2(
					&abs_srcBlock, &ibSrc, &ieSrc, &jbSrc, &jeSrc);
			iGlobal = CppPOPBlocksMod::POP_HaloMsgCreate_iGlobal;
		}
		if (dstBlock != 0) {
			int abs_dstBlock = abs(dstBlock);
			CppPOPBlocksMod::pop_blocks_get_block_info1(
					&abs_dstBlock, &ibDst, &ieDst, &jbDst, &jeDst);
		}
	}

	// if both blocks are local, create a local copy to fill halo
	if ((srcProc == (POP_myTask + 1)) && (dstProc == (POP_myTask + 1))) {

    // compute addresses based on direction
		msgIndx = halo.numLocalCopies;

    // if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
    //     msgIndx > size(halo%dstLocalAddr,dim=2)) then
  	// srclocaladdr[numlocalcopies][3]
 		// dstLocalAddr[numLocalCopies][3]
		if (msgIndx > halo.numLocalCopies) {
			printf("POP_HaloMsgCreate: msg count > array size\n");
			exit(0);
		}

		if (str_trim_cmp(direction, "east") == 0) {
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
        	// copy easternmost physical domain of src
        	// into westernmost halo of dst
					halo.srcLocalAddr[msgIndx][0] = ieSrc - POP_HALO_WIDTH + i;
					halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
					halo.srcLocalAddr[msgIndx][2] = srcLocalID;

					halo.dstLocalAddr[msgIndx][0] = i;
					halo.dstLocalAddr[msgIndx][1] = jbDst + j - 1;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "west") == 0) {
      // copy westernmost physical domain of src
      // into easternmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
					halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
					halo.srcLocalAddr[msgIndx][2] = srcLocalID;

					halo.dstLocalAddr[msgIndx][0] = ieDst + i;
					halo.dstLocalAddr[msgIndx][1] = jbDst + j - 1;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "north") == 0) {
      // copy northern physical domain of src
      // into southern halo of dst
			if (srcBlock > 0 && dstBlock > 0) {
				// normal north boundary
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
						halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
						halo.srcLocalAddr[msgIndx][2] = srcLocalID;

						halo.dstLocalAddr[msgIndx][0] = ibDst + i - 1;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			} else if (srcBlock > 0 && dstBlock < 0) {
      	// tripole grid - copy info into tripole buffer
      	// copy physical domain of top halo+1 rows
      	// into global buffer at src location

	      // perform an error check to make sure the
 	    	// block has enough points to perform a tripole
  	    // update
				if ((jeSrc - jbSrc) < POP_HALO_WIDTH) {
					printf ("POP_HaloMsgCreate: not enough points in block for tripole\n");
					exit(0);
				}
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
						halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
						halo.srcLocalAddr[msgIndx][2] = srcLocalID;

						halo.dstLocalAddr[msgIndx][0] = iGlobal[ibSrc + i - 2];
						halo.dstLocalAddr[msgIndx][1] = j;
						halo.dstLocalAddr[msgIndx][2] = - dstLocalID;
						msgIndx += 1;
					}
				}
			} else if (srcBlock < 0 && dstBlock > 0) {
      	// tripole grid - set up for copying out of 
      	// tripole buffer into ghost cell domains
      	// include e-w ghost cells
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc + POP_HALO_WIDTH; ++i) {
						halo.srcLocalAddr[msgIndx][0] = nxGlobal - iGlobal[i-1] + 1;
						halo.srcLocalAddr[msgIndx][1] = POP_HALO_WIDTH + 3 - j;
						halo.srcLocalAddr[msgIndx][2] = - srcLocalID;

						halo.dstLocalAddr[msgIndx][0] = i;
						halo.dstLocalAddr[msgIndx][1] = j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			}
		} else if (str_trim_cmp(direction, "south") == 0) {
      // copy southern physical domain of src
      // into northern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
					halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
					halo.srcLocalAddr[msgIndx][1] = jeSrc - POP_HALO_WIDTH + j;
					halo.srcLocalAddr[msgIndx][2] = srcLocalID;

					halo.dstLocalAddr[msgIndx][0] = ibDst + i - 1;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "northeast") == 0) {
			if (dstBlock > 0) {
        // normal northeast boundary - just copy NE corner
        // of physical domain into SW halo of NE nbr block
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.srcLocalAddr[msgIndx][0] = ieSrc - POP_HALO_WIDTH + i;
						halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
						halo.srcLocalAddr[msgIndx][2] = srcLocalID;

						halo.dstLocalAddr[msgIndx][0] = i;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			} else {
        // tripole grid - this local copy should already
        // have taken place for the north boundary
			}
		} else if (str_trim_cmp(direction, "northwest") == 0) {
			if (dstBlock > 0) {
        // normal northeast boundary - just copy NW corner
      	// of physical domain into SE halo of NW nbr block
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
						halo.srcLocalAddr[msgIndx][1] = jbSrc + j - 1;
						halo.srcLocalAddr[msgIndx][2] = srcLocalID;

						halo.dstLocalAddr[msgIndx][0] = ieDst + i;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			} else {
        // tripole grid - this local copy should already
        // have taken place for the north boundary
			}
		} else if (str_trim_cmp(direction, "southeast") == 0) {
      // copy southeastern corner of src physical domain
      // into northwestern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = ieSrc - POP_HALO_WIDTH + i;
					halo.srcLocalAddr[msgIndx][1] = jeSrc - POP_HALO_WIDTH + j;
					halo.srcLocalAddr[msgIndx][2] = srcLocalID;

					halo.dstLocalAddr[msgIndx][0] = i;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "southwest") == 0) {
      // copy southwestern corner of src physical domain
      // into northeastern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = ibSrc + i - 1;
					halo.srcLocalAddr[msgIndx][1] = jeSrc - POP_HALO_WIDTH + j;
					halo.srcLocalAddr[msgIndx][2] = srcLocalID;

					halo.dstLocalAddr[msgIndx][0] = ieDst + i;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else {
			printf("POP_HaloMsgCreate: unknown direction local copy\n");
			exit(0);
		}
		halo.numLocalCopies = msgIndx;
		// TODO
      // if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
      //     msgIndx > size(halo%dstLocalAddr,dim=2)) then
      //    call LICOM_ErrorSet(errorCode, &
      //       'POP_HaloMsgCreate: msg count > array size')
      //    return
      // endif
	// -----------------------------------------------
	// if dest block is local and source block does not exist, create a 
	// local copy to fill halo with a fill value
	// -----------------------------------------------
	} else if ((srcProc == 0) && (dstProc == (POP_myTask + 1))) {
		msgIndx = halo.numLocalCopies;
			// TODO
      // if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
      //     msgIndx > size(halo%dstLocalAddr,dim=2)) then
      //    call LICOM_ErrorSet(errorCode, &
      //       'POP_HaloMsgCreate: msg count > array size')
      //    return
      // endif

    // compute addresses based on direction
		if (str_trim_cmp(direction, "east") == 0) {
      // copy easternmost physical domain of src
      // into westernmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = 0;
					halo.srcLocalAddr[msgIndx][1] = 0;
					halo.srcLocalAddr[msgIndx][2] = 0;

					halo.dstLocalAddr[msgIndx][0] = i;
					halo.dstLocalAddr[msgIndx][1] = jbDst + j - i;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "west") == 0) {
      // copy westernmost physical domain of src
      // into easternmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = 0;
					halo.srcLocalAddr[msgIndx][1] = 0;
					halo.srcLocalAddr[msgIndx][2] = 0;

					halo.dstLocalAddr[msgIndx][0] = ieDst + i;
					halo.dstLocalAddr[msgIndx][1] = jbDst + j - i;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "north") == 0) {
      // copy northern physical domain of src
      // into southern halo of dst
			if (dstBlock > 0) {
				// normal north boundary
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.srcLocalAddr[msgIndx][0] = 0;
						halo.srcLocalAddr[msgIndx][1] = 0;
						halo.srcLocalAddr[msgIndx][2] = 0;

						halo.dstLocalAddr[msgIndx][0] = ibDst + i - 1;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			}
		} else if (str_trim_cmp(direction, "south") == 0) {
      // copy southern physical domain of src
      // into northern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
					halo.srcLocalAddr[msgIndx][0] = 0;
					halo.srcLocalAddr[msgIndx][1] = 0;
					halo.srcLocalAddr[msgIndx][2] = 0;

					halo.dstLocalAddr[msgIndx][0] = jbDst + i - 1;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "northeast") == 0) {
      // normal northeast boundary - just copy NE corner
      // of physical domain into SW halo of NE nbr block
			if (dstBlock > 0) {
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.srcLocalAddr[msgIndx][0] = 0;
						halo.srcLocalAddr[msgIndx][1] = 0;
						halo.srcLocalAddr[msgIndx][2] = 0;
     
						halo.dstLocalAddr[msgIndx][0] = i;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			}
		} else if (str_trim_cmp(direction, "northwest") == 0) {
      // normal northeast boundary - just copy NW corner
      // of physical domain into SE halo of NW nbr block
			if (dstBlock > 0) {
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.srcLocalAddr[msgIndx][0] = 0;
						halo.srcLocalAddr[msgIndx][1] = 0;
						halo.srcLocalAddr[msgIndx][2] = 0;
     
						halo.dstLocalAddr[msgIndx][0] = ieDst + i;
						halo.dstLocalAddr[msgIndx][1] = jeDst + j;
						halo.dstLocalAddr[msgIndx][2] = dstLocalID;
						msgIndx += 1;
					}
				}
			}
		} else if (str_trim_cmp(direction, "southeast") == 0) {
      // copy southeastern corner of src physical domain
      // into northwestern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = 0;
					halo.srcLocalAddr[msgIndx][1] = 0;
					halo.srcLocalAddr[msgIndx][2] = 0;
    
					halo.dstLocalAddr[msgIndx][0] = i;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else if (str_trim_cmp(direction, "southwest") == 0) {
      // copy southwestern corner of src physical domain
      // into northeastern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.srcLocalAddr[msgIndx][0] = 0;
					halo.srcLocalAddr[msgIndx][1] = 0;
					halo.srcLocalAddr[msgIndx][2] = 0;
    
					halo.dstLocalAddr[msgIndx][0] = ieDst + i;
					halo.dstLocalAddr[msgIndx][1] = j;
					halo.dstLocalAddr[msgIndx][2] = dstLocalID;
					msgIndx += 1;
				}
			}
		} else {
			printf ("POP_HaloMsgCreate: unknown direction local copy\n");
			exit(0);
		}
		halo.numLocalCopies = msgIndx;
			// TODO
      // if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
      //     msgIndx > size(halo%dstLocalAddr,dim=2)) then
      //    call LICOM_ErrorSet(errorCode, &
      //       'POP_HaloMsgCreate: msg count > array size')
      //    return
      // endif
	// -----------------------------------
	// if source block local and dest block remote, send a message
	// -----------------------------------
	} else if ((srcProc == (POP_myTask + 1)) && (dstProc != (POP_myTask + 1)) &&
						 (dstProc > 0)) {
    // first check to see if a message to this processor has
    // already been defined
    // if not, update counters and indices
		msgIndx = -1;
		for (int n = 1; n <= halo.numMsgSend; ++n) {
			if (halo.sendTask[n - 1] == (dstProc - 1)) {
				msgIndx = n;
				bufSize = halo.sizeSend[n - 1];
				break;
			}
		}
		if (msgIndx == -1) {
			msgIndx = halo.numMsgSend + 1;
			halo.numMsgSend = msgIndx;
			halo.sendTask[msgIndx - 1] = dstProc - 1;
			bufSize = 0;
		}
		int msgIndex = msgIndx - 1;
    // now compute message info based on msg direction
		if (str_trim_cmp (direction, "east") == 0) {
      // send easternmost physical domain of src
      // into westernmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
							= ieSrc - POP_HALO_WIDTH + i;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
							= jbSrc + j - 1;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
							= srcLocalID;
					bufSize += 1;
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "west") == 0) {
      // copy westernmost physical domain of src
      // into easternmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
							= ibSrc + i - 1;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
							= jbSrc + j - 1;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
							= srcLocalID;

					bufSize += 1;
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "north") == 0) {
			if (dstBlock > 0) {
        // copy northern physical domain of src
        // into southern halo of dst
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ibSrc + i - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - send top three rows of phys domain
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ibSrc + i - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "south") == 0) {
      // copy southern physical domain of src
      // into northern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
							= ibSrc + i - 1;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
							= jeSrc - POP_HALO_WIDTH + j;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
							= srcLocalID;
					bufSize += 1;
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "northeast") == 0) {
      // copy northeast corner of src physical domain
      // into southwestern halo of dst
			if (dstBlock > 0) {
				// normal northeast corner
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ieSrc - POP_HALO_WIDTH + i;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - send top three rows of phys domain
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ibSrc + i - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "northwest") == 0) {
      // copy northwest corner of src physical domain
      // into southeastern halo of dst
			if (dstBlock > 0) {
        // normal northwest corner
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ibSrc + i - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - send top three rows of phys domain
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
								= ibSrc + i - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
								= jbSrc + j - 1;
						halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
								= srcLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "southeast") == 0) {
      // copy southeastern corner of src physical domain
      // into northwestern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
							= ieSrc - POP_HALO_WIDTH + i;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
							= jeSrc - POP_HALO_WIDTH + j;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
							= srcLocalID;
					bufSize += 1;
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "southwest") == 0) {
      // copy southwestern corner of src physical domain
      // into northeastern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][0] 
							= ibSrc + i - 1;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][1] 
							= jeSrc - POP_HALO_WIDTH + j;
					halo.sendAddr[msgIndex*bufSizeSend + bufSize][2] 
							= srcLocalID;
					bufSize += 1;
				}
			}
			halo.sizeSend[msgIndex] = bufSize;
		}
	// -----------------------------------------
	// if source block remote and dest block local, recv a message
	// -----------------------------------------
	} else if ((dstProc == (POP_myTask + 1)) &&
						 (srcProc != (POP_myTask + 1)) &&
						 (srcProc > 0)) {
    // first check to see if a message from this processor has
    // already been defined
    // if not, update counters and indices
		msgIndx = -1;
		for (int n = 1; n <= halo.numMsgRecv; ++n) {
			if (halo.recvTask[n-1] == (srcProc - 1)) {
				msgIndx = n;
				bufSize = halo.sizeRecv[n-1];
				break;
			}
		}
		if (msgIndx == -1) {
			msgIndx = halo.numMsgRecv + 1;
			halo.numMsgRecv = msgIndx;
			halo.recvTask[msgIndx - 1] = srcProc - 1;
			bufSize = 0;
		}
		int msgIndex = msgIndx - 1;
    // now compute message info based on msg direction
		if (str_trim_cmp (direction, "east") == 0) {
      // send easternmost physical domain of src
      // into westernmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
							= i;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 							= jbDst + j - 1;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
							= dstLocalID;
					bufSize += 1;
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "west") == 0) {
      // copy westernmost physical domain of src
      // into easternmost halo of dst
			for (j = 1; j <= jeSrc - jbSrc + 1; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
							= ieDst + i;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 							= jbDst + j - 1;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
							= dstLocalID;
					bufSize += 1;
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "north") == 0) {
      // copy northern physical domain of src
      // into southern halo of dst
			if (dstBlock > 0) {
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= ieDst - ibDst + 1; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= ibDst + i - 1;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= jeDst + j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= dstLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - receive into tripole buffer
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= iGlobal[ibSrc + i - 2];
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= - dstLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "south") == 0) {
      // copy southern physical domain of src
      // into northern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
							= ibDst + i - 1;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 							= j;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
							= dstLocalID;
					bufSize += 1;
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "northeast") == 0) {
      // copy northeast physical domain into
      // into southwest halo of dst
			if (dstBlock > 0) {
        // normal northeast neighbor
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= i;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= jeDst + j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= dstLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - receive into tripole buffer
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - ibSrc + 1; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= iGlobal[ibSrc + i - 2];
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= - dstLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "northwest") == 0) {
      // copy northwest physical domain into
      // into southeast halo of dst
			if (dstBlock > 0) {
        // normal northwest neighbor
				for (j = 1; j <= POP_HALO_WIDTH; ++j) {
					for (i = 1; i <= POP_HALO_WIDTH; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= ieDst + i;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= jeDst + j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= dstLocalID;
						bufSize += 1;
					}
				}
			} else {
        // tripole block - receive into tripole buffer
				for (j = 1; j <= POP_HALO_WIDTH + 1; ++j) {
					for (i = 1; i <= ieSrc - jbSrc + 1; ++i) {
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
								= iGlobal[ibSrc + i - 2];
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 								= j;
						halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
								= - dstLocalID;
						bufSize += 1;
					}
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "southeast") == 0) {
      // copy southeastern corner of src physical domain
      // into northwestern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
							= i;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 							= j;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
							= dstLocalID;
					bufSize += 1;
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		} else if (str_trim_cmp(direction, "southwest") == 0) {
      // copy southwestern corner of src physical domain
      // into northeastern halo of dst
			for (j = 1; j <= POP_HALO_WIDTH; ++j) {
				for (i = 1; i <= POP_HALO_WIDTH; ++i) {
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][0] 
							= ieDst + i;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][1] 
 							= j;
					halo.recvAddr[msgIndex*bufSizeRecv + bufSize][2] 
							= dstLocalID;
					bufSize += 1;
				}
			}
			halo.sizeRecv[msgIndex] = bufSize;
		}
	}
	// ---------------------------------
	// if none of the cases above, no message info required for this
	// block pair
	// ---------------------------------
	return ;
}

void pop_halo_create (
    const CppPOPHaloMod::pop_halo &f_halo,
    CppPOPHaloMod::pop_halo &c_halo) {
  c_halo.communicator   = f_halo.communicator;
  c_halo.numMsgSend     = f_halo.numMsgSend;
  c_halo.numMsgRecv     = f_halo.numMsgRecv;
  c_halo.numLocalCopies = f_halo.numLocalCopies;

	c_halo.recvTask = f_halo.recvTask;
	c_halo.sendTask = f_halo.sendTask;

	c_halo.sizeSend = f_halo.sizeSend;
	c_halo.sizeRecv = f_halo.sizeRecv;

	c_halo.srcLocalAddr = (int(*)[3])malloc(c_halo.numLocalCopies * sizeof(int[3]));
	c_halo.dstLocalAddr = (int(*)[3])malloc(c_halo.numLocalCopies * sizeof(int[3]));

	c_halo.sendAddr = (int(*)[3])malloc(c_halo.numMsgSend * bufSizeSend * sizeof(int[3]));
	c_halo.recvAddr = (int(*)[3])malloc(c_halo.numMsgRecv * bufSizeRecv * sizeof(int[3]));
	return ;
}

void pop_halo_create (const CppPOPDistributionMod::POP_distrb &distrb,
		const char *nsBoundaryType, const char *ewBoundaryType,
				const int &nxGlobal, CppPOPHaloMod::pop_halo &halo) {
	// DESCRIPTION:
	// This routine creates a halo type with info necessary for
	// performing a halo (ghost cell) update. This info is computed
	// based on the input block distribution.
	// INPUT:
	// distrb: distribution of blocks across procs
	// nsBoundaryType: type of boundary to use in logical ns dir
	// ewBoundaryType: type of boundary to use in logical ew dir
	// nxGlobal: global grid extent for tripole grids
	// OUTPUT:
	// halo: a new halo type with info for halo updates

	using CppPOPCommMod::POP_myTask;
  using CppPOPBlocksMod::POP_HALO_WIDTH;
	using CppPOPBlocksMod::POP_NX_BLOCK;
	using CppPOPBlocksMod::POP_NY_BLOCK;
	using CppPOPBlocksMod::POP_numBlocks;
	using CppPOPBlocksMod::pop_blocks_get_nbr_id;
	using CppPOPDistributionMod::pop_distribution_get_block_loc;

	// allocate status flag
  // int istat;
	// block id east,  west neighbors
  int eastBlock, westBlock;
	// block id north, south neighbors
  int northBlock, southBlock;
	// block id northeast, northwest nbrs
  int neBlock, nwBlock;
	// block id southeast, southwest nbrs
  int seBlock, swBlock;
	// source, dest processor locations
  int srcProc, dstProc;
	// local block index of src, dst blocks
  int srcLocalID, dstLocalID;
	// max buffer sizes
  int maxSizeSend, maxSizeRecv;
	// number of messages for this halo
  int numMsgSend, numMsgRecv;

	// flag for identifying north tripole blocks
	bool tripoleBlock;

	// distribution of blocks across procs
	// auto &distrb = *p_distrb;

	//------------------------------------
	// Initialize some useful variables and return if this task 
	// not in the current distribution.
	//------------------------------------

	// num of processors involved
  int numProcs;
	// communicator for message passing
  int communicator;
	CppPOPDistributionMod::pop_distribution_get(distrb, numProcs, communicator);

	if (POP_myTask >= numProcs) return;

	halo.communicator = communicator;

	// size of default physical domain in X
	int blockSizeX = POP_NX_BLOCK - 2 * POP_HALO_WIDTH;
	// size of default physical domain in Y
	int blockSizeY = POP_NY_BLOCK - 2 * POP_HALO_WIDTH;

	// nominal sizes for e-w msgs
	int eastMsgSize = POP_HALO_WIDTH * blockSizeY;
	int westMsgSize = POP_HALO_WIDTH * blockSizeY;

	// nominal sizes for n-s msgs
	int northMsgSize = POP_HALO_WIDTH * blockSizeX;
	int southMsgSize = POP_HALO_WIDTH * blockSizeX;

	// nominal size for corner msg
  int cornerMsgSize = POP_HALO_WIDTH * POP_HALO_WIDTH;
	int msgSize = 0;
	// size for tripole messages
  int tripoleMsgSize = (POP_HALO_WIDTH + 1) * blockSizeX;
	// size for tripole messages
  int tripoleMsgSizeOut = (POP_HALO_WIDTH + 1) * POP_NX_BLOCK;

	// flag for allocating tripole buffers
	// bool tripoleFlag;
	if (str_trim_cmp(nsBoundaryType, "tripole") == 0) {
		// tripoleFlag = true;
  	// allocate tripole message buffers if not already done
		if (bufTripole == nullptr) {
			bufTripole = new double(nxGlobal * (POP_HALO_WIDTH + 1));
		}
	} else {
		// tripoleFlag = false;
	}
	// --------------------------------------------------
	// Count the number of messages to send/recv from each processor
	// and number of words in each message.  These quantities are
	// necessary for allocating future arrays.
	// --------------------------------------------------

	// count number of words to each proc
	int *sendCount = sendCount = new int[numProcs];
	int *recvCount = recvCount = new int[numProcs];

	memset(sendCount, 0, sizeof(int) * numProcs);
	memset(recvCount, 0, sizeof(int) * numProcs);


	// msgCountLoop
	for (int iblock = 1; iblock <= POP_numBlocks; ++iblock) {

		pop_distribution_get_block_loc(distrb, iblock, srcProc, srcLocalID);

  	// find north neighbor block and add to message count
  	// also set tripole block flag for later special cases
		using CppPOPBlocksMod::POP_BLOCKS_NORTH;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH, ewBoundaryType, 
				nsBoundaryType, northBlock);
		if (northBlock > 0) {
			tripoleBlock = false;
			msgSize = northMsgSize;
			pop_distribution_get_block_loc(distrb, northBlock, dstProc, dstLocalID);
		} else if (northBlock < 0) {
			// tripole north row, count block
			tripoleBlock = true;
			msgSize = tripoleMsgSize;
			pop_distribution_get_block_loc(distrb, abs(northBlock), dstProc, dstLocalID);
		} else {
			tripoleBlock = false;
			msgSize      = northMsgSize;
			dstProc 	   = 0;
			dstLocalID   = 0;
		}
		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, msgSize, numProcs);

    // if a tripole boundary block, also create a local
    // message into and out of tripole buffer 
		if (tripoleBlock) {
      // copy out of tripole buffer - includes halo
			pop_halo_increment_msg_count(sendCount, recvCount,
					srcProc, srcProc, tripoleMsgSizeOut, numProcs);

      // copy in only required if dstProc not same as srcProc
			if (dstProc != srcProc) {
				pop_halo_increment_msg_count(sendCount, recvCount,
						srcProc, srcProc, msgSize, numProcs);
			}
		}

    // find south neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH,
				ewBoundaryType, nsBoundaryType, southBlock);

		if (southBlock > 0) {
			pop_distribution_get_block_loc(distrb, southBlock,
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, southMsgSize, numProcs);

    // find east neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_EAST,
				ewBoundaryType, nsBoundaryType, eastBlock);
		if (eastBlock > 0) {
			pop_distribution_get_block_loc(distrb, eastBlock, 
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, eastMsgSize, numProcs);

    // if a tripole boundary block, non-local east neighbor
    // needs a chunk of the north boundary, so add a message
    // for that
		if (tripoleBlock && dstProc != srcProc) {
			pop_halo_increment_msg_count(sendCount, recvCount,
					srcProc, dstProc, tripoleMsgSize, numProcs);
		}

    // find west neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST, 
				ewBoundaryType, nsBoundaryType, westBlock);

		if (westBlock > 0) {
			pop_distribution_get_block_loc(distrb, westBlock,
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}
		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, westMsgSize, numProcs);

    // if a tripole boundary block, non-local west neighbor
    // needs a chunk of the north boundary, so add a message
    // for that
		if (tripoleBlock && dstProc != srcProc) {
			pop_halo_increment_msg_count(sendCount, recvCount,
					srcProc, dstProc, tripoleMsgSize, numProcs);
		}
    // find northeast neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_NORTH_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH_EAST,
				ewBoundaryType, nsBoundaryType, neBlock);
		if (neBlock > 0) {
			// normal corner message 
			msgSize = cornerMsgSize;
			pop_distribution_get_block_loc(distrb, neBlock,
					dstProc, dstLocalID);
		} else if (neBlock < 0) {
			// tripole north row
			// tripole needs whole top row of block
			msgSize = tripoleMsgSize;
			pop_distribution_get_block_loc(distrb, abs(neBlock),
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}
		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, msgSize, numProcs);

    // find northwest neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_NORTH_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH_WEST,
				ewBoundaryType, nsBoundaryType, nwBlock);
		if (nwBlock > 0) {
			// normal NE corner update
			msgSize = cornerMsgSize;
			pop_distribution_get_block_loc(distrb, nwBlock,
					dstProc, dstLocalID);
		} else if (nwBlock < 0) {
			// tripole north row, count block
			// tripole NE corner update - entire row needed
			msgSize = tripoleMsgSize;
			pop_distribution_get_block_loc(distrb, abs(nwBlock),
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}
		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, msgSize, numProcs);
  	// find southeast neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH_EAST,
				ewBoundaryType, nsBoundaryType, seBlock);
		if (seBlock > 0) {
			pop_distribution_get_block_loc(distrb, seBlock, 
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, cornerMsgSize, numProcs);

  	// find southwest neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH_WEST,
				ewBoundaryType, nsBoundaryType, swBlock);

		if (swBlock > 0) {
			pop_distribution_get_block_loc(distrb, swBlock,
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}
		pop_halo_increment_msg_count(sendCount, recvCount,
				srcProc, dstProc, cornerMsgSize, numProcs);

    // for tripole grids with padded domain, padding will
    // prevent tripole buffer from getting all the info
    // it needs - must extend footprint at top boundary
		if (tripoleBlock && ((nxGlobal%blockSizeX) != 0)) {
			// tripole padding
      // find east2 neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_EAST2;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_EAST2,
					ewBoundaryType, nsBoundaryType, eastBlock);
			if (eastBlock > 0) {
				pop_distribution_get_block_loc(distrb, eastBlock,
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_increment_msg_count(sendCount, recvCount,
						srcProc, dstProc, tripoleMsgSize, numProcs);
			}
      // find EastNorthEast neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_EAST_NORTH_EAST;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_EAST_NORTH_EAST,
					ewBoundaryType, nsBoundaryType, neBlock);
			if (neBlock < 0) {
				// tripole north row
				// tripole needs whole top row of block
				msgSize = tripoleMsgSize;
				pop_distribution_get_block_loc(distrb, abs(neBlock),
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_increment_msg_count(sendCount, recvCount,
						srcProc, dstProc, msgSize, numProcs);
			}
      // find west2 neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_WEST2;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST2,
					ewBoundaryType, nsBoundaryType, westBlock);
			if (westBlock > 0) {
				pop_distribution_get_block_loc(distrb, westBlock,
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_increment_msg_count(sendCount, recvCount,
						srcProc, dstProc, tripoleMsgSize, numProcs);
			}
      // find WestNorthWest neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_WEST_NORTH_WEST;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST_NORTH_WEST,
					ewBoundaryType, nsBoundaryType, nwBlock);
			if (nwBlock < 0) {
				// tripole north row
				// tripole needs whole top row of block
				msgSize = tripoleMsgSize;
				pop_distribution_get_block_loc(distrb, abs(nwBlock),
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_increment_msg_count(sendCount, recvCount,
						srcProc, dstProc, msgSize, numProcs);
			}
		}
	} // End msgCountLoop

	// -------------------------------------------------
	// if messages are received from the same processor, the message is 
	// actually a local copy - count them and reset to zero
	// -------------------------------------------------
	halo.numLocalCopies = recvCount[POP_myTask];
	sendCount[POP_myTask] = 0;
	recvCount[POP_myTask] = 0;

	// -------------------------------------------------
	// now count the number of actual messages to be sent and received
	// -------------------------------------------------
	numMsgSend = 0;
	numMsgRecv = 0;
	for (int n = 0; n < numProcs; ++n) {
		if (sendCount[n] != 0) numMsgSend += 1;
		if (recvCount[n] != 0) numMsgRecv += 1;
	}
	halo.numMsgSend = numMsgSend;
	halo.numMsgRecv = numMsgRecv;

	// -------------------------------------------------
	// allocate buffers for 2-d halo updates to save time later
	// if the buffers are already allocated by previous create call,
	// check to see if they need to be re-sized
	// -------------------------------------------------
	// temp for global maxval
  int maxTmp = 0;
	for (int n = 0; n < numProcs; ++n) {
		if (sendCount[n] >= maxTmp) maxTmp = sendCount[n];
	}
	using CppPOPReductionsMod::pop_global_maxval_scalar;
	pop_global_maxval_scalar(maxTmp, distrb, maxSizeSend);

  maxTmp = 0;
	for (int n = 0; n < numProcs; ++n) {
		if (recvCount[n] >= maxTmp) maxTmp = recvCount[n];
	}
	pop_global_maxval_scalar(maxTmp, distrb, maxSizeRecv);

	// TODO
	/*
	// flag for resizing buffers
	bool resize;
	if (bufSendR8 == nullptr && bufRecvR8 == nullptr) {
		bufSizeSend = maxSizeSend;
		bufSizeRecv = maxSizeRecv;

		bufSendR8 = new double[bufSizeSend * numMsgSend];
		bufRecvR8 = new double[bufSizeRecv * numMsgRecv];
	} else {
		resize = false;
		if (maxSizeSend > bufSizeSend) {
			resize = true;
			bufSizeSend = maxSizeSend;
		}
		if (maxSizeRecv > bufSizeRecv) {
			resize = true;
			bufSizeRecv = maxSizeRecv;
		}

    // if (numMsgSend > size(bufSendR8,dim=2)) resize = .true.
    // if (numMsgRecv > size(bufRecvR8,dim=2)) resize = .true.
		if (resize) {
			delete [] bufSendR8;
			delete [] bufRecvR8;

			bufSendR8 = nullptr;
			bufRecvR8 = nullptr;

			bufSendR8 = new double[bufSizeSend * numMsgSend];
			bufRecvR8 = new double[bufSizeRecv * numMsgRecv];
		}
	}
	*/
	// ----------------------------------------------------
	// allocate arrays for message information and initialize
	// ----------------------------------------------------
	halo.sendTask = new int[numMsgSend];
	halo.recvTask = new int[numMsgRecv];
	halo.sizeSend = new int[numMsgSend];
	halo.sizeRecv = new int[numMsgRecv];
	halo.sendAddr = new int[numMsgSend * bufSizeSend][3];
	halo.recvAddr = new int[numMsgRecv * bufSizeRecv][3];
	halo.srcLocalAddr = new int[halo.numLocalCopies][3];
	halo.dstLocalAddr = new int[halo.numLocalCopies][3];

	memset(halo.sendTask, 0, sizeof(int) * numMsgSend);
	memset(halo.recvTask, 0, sizeof(int) * numMsgSend);
	memset(halo.sizeSend, 0, sizeof(int) * numMsgSend);
	memset(halo.sizeRecv, 0, sizeof(int) * numMsgSend);

	memset(halo.sendAddr, 0, sizeof(int[3]) * numMsgSend * bufSizeSend);
	memset(halo.recvAddr, 0, sizeof(int[3]) * numMsgSend * bufSizeRecv);

	memset(halo.srcLocalAddr, 0, sizeof(int[3]) * halo.numLocalCopies);
	memset(halo.dstLocalAddr, 0, sizeof(int[3]) * halo.numLocalCopies);

	delete [] sendCount;
	delete [] recvCount;

	sendCount = nullptr;
	recvCount = nullptr;

	// ---------------------------------------------------
	// repeat loop through blocks but this time, determine all the
	// required message information for each message or local copy
	// ---------------------------------------------------

  // reset halo scalars to use as counters
	halo.numMsgSend     = 0;
	halo.numMsgRecv     = 0;
	halo.numLocalCopies = 0;

	// msgConfigLoop
	for (int iblock = 1; iblock <= POP_numBlocks; ++iblock) {

		pop_distribution_get_block_loc(distrb, iblock, srcProc, srcLocalID);

    // find north neighbor block and set msg info
    // also set tripole block flag for later special cases
		using CppPOPBlocksMod::POP_BLOCKS_NORTH;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH,
				ewBoundaryType, nsBoundaryType, northBlock);
		if (northBlock > 0) {
			tripoleBlock = false;
			pop_distribution_get_block_loc(distrb, northBlock,
					dstProc, dstLocalID);
		} else if (northBlock < 0) {
			// tripole north row, count block
			tripoleBlock = true;
			pop_distribution_get_block_loc(distrb, abs(northBlock),
					dstProc, dstLocalID);
		} else {
			tripoleBlock = false;
			dstProc      = 0;
			dstLocalID   = 0;
		}
		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				northBlock, dstProc, dstLocalID, nxGlobal, "north");

    // if a tripole boundary block, also create a local
    // message into and out of tripole buffer 
		if (tripoleBlock) {
      // copy out of tripole buffer - includes halo
			// TODO
			int n_iblock = - iblock;
			pop_halo_msg_create (halo, n_iblock, srcProc, srcLocalID,
					iblock, srcProc, srcLocalID, nxGlobal, "north");
			
      // copy in only required if dstProc not same as srcProc
			if (dstProc != srcProc) {
				pop_halo_msg_create (halo, iblock, srcProc, srcLocalID,
						n_iblock, srcProc, srcLocalID, nxGlobal, "north");
			}
		}
    // find south neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH, ewBoundaryType,
				nsBoundaryType, southBlock);

		if (southBlock > 0) {
			pop_distribution_get_block_loc(distrb, southBlock, 
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				southBlock, dstProc, dstLocalID, nxGlobal, "south");
		
  // find east neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_EAST, ewBoundaryType,
				nsBoundaryType, eastBlock);

		if (eastBlock > 0) {
			pop_distribution_get_block_loc(distrb, eastBlock, dstProc,
					dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				eastBlock, dstProc, dstLocalID, nxGlobal, "east");

  	// if a tripole boundary block, non-local east neighbor
    // needs a chunk of the north boundary, so add a message
    // for that
		if (tripoleBlock && dstProc != srcProc) {
			pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
					-eastBlock, dstProc, dstLocalID, nxGlobal, "north");
		}
    // find west neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST, ewBoundaryType,
				nsBoundaryType, westBlock);

		if (westBlock > 0) {
			pop_distribution_get_block_loc(distrb, westBlock, dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				westBlock, dstProc, dstLocalID, nxGlobal, "west");

    // if a tripole boundary block, non-local west neighbor
    // needs a chunk of the north boundary, so add a message
    // for that
		if (tripoleBlock && dstProc != srcProc) {
			pop_halo_msg_create(halo, iblock, srcProc, srcLocalID, 
					-westBlock, dstProc, dstLocalID, nxGlobal, "north");
		}
    // find northeast neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_NORTH_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH_EAST,
				ewBoundaryType, nsBoundaryType, neBlock);

		if (neBlock != 0) {
			pop_distribution_get_block_loc(distrb, abs(neBlock),
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				neBlock, dstProc, dstLocalID, nxGlobal, "northeast");
  	// find northwest neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_NORTH_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_NORTH_WEST,
				ewBoundaryType, nsBoundaryType, nwBlock);
		if (nwBlock != 0) {
			pop_distribution_get_block_loc(distrb, abs(nwBlock),
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				nwBlock, dstProc, dstLocalID, nxGlobal, "northwest");
    // find southeast neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH_EAST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH_EAST,
				ewBoundaryType, nsBoundaryType, seBlock);

		if (seBlock > 0) {
			pop_distribution_get_block_loc(distrb, seBlock,
					dstProc, dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}

		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				seBlock, dstProc, dstLocalID, nxGlobal, "southeast");
    // find southwest neighbor block and add to message count
		using CppPOPBlocksMod::POP_BLOCKS_SOUTH_WEST;
		pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH_WEST,
				ewBoundaryType, nsBoundaryType, swBlock);

		if (swBlock > 0) {
			pop_distribution_get_block_loc(distrb, swBlock, dstProc,
					dstLocalID);
		} else {
			dstProc    = 0;
			dstLocalID = 0;
		}
		pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
				swBlock, dstProc, dstLocalID, nxGlobal, "southwest");

    // for tripole grids with padded domain, padding will
    // prevent tripole buffer from getting all the info
    // it needs - must extend footprint at top boundary
		if (tripoleBlock && ((nxGlobal%blockSizeX) != 0)) {
      // find east2 neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_SOUTH_EAST2;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_SOUTH_EAST2,
					ewBoundaryType, nsBoundaryType, eastBlock);
			if (eastBlock > 0) {
				pop_distribution_get_block_loc(distrb, eastBlock,
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
						-eastBlock, dstProc, dstLocalID, nxGlobal, "north");
			}
      // find EastNorthEast neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_EAST_NORTH_EAST;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_EAST_NORTH_EAST,
					ewBoundaryType, nsBoundaryType, neBlock);
			if (neBlock < 0) {
				// tripole north row
				// tripole needs whole top row of block
				msgSize = tripoleMsgSize;
				pop_distribution_get_block_loc(distrb, abs(neBlock),
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
						neBlock, dstProc, dstLocalID, nxGlobal, "north");
			}
      // find west2 neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_WEST2;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST2, 
					ewBoundaryType, nsBoundaryType, westBlock);

			if (westBlock > 0) {
				pop_distribution_get_block_loc(distrb, westBlock,
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
						-westBlock, dstProc, dstLocalID, nxGlobal, "north");
			}
      // find WestNorthWest neighbor block and add to message count
			using CppPOPBlocksMod::POP_BLOCKS_WEST_NORTH_WEST;
			pop_blocks_get_nbr_id(iblock, POP_BLOCKS_WEST_NORTH_WEST, 
					ewBoundaryType, nsBoundaryType, nwBlock);
			if (nwBlock < 0) {
				// tripole north row
				// tripole needs whole top row of block
				msgSize = tripoleMsgSize;
				pop_distribution_get_block_loc(distrb, abs(nwBlock),
						dstProc, dstLocalID);
			} else {
				dstProc    = 0;
				dstLocalID = 0;
			}
			if (dstProc != srcProc) {
				pop_halo_msg_create(halo, iblock, srcProc, srcLocalID,
						nwBlock, dstProc, dstLocalID, nxGlobal, "north");
			}
		}
	} // msgConfigLoop
	return ;
}

// extern int* work_mod_mp_test_recvtask_;
// extern int* work_mod_mp_test_sendtask_;
// extern int* work_mod_mp_test_sizerecv_;
// extern int* work_mod_mp_test_sizesend_;
// extern int* work_mod_mp_test_srclocaladdr_;
// extern int* work_mod_mp_test_dstlocaladdr_;
// extern int* work_mod_mp_test_sendaddr_;
// extern int* work_mod_mp_test_recvaddr_;

void pop_halo_create_from_fortran(CppPOPHaloMod::POPHalo &halo) {
	int comm;
	pop_halo_create_from_fortran_get_1_(&comm, 
			&(halo.numMsgSend), &(halo.numMsgRecv), &(halo.numLocalCopies));

	halo.communicator = MPI_Comm_f2c(comm);

	if (halo.sendTask == nullptr) {
		halo.sendTask = new int[halo.numMsgSend];
	}
	if (halo.recvTask == nullptr) {
		halo.recvTask = new int[halo.numMsgRecv];
	}
	if (halo.sizeSend == nullptr) {
		halo.sizeSend = new int[halo.numMsgSend];
	}
	if (halo.sizeRecv == nullptr) {
		halo.sizeRecv = new int[halo.numMsgRecv];
	}
	if (halo.srcLocalAddr == nullptr) {
		halo.srcLocalAddr = new int[halo.numLocalCopies][3];
	}
	if (halo.dstLocalAddr == nullptr) {
		halo.dstLocalAddr = new int[halo.numLocalCopies][3];
	}
	if (halo.sendAddr == nullptr) {
		halo.sendAddr = new int[halo.numMsgSend * bufSizeSend][3];
	}
	if (halo.recvAddr == nullptr) {
		halo.recvAddr = new int[halo.numMsgRecv * bufSizeRecv][3];
	}
	// printf ("C %d %d %d %d %d %d %d %d\n", 
	// halo.sendTask, halo.recvTask,
	// 		halo.sizeSend, halo.sizeRecv, &halo.srcLocalAddr[0][0], &halo.dstLocalAddr[0][0],
	// 				&halo.sendAddr[0][0], &halo.recvAddr[0][0]
	// );
	pop_halo_create_from_fortran_get_2_ (halo.sendTask, halo.recvTask,
			halo.sizeSend, halo.sizeRecv, &halo.srcLocalAddr[0][0], &halo.dstLocalAddr[0][0],
					&halo.sendAddr[0][0], &halo.recvAddr[0][0]);
	// for (int i = 0; i < halo.numMsgSend; ++i) {
	// 	halo.sendTask[i] = POPHalo_bufSendTask[i];
	// 	halo.sizeSend[i] = POPHalo_bufSizeSend[i];
	// }
	// for (int i = 0; i < halo.numMsgRecv; ++i) {
	// 	halo.recvTask[i] = POPHalo_bufRecvTask[i];
	// 	halo.sizeRecv[i] = POPHalo_bufSizeRecv[i];
	// }
	// for (int j = 0; j < halo.numLocalCopies; ++j) {
	// 	for (int i = 0; i < 3; ++i) {
	// 		halo.srcLocalAddr[j][i] = POPHalo_bufSrcLocalAddr[j * 3 + i];
	// 		halo.dstLocalAddr[j][i] = POPHalo_bufDstLocalAddr[j * 3 + i];
	// 	}
	// }
	// for (int k = 0; k < halo.numMsgSend; ++k) {
	// 	for (int j = 0; j < bufSizeSend; ++j) {
	// 		for (int i = 0; i < 3; ++i) {
	// 			halo.sendAddr[k * bufSizeSend + j][i] = POPHalo_bufSendAddr[k*bufSizeSend*3 + j*3 + i];
	// 		}
	// 	}
	// }
	// for (int k = 0; k < halo.numMsgRecv; ++k) {
	// 	for (int j = 0; j < bufSizeRecv; ++j) {
	// 		for (int i = 0; i < 3; ++i) {
	// 			halo.recvAddr[k * bufSizeRecv + j][i] = POPHalo_bufRecvAddr[k*bufSizeRecv*3 + j*3 + i];
	// 		}
	// 	}
	// }
// if (CppParamMod::mytid == 0) {
//   printf("wjl debug pop_halo_create_from_fortran_free_buf_ 1\n");
// }
	// pop_halo_create_from_fortran_free_buf_();
// if (CppParamMod::mytid == 0) {
//   printf("wjl debug pop_halo_create_from_fortran_free_buf_ 2\n");
// }
	if (CppPOPHaloMod::sndRequest == nullptr) {
		CppPOPHaloMod::sndRequest = new MPI_Request[halo.numMsgSend];
	}
	if (CppPOPHaloMod::rcvRequest == nullptr) {
		CppPOPHaloMod::rcvRequest = new MPI_Request[halo.numMsgRecv];
	}
	if (CppPOPHaloMod::sndStatus == nullptr) {
		CppPOPHaloMod::sndStatus = new MPI_Status[halo.numMsgSend];
	}
	if (CppPOPHaloMod::rcvStatus == nullptr) {
		CppPOPHaloMod::rcvStatus = new MPI_Status[halo.numMsgRecv];
	}
	// malloc for buffers
	using CppParamMod::IMT_GLOBAL;
	using CppPOPBlocksMod::POP_HALO_WIDTH;
	// max 3D double in POP_HaloUpdate
	using CppParamMod::KM;
	if (bufTripole == nullptr) {
		bufTripole = new double[(POP_HALO_WIDTH + 1) * IMT_GLOBAL * KM];
	}
	if (bufSend == nullptr) {
		bufSend = new double[bufSizeSend * halo.numMsgSend * KM];
	}
	if (bufRecv == nullptr) {
		bufRecv = new double[bufSizeRecv * halo.numMsgRecv * KM];
	}

  // printf("c mytid %d, communicator %d, numMsgSend %d, numMsgRecv %d, numLocalCopies %d, bufSizeSend %d, bufSizeRecv %d\n", 
  // CppParamMod::mytid, halo.communicator, halo.numMsgSend, halo.numMsgRecv, halo.numLocalCopies, bufSizeSend, bufSizeRecv);
  // int recvtask = 0;
  // int sendtask = 0;
  // int sizerecv = 0;
  // int sizesend = 0;
  // int srclocaladdr = 0;
  // int dstlocaladdr = 0;
  // int sendaddr = 0;
  // int recvaddr = 0;
  // for (int i = 0; i < POP_haloClinic.numMsgSend; ++i) {
  //   sendtask += POP_haloClinic.sendTask[i] - work_mod_mp_test_sendtask_[i];
  //   sizesend += POP_haloClinic.sizeSend[i] - work_mod_mp_test_sizesend_[i];
  // }
  // for (int i = 0; i < POP_haloClinic.numMsgRecv; ++i) {
    // recvtask += POP_haloClinic.recvTask[i] - work_mod_mp_test_recvtask_[i];
    // sizerecv += POP_haloClinic.sizeRecv[i] - work_mod_mp_test_sizerecv_[i];
  // }
  // for (int i = 0; i < POP_haloClinic.numLocalCopies; ++i){
  //   for (int j = 0; j < 3; ++j) {
  //     srclocaladdr += POP_haloClinic.srcLocalAddr[i][j] - work_mod_mp_test_srclocaladdr_[i * 3 + j];
  //     dstlocaladdr += POP_haloClinic.dstLocalAddr[i][j] - work_mod_mp_test_dstlocaladdr_[i * 3 + j];
  //   }
  // }
  // using CppPOPHaloMod::bufSizeSend;
  // using CppPOPHaloMod::bufSizeRecv;
  // for (int i = 0; i < POP_haloClinic.numMsgSend; ++i){
  //   for (int k = 0; k < bufSizeSend; ++k) {
  //   for (int j = 0; j < 3; ++j) {
  //     sendaddr += POP_haloClinic.sendAddr[i * bufSizeSend + k][j] - 
  //     work_mod_mp_test_sendaddr_[i * bufSizeSend * 3 + k * 3+ j];
  //     recvaddr += POP_haloClinic.recvAddr[i * bufSizeSend + k][j] - 
  //     work_mod_mp_test_recvaddr_[i * bufSizeSend * 3 + k * 3+ j];
  //   }
  //   }
  // }
  // printf("mytid = %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
  // mytid, recvtask,sendtask,sizerecv,sizesend,
  // srclocaladdr,dstlocaladdr,sendaddr,recvaddr);
	return ;
}

// POP_HaloUpdate

// Cpp version pop haloupdate 2D
void pop_halo_update_2dr8(
    double (&array)[CppParamMod::JMT][CppParamMod::IMT],
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 2d horizontal arrays of double precision.

  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::mytid;
  using CppParamMod::IMT_GLOBAL;
  using CppPOPBlocksMod::POP_HALO_WIDTH;
  using CppPOPCommMod::POP_MPI_TAG_HALO;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_ANGLE;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_W_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_S_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER;

  // Halo info
  const MPI_Comm comm          = halo.communicator;
  const int  numMsgSend        = halo.numMsgSend;
  const int  numMsgRecv        = halo.numMsgRecv;
  const int  numLocalCopies    = halo.numLocalCopies;
  int* const recvTask          = halo.recvTask;
  int* const sendTask          = halo.sendTask;
  int* const sizeSend          = halo.sizeSend;
  int* const sizeRecv          = halo.sizeRecv;
  int const (*srcLocalAddr)[3] = halo.srcLocalAddr;
  int const (*dstLocalAddr)[3] = halo.dstLocalAddr;
  int const (*sendAddr)[3]     = halo.sendAddr;
  int const (*recvAddr)[3]     = halo.recvAddr;

  // dummy loop indices
  int i, n, nmsg; 
  //int i, j, n, nmsg; 
  // error or status flag for MPI,alloc
  //int ierr;
  // size of an individual message
  int msgSize;
  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  // current version srcBlock = 1
  int srcBlock;
  // local block number for destination
  // current version dstBlock = 1
  int dstBlock;
  // address shifts for tripole
  int ioffset(0), joffset(0);
  // sign factor for tripole grids
  int isign(0);

  const int fill(0);
  // scalars for enforcing symmetry at U pts
  double x1(0.0), x2(0.0), xavg(0.0);

	double* bufRecv2DR8    = static_cast<double*>(bufRecv);
	double* bufSend2DR8    = static_cast<double*>(bufSend);
	double* bufTripole2DR8 = static_cast<double*>(bufTripole);

  // MPI var
  // MPI request ids
  // MPI_Request* sndRequest = (MPI_Request*)malloc(numMsgSend * sizeof(MPI_Request));
  // MPI_Request* rcvRequest = (MPI_Request*)malloc(numMsgRecv * sizeof(MPI_Request));
  // MPI status flags
  // MPI_Status* sndStatus = (MPI_Status*)malloc(numMsgSend * sizeof(MPI_Status));
  // MPI_Status* rcvStatus = (MPI_Status*)malloc(numMsgRecv * sizeof(MPI_Status));

  // memset(bufTripole2DR8, fill, (POP_HALO_WIDTH + 1) * IMT_GLOBAL * sizeof(double));

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    msgSize = sizeRecv[nmsg];

    MPI_Irecv(&bufRecv2DR8[nmsg * bufSizeRecv], msgSize, MPI_DOUBLE,
        recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
            comm, &rcvRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numMsgSend; ++nmsg) {
    msgSize = sizeSend[nmsg];
    for (n = 0; n < msgSize; ++n) {
      iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
      jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
      // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];
      bufSend2DR8[nmsg * bufSizeSend + n] = array[jSrc][iSrc];
    }
    // for (n = msgSize; n < bufSizeSend; ++n) {
    //   bufSend2DR8[nmsg * bufSizeSend + n] = static_cast<double>(fill);
    // }

    MPI_Isend(&bufSend2DR8[nmsg * bufSizeSend], msgSize, MPI_DOUBLE,
        sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  do local copies while waiting for messages to complete
	// !  if srcBlock is zero, that denotes an eliminated land block or a 
	// !    closed boundary where ghost cell values are undefined
	// !  if srcBlock is less than zero, the message is a copy out of the
	// !    tripole buffer and will be treated later
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
    iSrc     = srcLocalAddr[nmsg][0] - 1;
    jSrc     = srcLocalAddr[nmsg][1] - 1;
    srcBlock = srcLocalAddr[nmsg][2];
    iDst     = dstLocalAddr[nmsg][0] - 1;
    jDst     = dstLocalAddr[nmsg][1] - 1;
    dstBlock = dstLocalAddr[nmsg][2];

    if (srcBlock > 0) {
      if (dstBlock > 0) {
        array[jDst][iDst] = array[jSrc][iSrc];
      } else if (dstBlock < 0) {
				// tripole copy into buffer
        bufTripole2DR8[jDst * nxGlobal + iDst] = array[jSrc][iSrc];
      }
    } else if (srcBlock == 0) {
      array[jDst][iDst] = static_cast<double>(fill);
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for receives to finish and then unpack the recv buffer into
	// !  ghost cells
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgRecv, rcvRequest, rcvStatus);

  for (nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    const int sizeRecv_tmp = sizeRecv[nmsg];
    for (n = 0; n < sizeRecv_tmp; ++n) {
      iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
      jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
      dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
      if (dstBlock > 0) {
        array[jDst][iDst] = bufRecv2DR8[nmsg * bufSizeRecv + n];
      } else if (dstBlock < 0) {
				// tripole
        bufTripole2DR8[jDst * nxGlobal + iDst] =  bufRecv2DR8[nmsg * bufSizeRecv + n];
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  take care of northern boundary in tripole case
	// !  bufTripole array contains the top haloWidth+1 rows of physical
	// !    domain for entire (global) top row 
	// !
	// !-----------------------------------------------------------------------

  if (nxGlobal > 0) {
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign = 1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_VECTOR) {
      isign = -1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_ANGLE) {
      isign = -1;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field kind\n");
        exit(0);
      }
    }

    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
      for (i = 0; i < (nxGlobal/2) - 1; ++i) {
        iDst = nxGlobal - i - 2;
        x1 = bufTripole2DR8[i];
        x2 = bufTripole2DR8[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2DR8[i]    = xavg;
          bufTripole2DR8[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          bufTripole2DR8[i]    = isign * sign(xavg, x2);
          bufTripole2DR8[iDst] = isign * sign(xavg, x1);
        }
      }
      bufTripole2DR8[nxGlobal - 1] *= isign;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
      for (i = 0; i < nxGlobal/2; ++i) {
        iDst = nxGlobal - i - 1;
        x1 = bufTripole2DR8[i];
        x2 = bufTripole2DR8[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2DR8[i]    = xavg;
          bufTripole2DR8[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          // bufTripole2DR8[i-1]    = isign * sign(xavg, x2);
          // bufTripole2DR8[iDst-1] = isign * sign(xavg, x1);
          bufTripole2DR8[i]    = isign * sign(xavg, x2);
          bufTripole2DR8[iDst] = isign * sign(xavg, x1);
        }
      }
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_S_FACE) {
      ioffset = 1;
      joffset = 2;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field location\n");
        exit(0);
      }
    }

    for (nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
      srcBlock = srcLocalAddr[nmsg][2];
      if (srcBlock < 0) {
        iSrc = srcLocalAddr[nmsg][0] - 1;
        jSrc = srcLocalAddr[nmsg][1] - 1;

        iDst = dstLocalAddr[nmsg][0] - 1;
        jDst = dstLocalAddr[nmsg][1] - 1;
        dstBlock = dstLocalAddr[nmsg][2];

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
          array[jDst][iDst] = isign * bufTripole2DR8[jSrc * nxGlobal + iSrc];
        }
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for sends to complete and deallocate arrays
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgSend, sndRequest, sndStatus);
  // free(sndRequest);
  // free(rcvRequest);
  // free(sndStatus);
  // free(rcvStatus);

  return;
}
void pop_halo_update_2dr8(double* const array, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 2d horizontal arrays of double precision.

  using CppParamMod::mytid;
  using CppParamMod::IMT_GLOBAL;
  using CppPOPBlocksMod::POP_HALO_WIDTH;
  using CppPOPCommMod::POP_MPI_TAG_HALO;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_ANGLE;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_W_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_S_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER;

  // Halo info
  const MPI_Comm comm          = halo.communicator;
  const int  numMsgSend        = halo.numMsgSend;
  const int  numMsgRecv        = halo.numMsgRecv;
  const int  numLocalCopies    = halo.numLocalCopies;
  int* const recvTask          = halo.recvTask;
  int* const sendTask          = halo.sendTask;
  int* const sizeSend          = halo.sizeSend;
  int* const sizeRecv          = halo.sizeRecv;
  int const (*srcLocalAddr)[3] = halo.srcLocalAddr;
  int const (*dstLocalAddr)[3] = halo.dstLocalAddr;
  int const (*sendAddr)[3]     = halo.sendAddr;
  int const (*recvAddr)[3]     = halo.recvAddr;

  // dummy loop indices
  int i, n, nmsg; 
  //int i, j, n, nmsg; 
  // error or status flag for MPI,alloc
  //int ierr;
  // size of an individual message
  int msgSize;
  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  // current version srcBlock = 1
  int srcBlock;
  // local block number for destination
  // current version dstBlock = 1
  int dstBlock;
  // address shifts for tripole
  int ioffset(0), joffset(0);
  // sign factor for tripole grids
  int isign(0);

  const int fill(0);
  // scalars for enforcing symmetry at U pts
  double x1(0.0), x2(0.0), xavg(0.0);

	double* bufRecv2DR8    = static_cast<double*>(bufRecv);
	double* bufSend2DR8    = static_cast<double*>(bufSend);
	double* bufTripole2DR8 = static_cast<double*>(bufTripole);

  // MPI var
  // MPI request ids
  // MPI_Request* sndRequest = (MPI_Request*)malloc(numMsgSend * sizeof(MPI_Request));
  // MPI_Request* rcvRequest = (MPI_Request*)malloc(numMsgRecv * sizeof(MPI_Request));
  // MPI status flags
  // MPI_Status* sndStatus = (MPI_Status*)malloc(numMsgSend * sizeof(MPI_Status));
  // MPI_Status* rcvStatus = (MPI_Status*)malloc(numMsgRecv * sizeof(MPI_Status));

  // memset(bufTripole2DR8, fill, (POP_HALO_WIDTH + 1) * IMT_GLOBAL * sizeof(double));

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    msgSize = sizeRecv[nmsg];

    MPI_Irecv(&bufRecv2DR8[nmsg * bufSizeRecv], msgSize, MPI_DOUBLE,
        recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
            comm, &rcvRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numMsgSend; ++nmsg) {
    msgSize = sizeSend[nmsg];
    for (n = 0; n < msgSize; ++n) {
      iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
      jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
      // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];
      bufSend2DR8[nmsg * bufSizeSend + n] = array[jSrc * I0 + iSrc];
    }
    // for (n = msgSize; n < bufSizeSend; ++n) {
    //   bufSend2DR8[nmsg * bufSizeSend + n] = static_cast<double>(fill);
    // }

    MPI_Isend(&bufSend2DR8[nmsg * bufSizeSend], msgSize, MPI_DOUBLE,
        sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  do local copies while waiting for messages to complete
	// !  if srcBlock is zero, that denotes an eliminated land block or a 
	// !    closed boundary where ghost cell values are undefined
	// !  if srcBlock is less than zero, the message is a copy out of the
	// !    tripole buffer and will be treated later
	// !
	// !-----------------------------------------------------------------------

  for (nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
    iSrc     = srcLocalAddr[nmsg][0] - 1;
    jSrc     = srcLocalAddr[nmsg][1] - 1;
    srcBlock = srcLocalAddr[nmsg][2];
    iDst     = dstLocalAddr[nmsg][0] - 1;
    jDst     = dstLocalAddr[nmsg][1] - 1;
    dstBlock = dstLocalAddr[nmsg][2];

    if (srcBlock > 0) {
      if (dstBlock > 0) {
        array[jDst * I0 + iDst] = array[jSrc * I0 + iSrc];
      } else if (dstBlock < 0) {
				// tripole copy into buffer
        bufTripole2DR8[jDst * nxGlobal + iDst] = array[jSrc * I0 + iSrc];
      }
    } else if (srcBlock == 0) {
      array[jDst * I0 + iDst] = static_cast<double>(fill);
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for receives to finish and then unpack the recv buffer into
	// !  ghost cells
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgRecv, rcvRequest, rcvStatus);

  for (nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    const int sizeRecv_tmp = sizeRecv[nmsg];
    for (n = 0; n < sizeRecv_tmp; ++n) {
      iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
      jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
      dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
      if (dstBlock > 0) {
        array[jDst * I0 + iDst] = bufRecv2DR8[nmsg * bufSizeRecv + n];
      } else if (dstBlock < 0) {
				// tripole
        bufTripole2DR8[jDst * nxGlobal + iDst] =  bufRecv2DR8[nmsg * bufSizeRecv + n];
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  take care of northern boundary in tripole case
	// !  bufTripole array contains the top haloWidth+1 rows of physical
	// !    domain for entire (global) top row 
	// !
	// !-----------------------------------------------------------------------

  if (nxGlobal > 0) {
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign = 1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_VECTOR) {
      isign = -1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_ANGLE) {
      isign = -1;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field kind\n");
        exit(0);
      }
    }

    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
      for (i = 0; i < (nxGlobal/2) - 1; ++i) {
        iDst = nxGlobal - i - 2;
        x1 = bufTripole2DR8[i];
        x2 = bufTripole2DR8[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2DR8[i]    = xavg;
          bufTripole2DR8[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          bufTripole2DR8[i]    = isign * sign(xavg, x2);
          bufTripole2DR8[iDst] = isign * sign(xavg, x1);
        }
      }
      bufTripole2DR8[nxGlobal - 1] *= isign;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
      for (i = 0; i < nxGlobal/2; ++i) {
        iDst = nxGlobal - i - 1;
        x1 = bufTripole2DR8[i];
        x2 = bufTripole2DR8[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2DR8[i]    = xavg;
          bufTripole2DR8[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          // bufTripole2DR8[i-1]    = isign * sign(xavg, x2);
          // bufTripole2DR8[iDst-1] = isign * sign(xavg, x1);
          bufTripole2DR8[i]    = isign * sign(xavg, x2);
          bufTripole2DR8[iDst] = isign * sign(xavg, x1);
        }
      }
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_S_FACE) {
      ioffset = 1;
      joffset = 2;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field location\n");
        exit(0);
      }
    }

    for (nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
      srcBlock = srcLocalAddr[nmsg][2];
      if (srcBlock < 0) {
        iSrc = srcLocalAddr[nmsg][0] - 1;
        jSrc = srcLocalAddr[nmsg][1] - 1;

        iDst = dstLocalAddr[nmsg][0] - 1;
        jDst = dstLocalAddr[nmsg][1] - 1;
        // dstBlock = dstLocalAddr[nmsg][2];

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
          array[jDst * I0 + iDst] = isign * bufTripole2DR8[jSrc * nxGlobal + iSrc];
        }
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for sends to complete and deallocate arrays
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgSend, sndRequest, sndStatus);
  // free(sndRequest);
  // free(rcvRequest);
  // free(sndStatus);
  // free(rcvStatus);

  return;
}
// End cpp version pop haloupdate 2D

// Cpp version pop haloupdate 3D
void pop_halo_update_3dr8(
    double (&array)[CppParamMod::KM][CppParamMod::JMT][CppParamMod::IMT],
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 3d horizontal arrays of double precision.

  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::KM;
  using CppParamMod::mytid;
  using CppParamMod::IMT_GLOBAL;
  using CppPOPBlocksMod::POP_HALO_WIDTH;
  using CppPOPCommMod::POP_MPI_TAG_HALO;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_ANGLE;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_W_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_S_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER;

  // Halo info
  const MPI_Comm comm          = halo.communicator;
  const int  numMsgSend        = halo.numMsgSend;
  const int  numMsgRecv        = halo.numMsgRecv;
  const int  numLocalCopies    = halo.numLocalCopies;
  int* const recvTask          = halo.recvTask;
  int* const sendTask          = halo.sendTask;
  int* const sizeSend          = halo.sizeSend;
  int* const sizeRecv          = halo.sizeRecv;
  int const (*srcLocalAddr)[3] = halo.srcLocalAddr;
  int const (*dstLocalAddr)[3] = halo.dstLocalAddr;
  int const (*sendAddr)[3]     = halo.sendAddr;
  int const (*recvAddr)[3]     = halo.recvAddr;

  // error or status flag for MPI,alloc
  //int ierr;
  // size of an individual message
  int msgSize;
  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  // current version srcBlock = 1
  int srcBlock;
  // local block number for destination
  // current version dstBlock = 1
  int dstBlock;
  // address shifts for tripole
  int ioffset(0), joffset(0);
  // sign factor for tripole grids
  int isign(0);

  const int fill(0);
  // scalars for enforcing symmetry at U pts
  double x1(0.0), x2(0.0), xavg(0.0);
	const int strideSend = bufSizeSend * KM;
	const int strideRecv = bufSizeRecv * KM;
	const int strideBuf  = (POP_HALO_WIDTH + 1) * nxGlobal;

	double* bufRecv3DR8    = static_cast<double*>(bufRecv);
	double* bufSend3DR8    = static_cast<double*>(bufSend);
	double* bufTripole3DR8 = static_cast<double*>(bufTripole);
  // MPI var
  // MPI request ids
  // MPI_Request* sndRequest = (MPI_Request*)malloc(numMsgSend * sizeof(MPI_Request));
  // MPI_Request* rcvRequest = (MPI_Request*)malloc(numMsgRecv * sizeof(MPI_Request));
  // MPI status flags
  // MPI_Status* sndStatus = (MPI_Status*)malloc(numMsgSend * sizeof(MPI_Status));
  // MPI_Status* rcvStatus = (MPI_Status*)malloc(numMsgRecv * sizeof(MPI_Status));

  // memset(bufTripole3DR8, fill, KM * (POP_HALO_WIDTH + 1) * nxGlobal * sizeof(double));

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    msgSize = sizeRecv[nmsg] * KM;
    MPI_Irecv(&bufRecv3DR8[nmsg * strideRecv], msgSize, MPI_DOUBLE,
        recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg], comm, &rcvRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
		int i = 0;
		const int size_Send = sizeSend[nmsg];
    for (int n = 0; n < size_Send; ++n) {
      iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
      jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
      // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];

			for (int k = 0; k < KM; ++k) {
      	bufSend3DR8[nmsg * strideSend + i] = array[k][jSrc][iSrc];
				i += 1;
			}
    }
    msgSize = sizeSend[nmsg] * KM;
    // for (int n = i; n < msgSize; ++n) {
    //   bufSend3DR8[nmsg * strideSend + n] = static_cast<double>(fill);
    // }

    MPI_Isend(&bufSend3DR8[nmsg * strideSend], msgSize, MPI_DOUBLE,
        sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  do local copies while waiting for messages to complete
	// !  if srcBlock is zero, that denotes an eliminated land block or a 
	// !    closed boundary where ghost cell values are undefined
	// !  if srcBlock is less than zero, the message is a copy out of the
	// !    tripole buffer and will be treated later
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
    iSrc     = srcLocalAddr[nmsg][0] - 1;
    jSrc     = srcLocalAddr[nmsg][1] - 1;
    srcBlock = srcLocalAddr[nmsg][2];
    iDst     = dstLocalAddr[nmsg][0] - 1;
    jDst     = dstLocalAddr[nmsg][1] - 1;
    dstBlock = dstLocalAddr[nmsg][2];

    if (srcBlock > 0) {
      if (dstBlock > 0) {
				for (int k = 0; k < KM; ++k) {
        	array[k][jDst][iDst] = array[k][jSrc][iSrc];
				}
      } else if (dstBlock < 0) {
				// tripole copy into buffer
				for (int k = 0; k < KM; ++k) {
        	bufTripole3DR8[k * strideBuf + jDst * nxGlobal + iDst] = array[k][jSrc][iSrc];
				}
      }
    } else if (srcBlock == 0) {
			for (int k = 0; k < KM; ++k) {
      	array[k][jDst][iDst] = static_cast<double>(fill);
			}
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for receives to finish and then unpack the recv buffer into
	// !  ghost cells
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgRecv, rcvRequest, rcvStatus);

  for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    const int sizeRecv_tmp = sizeRecv[nmsg];
		int i = 0;
    for (int n = 0; n < sizeRecv_tmp; ++n) {
      iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
      jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
      dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
      if (dstBlock > 0) {
				for (int k = 0; k < KM; ++k) {
        	array[k][jDst][iDst] = bufRecv3DR8[nmsg * strideRecv + i];
					i += 1;
				}
      } else if (dstBlock < 0) {
				// tripole
				for (int k = 0; k < KM; ++k) {
					bufTripole3DR8[k * strideBuf + jDst * nxGlobal + iDst] =  bufRecv3DR8[nmsg * strideRecv + i];
					i += 1;
				}
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  take care of northern boundary in tripole case
	// !  bufTripole array contains the top haloWidth+1 rows of physical
	// !    domain for entire (global) top row 
	// !
	// !-----------------------------------------------------------------------

  if (nxGlobal > 0) {
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign = 1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_VECTOR) {
      isign = -1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_ANGLE) {
      isign = -1;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field kind\n");
        exit(0);
      }
    }

    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
			const int iend = (nxGlobal / 2) - 1;
			for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < iend; ++i) {
          iDst = nxGlobal - i - 2;
          x1 = bufTripole3DR8[k * strideBuf + i];
          x2 = bufTripole3DR8[k * strideBuf + iDst];
          if (isign == 1) {
            xavg = 0.5 * (x1 + x2);
            bufTripole3DR8[k * strideBuf + i]    = xavg;
            bufTripole3DR8[k * strideBuf + iDst] = xavg;
          } else {
            xavg = 0.5 * (abs(x1) + abs(x2));
            bufTripole3DR8[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3DR8[k * strideBuf + iDst] = isign * sign(xavg, x1);
          }
        }
        bufTripole3DR8[k * strideBuf + nxGlobal - 1] *= isign;
			}
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
			const int iend = nxGlobal / 2;
			for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < iend; ++i) {
          iDst = nxGlobal - i - 1;
          x1 = bufTripole3DR8[k * strideBuf + i];
          x2 = bufTripole3DR8[k * strideBuf + iDst];
          if (isign == 1) {
            xavg = 0.5 * (x1 + x2);
            bufTripole3DR8[k * strideBuf + i]    = xavg;
            bufTripole3DR8[k * strideBuf + iDst] = xavg;
          } else {
            xavg = 0.5 * (abs(x1) + abs(x2));
            // bufTripole3DR8[i-1]    = isign * sign(xavg, x2);
            // bufTripole3DR8[iDst-1] = isign * sign(xavg, x1);
            bufTripole3DR8[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3DR8[k * strideBuf + iDst] = isign * sign(xavg, x1);
          }
        }
			}
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_S_FACE) {
      ioffset = 1;
      joffset = 2;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field location\n");
        exit(0);
      }
    }

    for (int nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
      srcBlock = srcLocalAddr[nmsg][2];
      if (srcBlock < 0) {
        iSrc = srcLocalAddr[nmsg][0] - 1;
        jSrc = srcLocalAddr[nmsg][1] - 1;

        iDst = dstLocalAddr[nmsg][0] - 1;
        jDst = dstLocalAddr[nmsg][1] - 1;
        // dstBlock = dstLocalAddr[nmsg][2];

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
					for (int k = 0; k < KM; ++k) {
          	array[k][jDst][iDst] = isign 
								* bufTripole3DR8[k * strideBuf + jSrc * nxGlobal + iSrc];
					}
        }
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for sends to complete and deallocate arrays
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgSend, sndRequest, sndStatus);
  // free(sndRequest);
  // free(rcvRequest);
  // free(sndStatus);
  // free(rcvStatus);

  return;
}

void pop_halo_update_3dr8(double* const array, 
		const int &I2, const int &I1, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 3d horizontal arrays of double precision.

  using CppParamMod::mytid;
  using CppParamMod::IMT_GLOBAL;
  using CppPOPBlocksMod::POP_HALO_WIDTH;
  using CppPOPCommMod::POP_MPI_TAG_HALO;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_ANGLE;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR;
  using CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_W_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_S_FACE;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER;
  using CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER;

  // Halo info
  const MPI_Comm comm          = halo.communicator;
  const int  numMsgSend        = halo.numMsgSend;
  const int  numMsgRecv        = halo.numMsgRecv;
  const int  numLocalCopies    = halo.numLocalCopies;
  int* const recvTask          = halo.recvTask;
  int* const sendTask          = halo.sendTask;
  int* const sizeSend          = halo.sizeSend;
  int* const sizeRecv          = halo.sizeRecv;
  int const (*srcLocalAddr)[3] = halo.srcLocalAddr;
  int const (*dstLocalAddr)[3] = halo.dstLocalAddr;
  int const (*sendAddr)[3]     = halo.sendAddr;
  int const (*recvAddr)[3]     = halo.recvAddr;

  // error or status flag for MPI,alloc
  //int ierr;
  // size of an individual message
  int msgSize;
  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  // current version srcBlock = 1
  int srcBlock;
  // local block number for destination
  // current version dstBlock = 1
  int dstBlock;
  // address shifts for tripole
  int ioffset(0), joffset(0);
  // sign factor for tripole grids
  int isign(0);

  const int fill(0);
  // scalars for enforcing symmetry at U pts
  double x1(0.0), x2(0.0), xavg(0.0);
	const int strideSend = bufSizeSend * I2;
	const int strideRecv = bufSizeRecv * I2;
	const int strideBuf  = (POP_HALO_WIDTH + 1) * nxGlobal;
	const int strideArr  = I0 * I1;

	double* bufRecv3DR8    = static_cast<double*>(bufRecv);
	double* bufSend3DR8    = static_cast<double*>(bufSend);
	double* bufTripole3DR8 = static_cast<double*>(bufTripole);

  // MPI var
  // MPI request ids
  // MPI_Request* sndRequest = (MPI_Request*)malloc(numMsgSend * sizeof(MPI_Request));
  // MPI_Request* rcvRequest = (MPI_Request*)malloc(numMsgRecv * sizeof(MPI_Request));
  // MPI status flags
  // MPI_Status* sndStatus = (MPI_Status*)malloc(numMsgSend * sizeof(MPI_Status));
  // MPI_Status* rcvStatus = (MPI_Status*)malloc(numMsgRecv * sizeof(MPI_Status));

  // memset(bufTripole3DR8, fill, KM * (POP_HALO_WIDTH + 1) * nxGlobal * sizeof(double));

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    msgSize = sizeRecv[nmsg] * I2;
    MPI_Irecv(&bufRecv3DR8[nmsg * strideRecv], msgSize, MPI_DOUBLE,
        recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg], comm, &rcvRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
		int i = 0;
		const int size_Send = sizeSend[nmsg];
    for (int n = 0; n < size_Send; ++n) {
      iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
      jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
      // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];

			for (int k = 0; k < I2; ++k) {
      	bufSend3DR8[nmsg * strideSend + i] = array[k*strideArr + jSrc*I0 + iSrc];
				i += 1;
			}
    }
    msgSize = sizeSend[nmsg] * I2;
    // for (int n = i; n < msgSize; ++n) {
    //   bufSend3DR8[nmsg * strideSend + n] = static_cast<double>(fill);
    // }

    MPI_Isend(&bufSend3DR8[nmsg * strideSend], msgSize, MPI_DOUBLE,
        sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
  }

	// !-----------------------------------------------------------------------
	// !
	// !  do local copies while waiting for messages to complete
	// !  if srcBlock is zero, that denotes an eliminated land block or a 
	// !    closed boundary where ghost cell values are undefined
	// !  if srcBlock is less than zero, the message is a copy out of the
	// !    tripole buffer and will be treated later
	// !
	// !-----------------------------------------------------------------------

  for (int nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
    iSrc     = srcLocalAddr[nmsg][0] - 1;
    jSrc     = srcLocalAddr[nmsg][1] - 1;
    srcBlock = srcLocalAddr[nmsg][2];
    iDst     = dstLocalAddr[nmsg][0] - 1;
    jDst     = dstLocalAddr[nmsg][1] - 1;
    dstBlock = dstLocalAddr[nmsg][2];

    if (srcBlock > 0) {
      if (dstBlock > 0) {
				for (int k = 0; k < I2; ++k) {
        	array[k*strideArr + jDst*I0 + iDst] = array[k*strideArr + jSrc*I0 + iSrc];
				}
      } else if (dstBlock < 0) {
				// tripole copy into buffer
				for (int k = 0; k < I2; ++k) {
        	bufTripole3DR8[k * strideBuf + jDst * nxGlobal + iDst] = array[k*strideArr + jSrc*I0 + iSrc];
				}
      }
    } else if (srcBlock == 0) {
			for (int k = 0; k < I2; ++k) {
      	array[k*strideArr + jDst*I0 + iDst] = static_cast<double>(fill);
			}
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for receives to finish and then unpack the recv buffer into
	// !  ghost cells
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgRecv, rcvRequest, rcvStatus);

  for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
    const int sizeRecv_tmp = sizeRecv[nmsg];
		int i = 0;
    for (int n = 0; n < sizeRecv_tmp; ++n) {
      iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
      jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
      dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
      if (dstBlock > 0) {
				for (int k = 0; k < I2; ++k) {
        	array[k*strideArr + jDst*I0 + iDst] = bufRecv3DR8[nmsg * strideRecv + i];
					i += 1;
				}
      } else if (dstBlock < 0) {
				// tripole
				for (int k = 0; k < I2; ++k) {
					bufTripole3DR8[k * strideBuf + jDst * nxGlobal + iDst] =  bufRecv3DR8[nmsg * strideRecv + i];
					i += 1;
				}
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  take care of northern boundary in tripole case
	// !  bufTripole array contains the top haloWidth+1 rows of physical
	// !    domain for entire (global) top row 
	// !
	// !-----------------------------------------------------------------------

  if (nxGlobal > 0) {
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign = 1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_VECTOR) {
      isign = -1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_ANGLE) {
      isign = -1;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field kind\n");
        exit(0);
      }
    }

    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
			const int iend = (nxGlobal / 2) - 1;
			for (int k = 0; k < I2; ++k) {
        for (int i = 0; i < iend; ++i) {
          iDst = nxGlobal - i - 2;
          x1 = bufTripole3DR8[k * strideBuf + i];
          x2 = bufTripole3DR8[k * strideBuf + iDst];
          if (isign == 1) {
            xavg = 0.5 * (x1 + x2);
            bufTripole3DR8[k * strideBuf + i]    = xavg;
            bufTripole3DR8[k * strideBuf + iDst] = xavg;
          } else {
            xavg = 0.5 * (abs(x1) + abs(x2));
            bufTripole3DR8[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3DR8[k * strideBuf + iDst] = isign * sign(xavg, x1);
          }
        }
        bufTripole3DR8[k * strideBuf + nxGlobal - 1] *= isign;
			}
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
			const int iend = nxGlobal / 2;
			for (int k = 0; k < I2; ++k) {
        for (int i = 0; i < iend; ++i) {
          iDst = nxGlobal - i - 1;
          x1 = bufTripole3DR8[k * strideBuf + i];
          x2 = bufTripole3DR8[k * strideBuf + iDst];
          if (isign == 1) {
            xavg = 0.5 * (x1 + x2);
            bufTripole3DR8[k * strideBuf + i]    = xavg;
            bufTripole3DR8[k * strideBuf + iDst] = xavg;
          } else {
            xavg = 0.5 * (abs(x1) + abs(x2));
            // bufTripole3DR8[i-1]    = isign * sign(xavg, x2);
            // bufTripole3DR8[iDst-1] = isign * sign(xavg, x1);
            bufTripole3DR8[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3DR8[k * strideBuf + iDst] = isign * sign(xavg, x1);
          }
        }
			}
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_S_FACE) {
      ioffset = 1;
      joffset = 2;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate2DR8: Unknown field location\n");
        exit(0);
      }
    }

    for (int nmsg = 0; nmsg < numLocalCopies; ++nmsg) {
      srcBlock = srcLocalAddr[nmsg][2];
      if (srcBlock < 0) {
        iSrc = srcLocalAddr[nmsg][0] - 1;
        jSrc = srcLocalAddr[nmsg][1] - 1;

        iDst = dstLocalAddr[nmsg][0] - 1;
        jDst = dstLocalAddr[nmsg][1] - 1;
        // dstBlock = dstLocalAddr[nmsg][2];

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
					for (int k = 0; k < I2; ++k) {
          	array[k*strideArr + jDst*I0 + iDst] = isign 
								* bufTripole3DR8[k * strideBuf + jSrc * nxGlobal + iSrc];
					}
        }
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  wait for sends to complete and deallocate arrays
	// !
	// !-----------------------------------------------------------------------

  MPI_Waitall(numMsgSend, sndRequest, sndStatus);
  // free(sndRequest);
  // free(rcvRequest);
  // free(sndStatus);
  // free(rcvStatus);

  return;
}

} // namespace CppPOPHaloMod

#endif // LICOM_ENABLE_FORTRAN
