#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_POP_HALO_MOD_HPP_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_POP_HALO_MOD_HPP_

#include "def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "cpp_param_mod.h"
#include "cpp_pop_blocks_mod.h"
#include "cpp_pop_comm_mod.h"
#include "cpp_pop_grid_horz_mod.h"
#include "cpp_pop_distribution_mod.h"

#include <mpi.h>

#include <cstdlib>
#include <type_traits>

namespace CppPOPHaloMod {
struct POPHalo {
  // communicator to use for update messages
  MPI_Comm communicator;
  // number of messages to send halo update
  int numMsgSend;
  // number of messages to recv halo update
  int numMsgRecv;
  // num local copies for halo update
  int numLocalCopies;
  // task from which to recv each msg
  // recvTask[numMsgRecv]
  int *recvTask = nullptr;
  // task from which to send each msg
  // sendTask[numMsgSend]
  int *sendTask = nullptr;
  // size of each sent message
  // sizeSend[numMsgSend]
  int *sizeSend = nullptr;
  // size of each recv message
  // sizeRecv[numMsgSend]
  int *sizeRecv = nullptr;
  // src addresses for each local copy
  // srcLocalAddr[numLocalCopies][3]
  int (*srcLocalAddr)[3] = nullptr;
  // dst addresses for each local copy
  // dstLocalAddr[numLocalCopies][3]
  int (*dstLocalAddr)[3] = nullptr;
  // src addresses for each sent message
  // sendAddr[numMsgSend][bufSizeSend][3]
  int (*sendAddr)[3] = nullptr;
  // dst addresses for each recvd message
  // recvAddr[numMsgRecv][bufSizeRecv][3]
  int (*recvAddr)[3] = nullptr;
}; // struct POPHalo
struct pop_halo {
  // communicator to use for update messages
  int communicator;
  // number of messages to send halo update
  int numMsgSend;
  // number of messages to recv halo update
  int numMsgRecv;
  // num local copies for halo update
  int numLocalCopies;
  // task from which to recv each msg
  // recvTask[numMsgRecv]
  int *recvTask = nullptr;
  // task from which to send each msg
  // sendTask[numMsgSend]
  int *sendTask = nullptr;
  // size of each sent message
  // sizeSend[numMsgSend]
  int *sizeSend = nullptr;
  // size of each recv message
  // sizeRecv[numMsgSend]
  int *sizeRecv = nullptr;
  // src addresses for each local copy
  // srcLocalAddr[numLocalCopies][3]
  int (*srcLocalAddr)[3] = nullptr;
  // dst addresses for each local copy
  // dstLocalAddr[numLocalCopies][3]
  int (*dstLocalAddr)[3] = nullptr;
  // src addresses for each sent message
  // sendAddr[numMsgSend][bufSizeSend][3]
  int (*sendAddr)[3] = nullptr;
  // dst addresses for each recvd message
  // recvAddr[numMsgRecv][bufSizeRecv][3]
  int (*recvAddr)[3] = nullptr;
}; // struct pop_halo

extern double* arrCommPriorK;

extern int &bufSizeSend;
extern int &bufSizeRecv;

// extern int*    &bufSendI4;
// extern int*    &bufRecvI4;
                     
// extern float*  &bufSendR4;
// extern float*  &bufRecvR4;
                    
// extern double* bufSendR8;
// extern double* bufRecvR8;

// extern int*    &bufTripoleI4;
// extern float*  &bufTripoleR4;
// extern double* bufTripoleR8;

extern void* bufSend;
extern void* bufRecv;
extern void* bufTripole;

// MPI request ids
extern MPI_Request* sndRequest;
extern MPI_Request* rcvRequest;
// MPI status flags
extern MPI_Status* sndStatus;
extern MPI_Status* rcvStatus;

// extern int* (&POPHalo_bufRecvTask);
// extern int* (&POPHalo_bufSendTask);
// extern int* (&POPHalo_bufSizeSend);
// extern int* (&POPHalo_bufSizeRecv);
// extern int* (&POPHalo_bufSrcLocalAddr);
// extern int* (&POPHalo_bufDstLocalAddr);
// extern int* (&POPHalo_bufRecvAddr);
// extern int* (&POPHalo_bufSendAddr);

extern void pop_halo_create(const CppPOPHaloMod::pop_halo &, CppPOPHaloMod::pop_halo &);

extern void pop_halo_create (const CppPOPDistributionMod::POP_distrb &distrb,
		const char *nsBoundaryType, const char *ewBoundaryType,
				const int &nxGlobal, CppPOPHaloMod::pop_halo &halo);
  
extern void pop_halo_create_from_fortran(CppPOPHaloMod::POPHalo &halo);

extern void pop_halo_update_2dr8(
    double (&array)[CppParamMod::JMT][CppParamMod::IMT],
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind);

extern void pop_halo_update_2dr8(
    double* const array, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind);

extern void pop_halo_update_3dr8(
    double (&array)[CppParamMod::KM][CppParamMod::JMT][CppParamMod::IMT],
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind);

extern void pop_halo_update_3dr8(double* const array, 
    const int &I2, const int &I1, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind);

template<typename T>
inline T sign (const T &x, const T &y) {
  return y >= static_cast<T>(0) ? std::abs(x) : -std::abs(x);
}

// POP_HaloUpdate2D
template<typename T>
void pop_halo_update(T* const array, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 2d horizontal arrays.

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

  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  int srcBlock;
  // local block number for destination
  int dstBlock;
  // address shifts for tripole
  int ioffset(0), joffset(0);
  // sign factor for tripole grids
  int isign(0);

  // scalars for enforcing symmetry at U pts
  T x1(0.0), x2(0.0), xavg(0.0);

	T* bufRecv2D    = static_cast<T *>(bufRecv);
	T* bufSend2D    = static_cast<T *>(bufSend);
	T* bufTripole2D = static_cast<T *>(bufTripole);

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------

	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg];
      MPI_Irecv(&bufRecv2D[nmsg * bufSizeRecv], msgSize, MPI_INT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg];
      MPI_Irecv(&bufRecv2D[nmsg * bufSizeRecv], msgSize, MPI_FLOAT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg];
      MPI_Irecv(&bufRecv2D[nmsg * bufSizeRecv], msgSize, MPI_DOUBLE,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	}

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeSend[nmsg];
      for (int n = 0; n < msgSize; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
        bufSend2D[nmsg * bufSizeSend + n] = array[jSrc * I0 + iSrc];
      }
      MPI_Isend(&bufSend2D[nmsg * bufSizeSend], msgSize, MPI_INT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeSend[nmsg];
      for (int n = 0; n < msgSize; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
        bufSend2D[nmsg * bufSizeSend + n] = array[jSrc * I0 + iSrc];
      }
      MPI_Isend(&bufSend2D[nmsg * bufSizeSend], msgSize, MPI_FLOAT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeSend[nmsg];
      for (int n = 0; n < msgSize; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
        bufSend2D[nmsg * bufSizeSend + n] = array[jSrc * I0 + iSrc];
      }
      MPI_Isend(&bufSend2D[nmsg * bufSizeSend], msgSize, MPI_DOUBLE,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
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
        array[jDst * I0 + iDst] = array[jSrc * I0 + iSrc];
      } else if (dstBlock < 0) {
				// tripole copy into buffer
        bufTripole2D[jDst * nxGlobal + iDst] = array[jSrc * I0 + iSrc];
      }
    } else if (srcBlock == 0) {
      array[jDst * I0 + iDst] = static_cast<T>(0.0);
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
    for (int n = 0; n < sizeRecv_tmp; ++n) {
      iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
      jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
      dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
      if (dstBlock > 0) {
        array[jDst * I0 + iDst] = bufRecv2D[nmsg * bufSizeRecv + n];
      } else if (dstBlock < 0) {
				// tripole
        bufTripole2D[jDst * nxGlobal + iDst] = bufRecv2D[nmsg * bufSizeRecv + n];
      }
    }
  }

	// !-----------------------------------------------------------------------
	// !
	// !  take care of northern boundary in tripole case
	// !  bufTripole array contains the top haloWidth+1 rows of physical
	// !  domain for entire (global) top row 
	// !
	// !-----------------------------------------------------------------------

  if (nxGlobal > 0) {
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign =  1;
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
      const int end_for = (nxGlobal >> 1) - 1;
      for (int i = 0; i < end_for; ++i) {
        iDst = nxGlobal - i - 2;
        x1 = bufTripole2D[i];
        x2 = bufTripole2D[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2D[i]    = xavg;
          bufTripole2D[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          bufTripole2D[i]    = isign * sign(xavg, x2);
          bufTripole2D[iDst] = isign * sign(xavg, x1);
        }
      }
      bufTripole2D[nxGlobal - 1] *= isign;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
      const int end_for = nxGlobal >> 1;
      for (int i = 0; i < end_for; ++i) {
        iDst = nxGlobal - i - 1;
        x1 = bufTripole2D[i];
        x2 = bufTripole2D[iDst];
        if (isign == 1) {
          xavg = 0.5 * (x1 + x2);
          bufTripole2D[i]    = xavg;
          bufTripole2D[iDst] = xavg;
        } else {
          xavg = 0.5 * (abs(x1) + abs(x2));
          bufTripole2D[i]    = isign * sign(xavg, x2);
          bufTripole2D[iDst] = isign * sign(xavg, x1);
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

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
          array[jDst * I0 + iDst] = isign * bufTripole2D[jSrc * nxGlobal + iSrc];
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
  return;
}

// POP_HaloUpdate3D
template<typename T>
void pop_halo_update(T* const array, 
		const int &I2, const int &I1, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 3d horizontal arrays.

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

  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  int srcBlock;
  // local block number for destination
  int dstBlock;

	const int strideSend = bufSizeSend * I2;
	const int strideRecv = bufSizeRecv * I2;
	const int strideBuf  = (POP_HALO_WIDTH + 1) * nxGlobal;
	const int strideArr  = I0 * I1;

	T* bufRecv3D    = static_cast<T *>(bufRecv);
	T* bufSend3D    = static_cast<T *>(bufSend);
	T* bufTripole3D = static_cast<T *>(bufTripole);

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------
	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_INT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_FLOAT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_DOUBLE,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	}

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg * strideSend + i] = array[k*strideArr + jSrc*I0 + iSrc];
						i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_INT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg * strideSend + i] = array[k*strideArr + jSrc*I0 + iSrc];
						i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_FLOAT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg * strideSend + i] = array[k*strideArr + jSrc*I0 + iSrc];
						i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_DOUBLE,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
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
        	bufTripole3D[k * strideBuf + jDst * nxGlobal + iDst] = array[k*strideArr + jSrc*I0 + iSrc];
				}
      }
    } else if (srcBlock == 0) {
			for (int k = 0; k < I2; ++k) {
      	array[k*strideArr + jDst*I0 + iDst] = static_cast<T>(0.0);
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
        	array[k*strideArr + jDst*I0 + iDst] = bufRecv3D[nmsg * strideRecv + i];
					i += 1;
				}
      } else if (dstBlock < 0) {
				// tripole
				for (int k = 0; k < I2; ++k) {
					bufTripole3D[k * strideBuf + jDst * nxGlobal + iDst] =  bufRecv3D[nmsg * strideRecv + i];
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
    // sign factor for tripole grids
    int isign(0);
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign =  1;
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

    // scalars for enforcing symmetry at U pts
    T x1(0.0), x2(0.0), xavg(0.0);
    // address shifts for tripole
    int ioffset(0), joffset(0);
    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
      if (isign == 1) {
			  const int iend = (nxGlobal >> 1) - 1;
			  for (int k = 0; k < I2; ++k) {
          for (int i = 0; i < iend; ++i) {
            iDst = nxGlobal - i - 2;
            x1 = bufTripole3D[k * strideBuf + i];
            x2 = bufTripole3D[k * strideBuf + iDst];
            xavg = 0.5 * (x1 + x2);
            bufTripole3D[k * strideBuf + i]    = xavg;
            bufTripole3D[k * strideBuf + iDst] = xavg;
          }
          // bufTripole3D[k * strideBuf + nxGlobal - 1] *= isign;
			  }
      } else {
			  const int iend = (nxGlobal >> 1) - 1;
			  for (int k = 0; k < I2; ++k) {
          for (int i = 0; i < iend; ++i) {
            iDst = nxGlobal - i - 2;
            x1 = bufTripole3D[k * strideBuf + i];
            x2 = bufTripole3D[k * strideBuf + iDst];
            xavg = 0.5 * (abs(x1) + abs(x2));
            bufTripole3D[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3D[k * strideBuf + iDst] = isign * sign(xavg, x1);
          }
          bufTripole3D[k * strideBuf + nxGlobal - 1] *= isign;
			  }
      }
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
      if (isign == 1) {
			  const int iend = nxGlobal >> 1;
			  for (int k = 0; k < I2; ++k) {
          for (int i = 0; i < iend; ++i) {
            iDst = nxGlobal - i - 1;
            x1 = bufTripole3D[k * strideBuf + i];
            x2 = bufTripole3D[k * strideBuf + iDst];
            xavg = 0.5 * (x1 + x2);
            bufTripole3D[k * strideBuf + i]    = xavg;
            bufTripole3D[k * strideBuf + iDst] = xavg;
          }
			  }
      } else {
			  const int iend = nxGlobal >> 1;
			  for (int k = 0; k < I2; ++k) {
          for (int i = 0; i < iend; ++i) {
            iDst = nxGlobal - i - 1;
            x1 = bufTripole3D[k * strideBuf + i];
            x2 = bufTripole3D[k * strideBuf + iDst];
            xavg = 0.5 * (abs(x1) + abs(x2));
            bufTripole3D[k * strideBuf + i]    = isign * sign(xavg, x2);
            bufTripole3D[k * strideBuf + iDst] = isign * sign(xavg, x1);
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
        iSrc  = srcLocalAddr[nmsg][0] - 1;
        jSrc  = srcLocalAddr[nmsg][1] - 1;
             
        iDst  = dstLocalAddr[nmsg][0] - 1;
        jDst  = dstLocalAddr[nmsg][1] - 1;

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
					for (int k = 0; k < I2; ++k) {
          	array[k*strideArr + jDst*I0 + iDst] = isign 
								* bufTripole3D[k * strideBuf + jSrc * nxGlobal + iSrc];
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
  return;
}

// 3D I2-priority, I1 I0 I2
template<typename T>
void pop_halo_update_priority_k (T* const array,
		const int &I2, const int &I1, const int &I0,
    const CppPOPHaloMod::POPHalo &halo,
    const int &fieldLoc, const int &fieldKind) {

	// DESCRIPTION:
	// This routine updates ghost cells for an input array and is a
	// member of a group of routines under the generic interface
	// POP\_HaloUpdate.  This routine is the specific interface
	// for 3d horizontal arrays.

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

  // global domain size in x (tripole)
  const int nxGlobal(IMT_GLOBAL);
  // source addresses for message
  int iSrc, jSrc;
  // dest   addresses for message
  int iDst, jDst;
  // local block number for source
  int srcBlock;
  // local block number for destination
  int dstBlock;

	const int strideSend = bufSizeSend * I2;
	const int strideRecv = bufSizeRecv * I2;
	const int strideBuf  = nxGlobal    * I2;
  const int strideArr  = I0          * I2;

	T* bufRecv3D    = static_cast<T *>(bufRecv);
	T* bufSend3D    = static_cast<T *>(bufSend);
	T* bufTripole3D = static_cast<T *>(bufTripole);

	// !-----------------------------------------------------------------------
	// !
	// !  post receives
	// !
	// !-----------------------------------------------------------------------
	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_INT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_FLOAT,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgRecv; ++nmsg) {
      // size of an individual message
      const int msgSize = sizeRecv[nmsg] * I2;
      MPI_Irecv(&bufRecv3D[nmsg * strideRecv], msgSize, MPI_DOUBLE,
          recvTask[nmsg], POP_MPI_TAG_HALO + recvTask[nmsg],
              comm, &rcvRequest[nmsg]);
	  }
	}

	// !-----------------------------------------------------------------------
	// !
	// !  fill send buffer and post sends
	// !
	// !-----------------------------------------------------------------------

	if (std::is_same<typename std::decay<T>::type, int>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg*strideSend + i] = array[jSrc*strideArr + iSrc*I2 + k];
					i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_INT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, float>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg * strideSend + i] = array[jSrc*strideArr + iSrc*I2 + k];
					i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_FLOAT,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
	} else if (std::is_same<typename std::decay<T>::type, double>::value) {
    for (int nmsg = 0; nmsg < numMsgSend; ++nmsg) {
			int i = 0;
			const int size_Send = sizeSend[nmsg];
      for (int n = 0; n < size_Send; ++n) {
        iSrc = sendAddr[nmsg * bufSizeSend + n][0] - 1;
        jSrc = sendAddr[nmsg * bufSizeSend + n][1] - 1;
 
				for (int k = 0; k < I2; ++k) {
        	bufSend3D[nmsg * strideSend + i] = array[jSrc*strideArr + iSrc*I2 + k];
					i += 1;
				}
      }
      // size of an individual message
      const int msgSize = sizeSend[nmsg] * I2;
 
      MPI_Isend(&bufSend3D[nmsg * strideSend], msgSize, MPI_DOUBLE,
          sendTask[nmsg], POP_MPI_TAG_HALO + mytid, comm, &sndRequest[nmsg]);
    }
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
          array[jDst*strideArr + iDst*I2 + k] = array[jSrc*strideArr + iSrc*I2 + k];
				}
      } else if (dstBlock < 0) {
				// tripole copy into buffer
				for (int k = 0; k < I2; ++k) {
        	bufTripole3D[jDst*strideBuf + iDst*I2 + k] = array[jSrc*strideArr + iSrc*I2 + k];
				}
      }
    } else if (srcBlock == 0) {
			for (int k = 0; k < I2; ++k) {
      	array[jDst*strideArr + iDst*I2 + k] = static_cast<T>(0.0);
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
        	array[jDst*strideArr + iDst*I2 + k] = bufRecv3D[nmsg * strideRecv + i];
					i += 1;
				}
      } else if (dstBlock < 0) {
				// tripole
				for (int k = 0; k < I2; ++k) {
					bufTripole3D[jDst*strideBuf + iDst*I2 + k] =  bufRecv3D[nmsg * strideRecv + i];
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
    // sign factor for tripole grids
    int isign(0);
    if (fieldKind == FLAG_POP_FIELD_KIND_SCALAR) {
      isign =  1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_VECTOR) {
      isign = -1;
    } else if (fieldKind == FLAG_POP_FIELD_KIND_ANGLE) {
      isign = -1;
    } else {
      if (mytid == 0) {
        printf("POP_HaloUpdate3D: Unknown field kind\n");
        exit(0);
      }
    }

    // scalars for enforcing symmetry at U pts
    T x1(0.0), x2(0.0), xavg(0.0);
    // address shifts for tripole
    int ioffset(0), joffset(0);
    if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_CENTER) {
			// cell center location
      ioffset = 1;
      joffset = 1;
      if (isign == 1) {
			  const int iend = (nxGlobal >> 1) - 1;
        for (int i = 0; i < iend; ++i) {
			    for (int k = 0; k < I2; ++k) {
            iDst = nxGlobal - i - 2;
            x1 = bufTripole3D[i   *I2 + k];
            x2 = bufTripole3D[iDst*I2 + k];
            xavg = 0.5 * (x1 + x2);
            bufTripole3D[i   *I2 + k] = xavg;
            bufTripole3D[iDst*I2 + k] = xavg;
          }
			  }
      } else {
			  const int iend = (nxGlobal >> 1) - 1;
        for (int i = 0; i < iend; ++i) {
			    for (int k = 0; k < I2; ++k) {
            iDst = nxGlobal - i - 2;
            x1 = bufTripole3D[i   *I2 + k];
            x2 = bufTripole3D[iDst*I2 + k];
            xavg = 0.5 * (std::abs(x1) + std::abs(x2));
            bufTripole3D[i   *I2 + k] = isign * sign(xavg, x2);
            bufTripole3D[iDst*I2 + k] = isign * sign(xavg, x1);
          }
			  }
        for (int k = 0; k < I2; ++k) {
          bufTripole3D[(nxGlobal-1)*I2 + k] *= isign;
        }
      }
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_SW_CORNER) {
      ioffset = 0;
      joffset = 2;
    } else if (fieldLoc == FLAG_POP_GRID_HORZ_LOC_W_FACE) {
      ioffset = 0;
      joffset = 1;
      if (isign == 1) {
			  const int iend = nxGlobal >> 1;
        for (int i = 0; i < iend; ++i) {
			    for (int k = 0; k < I2; ++k) {
            iDst = nxGlobal - i - 1;
            x1 = bufTripole3D[i    * I2 + k];
            x2 = bufTripole3D[iDst * I2 + k];
            xavg = 0.5 * (x1 + x2);
            bufTripole3D[i    * I2 + k] = xavg;
            bufTripole3D[iDst * I2 + k] = xavg;
          }
			  }
      } else {
			  const int iend = nxGlobal >> 1;
        for (int i = 0; i < iend; ++i) {
			    for (int k = 0; k < I2; ++k) {
            iDst = nxGlobal - i - 1;
            x1 = bufTripole3D[i   *I2 + k];
            x2 = bufTripole3D[iDst*I2 + k];
            if (isign == 1) {
              xavg = 0.5 * (x1 + x2);
              bufTripole3D[i   *I2 + k] = xavg;
              bufTripole3D[iDst*I2 + k] = xavg;
            } else {
              xavg = 0.5 * (abs(x1) + abs(x2));
              bufTripole3D[i   *I2 + k] = isign * sign(xavg, x2);
              bufTripole3D[iDst*I2 + k] = isign * sign(xavg, x1);
            }
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
        iSrc  = srcLocalAddr[nmsg][0] - 1;
        jSrc  = srcLocalAddr[nmsg][1] - 1;
             
        iDst  = dstLocalAddr[nmsg][0] - 1;
        jDst  = dstLocalAddr[nmsg][1] - 1;

        iSrc -= ioffset;
        jSrc -= joffset;
        if (iSrc == -1) {
          iSrc = nxGlobal - 1;
        }
        if (jSrc > -1) {
					for (int k = 0; k < I2; ++k) {
          	array[jDst*strideArr + iDst*I2 + k] = isign 
								* bufTripole3D[jSrc*strideBuf + iSrc*I2 + k];
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
  return;
}

} // namespace CppPOPHaloMod

#endif // LICOM_ENABLE_FORTRAN
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_POP_HALO_MOD_HPP_
