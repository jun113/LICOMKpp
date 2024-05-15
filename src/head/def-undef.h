#ifndef LICOM3_KOKKOS_SRC_HEAD_DEF_UNDEF_H_
#define LICOM3_KOKKOS_SRC_HEAD_DEF_UNDEF_H_

#define SPMD
#define SYNCH
#undef  FRC_ANN
#define CDFIN
#undef  FRC_DAILY
#define FRC_CORE
#define SOLAR
#define ACOS
#define BIHAR
#undef  SMAG_FZ
#undef  SMAG_OUT
#define NETCDF
#undef  BOUNDARY
#define NODIAG
#undef  ICE
#undef  SHOW_TIME
#undef  DEBUG
#undef  COUP
#undef  ISO
#define D_PRECISION
#define CANUTO
#undef  SOLARCHLORO
#undef  BCKMEX
#define TIDEMIX
#define SSSNORM
#undef  LDD97
#undef  SMAG
#define TEST_IDEAL
#define LOWRES 
#undef  HIGHRES 
#undef  SUPHIGH 
#undef  DAILYACC
#undef  DAILYBUGDET
#define BLCKX   180
#define BLCKY	110
#define NJMT	218 
#define NIMT	360
#define NKM     30 
#define MXBLCKS 1

#define LICOM_RES "100km"

#undef LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_TIME

#undef LICOM_ENABLE_GPTL

#endif

#define LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#define LICOM_ENABLE_KOKKOS

#endif
