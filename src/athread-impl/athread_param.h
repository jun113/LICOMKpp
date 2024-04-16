#ifndef LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_
#define LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_
#if (defined __sw_host__) || (defined __sw_slave__)

struct HaloTransposeDouble3D {
  int     startB[4];
  int     endB[4];
  int     startC[4];
  int     endC[4];
  double* arrSrc;
  double* arrObj;
  int     lenA;
  int     lenB;
  int     lenC;
};

#endif // (defined __sw_host__) || (defined __sw_slave__)
#endif // LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_