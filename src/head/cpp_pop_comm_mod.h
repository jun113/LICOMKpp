#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_POP_COMM_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_POP_COMM_MOD_H_
#include "def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

namespace CppPOPCommMod {

constexpr int POP_MPI_TAG_HALO   = 1;
constexpr int POP_MPI_TAG_REDIST = 1000;

extern int &POP_communicator;
extern int &POP_myTask;
extern int &POP_masterTask;

} // namespace CppPOPCommMod

#endif // LICOM_ENABLE_FORTRAN

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_POP_COMM_MOD_H_