#ifndef LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_POP_REDUCTIONS_MOD_HPP_
#define LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_POP_REDUCTIONS_MOD_HPP_
#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_pop_comm_mod.h"
#include "../head/cpp_pop_distribution_mod.h"

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <type_traits>

namespace CppPOPReductionsMod {

template<typename T>
void pop_global_maxval_scalar(const T &scalar,
		const CppPOPDistributionMod::POP_distrb &dist, T &globalMaxval) {
	// DESCRIPTION:
	// Computes the global maximum value of a scalar value across
	// a distributed machine.
	// 
	// REVISION HISTORY:
	// same as module
	// 
	// REMARKS:
	// This is actually the specific interface for the generic 
	// POP_GlobalMaxval function corresponding to double precision scalars.

	// number of processor participating
	int numProcs;
	// communicator for this distribution
	int communicator;

	// ----------------------------------------------
	// no operations required for serial execution - return input value
	// ----------------------------------------------
	CppPOPDistributionMod::pop_distribution_get(dist, numProcs, communicator);

  const MPI_Comm comm = MPI_Comm_f2c(communicator);
	
	// ----------------------------------------------------
	// now use MPI global reduction to reduce local maxval to global maxval
	// ----------------------------------------------------
	
	if (CppPOPCommMod::POP_myTask < numProcs) {
		if (std::is_same<typename std::decay<T>::type, int>::value) {
			MPI_Allreduce(&scalar, &globalMaxval, 1, MPI_INT, MPI_MAX, comm);
		} else if (std::is_same<typename std::decay<T>::type, float>::value) {
			MPI_Allreduce(&scalar, &globalMaxval, 1, MPI_FLOAT, MPI_MAX, comm);
		} else if (std::is_same<typename std::decay<T>::type, double>::value) {
			MPI_Allreduce(&scalar, &globalMaxval, 1, MPI_DOUBLE, MPI_MAX, comm);
		} else {
			printf("pop_global_maxval_scalar error data type\n");
			exit(0);
		}
	}
  return ;
}

void pop_global_maxval_scalar_r8(const double &scalar,
		const CppPOPDistributionMod::POP_distrb &dist,
				double &globalMaxval) {
	// DESCRIPTION:
	// Computes the global maximum value of a scalar value across
	// a distributed machine.
	// 
	// REVISION HISTORY:
	// same as module
	// 
	// REMARKS:
	// This is actually the specific interface for the generic 
	// POP_GlobalMaxval function corresponding to double precision scalars.

	// number of processor participating
	int numProcs;
	// communicator for this distribution
	int communicator;

	// ----------------------------------------------
	// no operations required for serial execution - return input value
	// ----------------------------------------------
	CppPOPDistributionMod::pop_distribution_get(dist,
			numProcs, communicator);

  const MPI_Comm comm = MPI_Comm_f2c(communicator);
	
	// ----------------------------------------------------
	// now use MPI global reduction to reduce local maxval to global maxval
	// ----------------------------------------------------
	if (CppPOPCommMod::POP_myTask < numProcs) {
		MPI_Allreduce(&scalar, &globalMaxval, 1, MPI_DOUBLE, MPI_MAX, comm);
	}
  return ;
}

} // namespace CppPOPReductionsMod

#endif // LICOM_ENABLE_FORTRAN
#endif // LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_POP_REDUCTIONS_MOD_HPP_
