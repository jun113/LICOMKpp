#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_POP_DISTRIBUTION_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_POP_DISTRIBUTION_MOD_H_
#include "def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

namespace CppPOPDistributionMod {

struct POP_distrb {
	// number of processors in this dist
  int numProcs;
	// communicator to use in this dist
  int communicator;
	// number of blocks distributed to this local processor
  int numLocalBlocks;

  // processor location for all blocks
	// blockLocation[POP_numBlocks]
	int *blockLocation = nullptr;
	// local block id for all blocks
	// blockLocalID[POP_numBlocks]
	int *blockLocalID  = nullptr;
	// global block id for each local block
	// blockGlobalID[Distrb.numLocalBlocks]
	int *blockGlobalID = nullptr;
};

constexpr int POP_distribMethodNone       = 0;
constexpr int POP_distribMethodCartesian  = 1;
constexpr int POP_distribMethodRake       = 2;
constexpr int POP_distribMethodSpacecurve = 3;

extern void pop_distribution_create(CppPOPDistributionMod::POP_distrb &newDistrb);

extern void pop_distribution_get(const CppPOPDistributionMod::POP_distrb &, int &, int &);

extern void pop_distribution_get_block_loc(
		const CppPOPDistributionMod::POP_distrb &,
				const int &, int &, int &);

} // namespace CppPOPDistributionMod

#endif // LICOM_ENABLE_FORTRAN

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_POP_DISTRIBUTION_MOD_H_