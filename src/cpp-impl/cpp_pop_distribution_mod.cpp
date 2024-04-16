#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_pop_blocks_mod.h"
#include "../head/cpp_pop_distribution_mod.h"
#include "../head/fortran_extern_functions.h"
#include "../head/cpp_param_mod.h"

#include <cstdio>
#include <cstdlib>
namespace CppPOPDistributionMod {

void pop_distribution_create(CppPOPDistributionMod::POP_distrb &newDistrb) {

	get_f_pop_distrb_info_(&(newDistrb.numProcs), &(newDistrb.communicator),
			&(newDistrb.numLocalBlocks));

	if (newDistrb.blockLocation == nullptr) {
		// newDistrb.blockLocation = new int[CppPOPBlocksMod::POP_numBlocks];
		newDistrb.blockLocation = CppDomain::POP_distrbClinic_blockLocation;
	}
	if (newDistrb.blockLocalID == nullptr) {
		// newDistrb.blockLocalID = new int[CppPOPBlocksMod::POP_numBlocks];
		newDistrb.blockLocalID = CppDomain::POP_distrbClinic_blockLocalID;
	}
	if (newDistrb.blockGlobalID == nullptr) {
		// newDistrb.blockGlobalID = new int[newDistrb.numLocalBlocks];
		newDistrb.blockGlobalID = CppDomain::POP_distrbClinic_blockGlobalID;
	}

	if (newDistrb.blockLocation == nullptr || newDistrb.blockLocalID  == nullptr ||
			newDistrb.blockGlobalID == nullptr) {
		printf ("error in pop_distribution_create\n");
		exit(0);
	}
	return ;
}

void pop_distribution_get (const CppPOPDistributionMod::POP_distrb &distribution,
		int &numProcs, int &communicator) {
	numProcs     = distribution.numProcs;
	communicator = distribution.communicator;
	return ;
}

// IROUTINE: POP_DistributionGetBlockLoc
void pop_distribution_get_block_loc(
		const CppPOPDistributionMod::POP_distrb &distribution,
				const int &blockId, int &processor, int &localId) {
	// DESCRIPTION:
	// Given a distribution of blocks and a global block ID, return
	// the processor and local index for the block.  A zero for both
	// is returned in the case that the block has been eliminated from
	// the distribution (i.e. has no active points).
	// ------------------------------------------------------------
	// extract the location from the distribution data structure
	// ------------------------------------------------------------
	processor = distribution.blockLocation[blockId - 1];
	localId   = distribution.blockLocalID[blockId - 1];
	return;
}

} // namespace CppPOPDistributionMod

#endif // LICOM_ENABLE_FORTRAN