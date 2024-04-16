#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_DOMAIN_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_DOMAIN_H_

#include "def-undef.h"

#ifndef LICOM_ENABLE_FORTRAN
#include "cpp_precision_mod.h"
#include "cpp_pop_halo_mod.hpp"
#include "cpp_pop_distribution_mod.h"

namespace CppDomain {

extern CppPOPDistributionMod::POP_distrb POP_distrbClinic;
extern CppPOPHaloMod::pop_halo POP_haloClinic;
extern CppPOPHaloMod::POPHalo POP_haloClinic_C;

extern int &nblocks_clinic;
extern int &nblocks_tropic;

extern int* &blocks_clinic;
extern int* &blocks_tropic;

// extern CppPOPDistributionMod::POP_distrb &POP_distrbClinic;
// extern CppPOPDistributionMod::POP_distrb &POP_distrbTropic;

// extern CppPOPHaloMod::pop_halo &f_POP_haloClinic;
// extern CppPOPHaloMod::pop_halo &f_POP_haloTropic;
extern int *(&POP_distrbClinic_blockLocation);
extern int *(&POP_distrbClinic_blockLocalID);
extern int *(&POP_distrbClinic_blockGlobalID);

extern int &nblocks_land;

extern int &clinicDistributionMethod;
extern int &tropicDistributionMethod;

extern char (&clinic_distribution_type)[CppPrecisionMod::CHAR_LEN];
extern char (&tropic_distribution_type)[CppPrecisionMod::CHAR_LEN];
extern char (&ew_boundary_type)[CppPrecisionMod::CHAR_LEN];
extern char (&ns_boundary_type)[CppPrecisionMod::CHAR_LEN];

extern int &nprocs_clinic;
extern int &nprocs_tropic;

extern bool &profile_barrier;
} // namespace CppDomain

#endif // LICOM_ENABLE_FORTRAN
#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_DOMAIN_H_
