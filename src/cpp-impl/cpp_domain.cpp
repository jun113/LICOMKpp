#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/fortran_domain.h"
#include "../head/cpp_precision_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_distribution_mod.h"

namespace CppDomain {

using namespace CppPrecisionMod;

CppPOPDistributionMod::POP_distrb POP_distrbClinic;

CppPOPHaloMod::pop_halo POP_haloClinic;
CppPOPHaloMod::POPHalo POP_haloClinic_C;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &nblocks_clinic                        = domain_mp_nblocks_clinic_;
int &nblocks_tropic                        = domain_mp_nblocks_tropic_;

int (*(&blocks_clinic)) = domain_mp_blocks_clinic_;
int (*(&blocks_tropic)) = domain_mp_blocks_tropic_;

// CppPOPDistributionMod::POP_distrb &POP_distrbClinic = domain_mp_pop_distrbclinic_;
// CppPOPDistributionMod::POP_distrb &POP_distrbTropic = domain_mp_pop_distrbtropic_;

// CppPOPHaloMod::pop_halo &f_POP_haloClinic = domain_mp_pop_haloclinic_;
// CppPOPHaloMod::pop_halo &f_POP_haloTropic = domain_mp_pop_halotropic_;
int (*(&POP_distrbClinic_blockLocation)) = domain_mp_pop_distrbclinic_blocklocation_;
int (*(&POP_distrbClinic_blockLocalID))  = domain_mp_pop_distrbclinic_blocklocalid_;
int (*(&POP_distrbClinic_blockGlobalID)) = domain_mp_pop_distrbclinic_blockglobalid_;

// distrb_clinic
// distrb_tropic

int &nblocks_land = domain_mp_nblocks_land_;

// blocks_land

// distrb_land

// pop_haloclinic
// pop_halotropic

// ltripole_grid

int &clinicDistributionMethod              = domain_mp_clinicdistributionmethod_;
int &tropicDistributionMethod              = domain_mp_tropicdistributionmethod_;

char (&clinic_distribution_type)[CHAR_LEN] = domain_mp_clinic_distribution_type_;
char (&tropic_distribution_type)[CHAR_LEN] = domain_mp_tropic_distribution_type_;
char (&ew_boundary_type)[CHAR_LEN]         = domain_mp_ew_boundary_type_;
char (&ns_boundary_type)[CHAR_LEN]         = domain_mp_ns_boundary_type_;

int &nprocs_clinic                         = domain_mp_nprocs_clinic_;
int &nprocs_tropic                         = domain_mp_nprocs_tropic_;

bool &profile_barrier                      = domain_mp_profile_barrier_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &nblocks_clinic                        = __domain_MOD_nblocks_clinic;
int &nblocks_tropic                        = __domain_MOD_nblocks_tropic;

int (*(&blocks_clinic)) 									 = __domain_MOD_blocks_clinic;
int (*(&blocks_tropic))                    = __domain_MOD_blocks_tropic;

// CppPOPDistributionMod::POP_distrb &POP_distrbClinic = __domain_MOD_pop_distrbclinic;
// CppPOPDistributionMod::POP_distrb &POP_distrbTropic = __domain_MOD_pop_distrbtropic;

// CppPOPHaloMod::pop_halo &f_POP_haloClinic = __domain_MOD_pop_haloclinic;
// CppPOPHaloMod::pop_halo &f_POP_haloTropic = __domain_MOD_pop_halotropic;
int (*(&POP_distrbClinic_blockLocation)) = __domain_MOD_pop_distrbclinic_blocklocation;
int (*(&POP_distrbClinic_blockLocalID))  = __domain_MOD_pop_distrbclinic_blocklocalid;
int (*(&POP_distrbClinic_blockGlobalID)) = __domain_MOD_pop_distrbclinic_blockglobalid;

// distrb_clinic
// distrb_tropic

int &nblocks_land                          = __domain_MOD_nblocks_land;

// blocks_land

// distrb_land

// pop_haloclinic
// pop_halotropic

// ltripole_grid

int &clinicDistributionMethod              = __domain_MOD_clinicdistributionmethod;
int &tropicDistributionMethod              = __domain_MOD_tropicdistributionmethod;

char (&clinic_distribution_type)[CHAR_LEN] = __domain_MOD_clinic_distribution_type;
char (&tropic_distribution_type)[CHAR_LEN] = __domain_MOD_tropic_distribution_type;
char (&ew_boundary_type)[CHAR_LEN]         = __domain_MOD_ew_boundary_type;
char (&ns_boundary_type)[CHAR_LEN]         = __domain_MOD_ns_boundary_type;

int &nprocs_clinic                         = __domain_MOD_nprocs_clinic;
int &nprocs_tropic                         = __domain_MOD_nprocs_tropic;

bool &profile_barrier                      = __domain_MOD_profile_barrier;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppDomain

#endif // LICOM_ENABLE_FORTRAN
