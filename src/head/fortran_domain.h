#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DOMAIN_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DOMAIN_H_

#include "def-undef.h"
#include "cpp_pop_halo_mod.hpp"
#include "cpp_pop_distribution_mod.h"
#include "cpp_precision_mod.h"

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int domain_mp_nblocks_clinic_;
extern int domain_mp_nblocks_tropic_;

extern int *domain_mp_blocks_clinic_;
extern int *domain_mp_blocks_tropic_;

// extern CppPOPDistributionMod::POP_distrb domain_mp_pop_distrbclinic_;
// extern CppPOPDistributionMod::POP_distrb domain_mp_pop_distrbtropic_;

extern int domain_mp_nblocks_land_;

// blocks_land

// distrb_land

// extern CppPOPHaloMod::pop_halo domain_mp_pop_haloclinic_;
// extern CppPOPHaloMod::pop_halo domain_mp_pop_halotropic_;
extern int *domain_mp_pop_distrbclinic_blocklocation_;
extern int *domain_mp_pop_distrbclinic_blocklocalid_;
extern int *domain_mp_pop_distrbclinic_blockglobalid_;

// ltripole_grid

extern int domain_mp_clinicdistributionmethod_;
extern int domain_mp_tropicdistributionmethod_;

extern char domain_mp_clinic_distribution_type_[CppPrecisionMod::CHAR_LEN];
extern char domain_mp_tropic_distribution_type_[CppPrecisionMod::CHAR_LEN];
extern char domain_mp_ew_boundary_type_[CppPrecisionMod::CHAR_LEN];
extern char domain_mp_ns_boundary_type_[CppPrecisionMod::CHAR_LEN];

extern int domain_mp_nprocs_clinic_;
extern int domain_mp_nprocs_tropic_;

extern bool domain_mp_profile_barrier_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __domain_MOD_nblocks_clinic;
extern int __domain_MOD_nblocks_tropic;

extern int *__domain_MOD_blocks_clinic;
extern int *__domain_MOD_blocks_tropic;

// extern CppPOPDistributionMod::POP_distrb __domain_MOD_pop_distrbclinic;
// extern CppPOPDistributionMod::POP_distrb __domain_MOD_pop_distrbtropic;
extern int *__domain_MOD_pop_distrbclinic_blocklocation;
extern int *__domain_MOD_pop_distrbclinic_blocklocalid;
extern int *__domain_MOD_pop_distrbclinic_blockglobalid;

// distrb_clinic
// distrb_tropic

extern int __domain_MOD_nblocks_land;

// blocks_land

// distrb_land

// extern CppPOPHaloMod::pop_halo __domain_MOD_pop_haloclinic;
// extern CppPOPHaloMod::pop_halo __domain_MOD_pop_halotropic;

// ltripole_grid

extern int __domain_MOD_clinicdistributionmethod;
extern int __domain_MOD_tropicdistributionmethod;

extern char __domain_MOD_clinic_distribution_type[CppPrecisionMod::CHAR_LEN];
extern char __domain_MOD_tropic_distribution_type[CppPrecisionMod::CHAR_LEN];
extern char __domain_MOD_ew_boundary_type[CppPrecisionMod::CHAR_LEN];
extern char __domain_MOD_ns_boundary_type[CppPrecisionMod::CHAR_LEN];

extern int __domain_MOD_nprocs_clinic;
extern int __domain_MOD_nprocs_tropic;

extern bool __domain_MOD_profile_barrier;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DOMAIN_H_
