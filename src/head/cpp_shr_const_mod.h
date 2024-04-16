#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_SHR_CONST_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_SHR_CONST_MOD_H_

#include "def-undef.h"

namespace CppShrConstMod { 

constexpr double SHR_CONST_PI      = 3.14159265358979323846;    
constexpr double SHR_CONST_CDAY    = 86400.0;
constexpr double SHR_CONST_SDAY    = 86164.0;
constexpr double SHR_CONST_OMEGA   = 2.0 * SHR_CONST_PI / SHR_CONST_SDAY;
constexpr double SHR_CONST_REARTH  = 6.37122e6;
constexpr double SHR_CONST_G       = 9.80616;
constexpr double SHR_CONST_STEBOL  = 5.67e-8;
constexpr double SHR_CONST_BOLTZ   = 1.38065e-23; 
constexpr double SHR_CONST_AVOGAD  = 6.02214e26;
constexpr double SHR_CONST_RGAS    = SHR_CONST_AVOGAD * SHR_CONST_BOLTZ ;
constexpr double SHR_CONST_MWDAIR  = 28.966;
constexpr double SHR_CONST_MWWV    = 18.016;
constexpr double SHR_CONST_RDAIR   = SHR_CONST_RGAS / SHR_CONST_MWDAIR;
constexpr double SHR_CONST_RWV     = SHR_CONST_RGAS / SHR_CONST_MWWV;
constexpr double SHR_CONST_ZVIR    = (SHR_CONST_RWV / SHR_CONST_RDAIR) - 1.0;
constexpr double SHR_CONST_KARMAN  = 0.4; 
constexpr double SHR_CONST_PSTD    = 101325.0;
constexpr double SHR_CONST_PDB     = 0.0112372;
constexpr double SHR_CONST_TKTRIP  = 273.16;
constexpr double SHR_CONST_TKFRZ   = 273.15;
constexpr double SHR_CONST_TKFRZSW = SHR_CONST_TKFRZ - 1.8;
constexpr double SHR_CONST_RHODAIR = SHR_CONST_PSTD / (SHR_CONST_RDAIR * SHR_CONST_TKFRZ);
constexpr double SHR_CONST_RHOFW   = 1.000e3;
constexpr double SHR_CONST_RHOSW   = 1.026e3;
constexpr double SHR_CONST_RHOICE  = 0.917e3;
constexpr double SHR_CONST_CPDAIR  = 1.00464e3;
constexpr double SHR_CONST_CPWV    = 1.810e3;
constexpr double SHR_CONST_CPVIR   = (SHR_CONST_CPWV / SHR_CONST_CPDAIR) - 1.0;
constexpr double SHR_CONST_CPFW    = 4.188e3;
constexpr double SHR_CONST_CPSW    = 3.996e3;
constexpr double SHR_CONST_CPICE   = 2.11727e3;
constexpr double SHR_CONST_LATICE  = 3.337e5;
constexpr double SHR_CONST_LATVAP  = 2.501e6;
constexpr double SHR_CONST_LATSUB  = SHR_CONST_LATICE + SHR_CONST_LATVAP;
constexpr double SHR_CONST_OCN_REF_SAL = 34.7;
constexpr double SHR_CONST_ICE_REF_SAL =  4.0;
constexpr double SHR_CONST_SPVAL   = 1.0e30;
} // namespace CppConstantMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_CONSTANT_MOD_H_
