#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_isopyc_mod.h"

namespace CppIsopycMod {
#ifdef ISO

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&dptlim)[NRPL+1]                                = isopyc_mod_mp_dptlim_;

double (&fzisop)[KM]                                    = isopyc_mod_mp_fzisop_;

double &slmxr                                           = isopyc_mod_mp_slmxr_;

double (&ahisop)[MAX_BLOCKS_CLINIC][JMT][IMT]           = isopyc_mod_mp_ahisop_;
double (&athkdf)[MAX_BLOCKS_CLINIC][JMT][IMT]           = isopyc_mod_mp_athkdf_;

// double (*(&e))[3][JMT][KMP1][IMT]                       = isopyc_mod_mp_e_;
// double (*(&rhoi))[NRPL][JMT][KM+1][IMT]                 = isopyc_mod_mp_rhoi_;

// double (*(&k1))[1][JMT][KM+1][IMT]                      = isopyc_mod_mp_k1_;
// double (*(&k2))[1][JMT][KM+1][IMT]                      = isopyc_mod_mp_k2_;
// double (*(&k3))[3][JMT][KM+1][IMT]                      = isopyc_mod_mp_k3_;

// double (*(&adv_vetiso))[JMT][KM][IMT]                       = isopyc_mod_mp_adv_vetiso_;

double (&adv_vbtiso)[MAX_BLOCKS_CLINIC][JMT][KM+1][IMT] = isopyc_mod_mp_adv_vbtiso_;
double (&adv_vntiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]   = isopyc_mod_mp_adv_vntiso_;

#ifdef isopycmixspatialvar
float (&dciso1)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = isopyc_mod_mp_dciso1_;
float (&dciso2)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = isopyc_mod_mp_dciso2_;

float &dslope                                           = isopyc_mod_mp_dslope_;
float &slopec                                           = isopyc_mod_mp_slopec_;
#endif // isopycmixspatialvar

double (&kisrpl)[KM]                                    = isopyc_mod_mp_kisrpl_;

int (&krplin)[NRPL]                                     = isopyc_mod_mp_krplin_;

double (&zt)[KM]                                        = isopyc_mod_mp_zt_;

double (&dzw)[KM+1]                                     = isopyc_mod_mp_dzw_;
double (&dzwr)[KM+1]                                    = isopyc_mod_mp_dzwr_;

double (&dzr)[KM]                                       = isopyc_mod_mp_dzr_;

double (&tmask)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = isopyc_mod_mp_tmask_;
double (&f3)[MAX_BLOCKS_CLINIC][JMT][IMT]               = isopyc_mod_mp_f3_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&dptlim)[NRPL+1]                                = __isopyc_mod_MOD_dptlim;

double (&fzisop)[KM]                                    = __isopyc_mod_MOD_fzisop;

double &slmxr                                           = __isopyc_mod_MOD_slmxr;

double (&ahisop)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __isopyc_mod_MOD_ahisop;
double (&athkdf)[MAX_BLOCKS_CLINIC][JMT][IMT]           = __isopyc_mod_MOD_athkdf;

// double (*(&e))[3][JMT][KMP1][IMT]                       = __isopyc_mod_MOD_e;
// double (*(&rhoi))[NRPL][JMT][KM+1][IMT]                 = __isopyc_mod_MOD_rhoi;

// double (*(&k1))[1][JMT][KM+1][IMT]                      = __isopyc_mod_MOD_k1;
// double (*(&k2))[1][JMT][KM+1][IMT]                      = __isopyc_mod_MOD_k2;
// double (*(&k3))[3][JMT][KM+1][IMT]                      = __isopyc_mod_MOD_k3;

// double (*(&adv_vetiso))[JMT][KM][IMT]                       = __isopyc_mod_MOD_adv_vetiso;

double (&adv_vbtiso)[MAX_BLOCKS_CLINIC][JMT][KM+1][IMT] = __isopyc_mod_MOD_adv_vbtiso;
double (&adv_vntiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]   = __isopyc_mod_MOD_adv_vntiso;

#ifdef isopycmixspatialvar
float (&dciso1)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = __isopyc_mod_MOD_dciso1;
float (&dciso2)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = __isopyc_mod_MOD_dciso2;

float &dslope                                           = __isopyc_mod_MOD_dslope;
float &slopec                                           = __isopyc_mod_MOD_slopec;
#endif // isopycmixspatialvar

double (&kisrpl)[KM]                                    = __isopyc_mod_MOD_kisrpl;

int (&krplin)[NRPL]                                     = __isopyc_mod_MOD_krplin;

double (&zt)[KM]                                        = __isopyc_mod_MOD_zt;

double (&dzw)[KM+1]                                     = __isopyc_mod_MOD_dzw;
double (&dzwr)[KM+1]                                    = __isopyc_mod_MOD_dzwr;

double (&dzr)[KM]                                       = __isopyc_mod_MOD_dzr;

double (&tmask)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]        = __isopyc_mod_MOD_tmask;
double (&f3)[MAX_BLOCKS_CLINIC][JMT][IMT]               = __isopyc_mod_MOD_f3;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // ISO

} // namespace CppIsopycMod
#endif // LICOM_ENABLE_FORTRAN
