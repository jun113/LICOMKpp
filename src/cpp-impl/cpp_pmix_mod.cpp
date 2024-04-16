#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_pmix_mod.h"

namespace CppPmixMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMM1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &rtst                                                 = pmix_mod_mp_rtst_;
int &rtend                                                = pmix_mod_mp_rtend_;
int &rust                                                 = pmix_mod_mp_rust_;
int &ruend                                                = pmix_mod_mp_ruend_;

double &wndmix                                            = pmix_mod_mp_wndmix_;
double &fricmx                                            = pmix_mod_mp_fricmx_;

double &diff_cbt_back                                     = pmix_mod_mp_diff_cbt_back_;
double &diff_cbt_limit                                    = pmix_mod_mp_diff_cbt_limit_;

double &visc_cbu_back                                     = pmix_mod_mp_visc_cbu_back_;
double &visc_cbu_limit                                    = pmix_mod_mp_visc_cbu_limit_;

// double (*(&ric))[KMM1][JMT][IMT]          = pmix_mod_mp_ric_;
// double (*(&rict))[KMM1][JMT][IMT]         = pmix_mod_mp_rict_;
// double (*(&rict_replace))[KMM1][JMT][IMT] = pmix_mod_mp_rict_replace_;


// double (*(&rict_ref))[JMT][IMT]                           = pmix_mod_mp_rict_ref_;

// double (*(&rit))[KMM1][JMT][IMT]                          = pmix_mod_mp_rit_;
// double (*(&riu))[KM+1][JMT][IMT]                          = pmix_mod_mp_riu_;

double (&ricdt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]        = pmix_mod_mp_ricdt_;
double (&ricdttms)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]     = pmix_mod_mp_ricdttms_;

double (&ridt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]         = pmix_mod_mp_ridt_;

double (&s2u)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]          = pmix_mod_mp_s2u_;
double (&s2t)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]          = pmix_mod_mp_s2t_;

#ifdef SOLAR
double (&pen)[KMM1]                                       = pmix_mod_mp_pen_;
#endif // SOLAR

#ifdef SOLARCHLORO
double (&pen_chl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]        = pmix_mod_mp_pen_chl_;
#endif // SOLARCHLORO
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &rtst                                                 = __pmix_mod_MOD_rtst;
int &rtend                                                = __pmix_mod_MOD_rtend;
int &rust                                                 = __pmix_mod_MOD_rust;
int &ruend                                                = __pmix_mod_MOD_ruend;

double &wndmix                                            = __pmix_mod_MOD_wndmix;
double &fricmx                                            = __pmix_mod_MOD_fricmx;

double &diff_cbt_back                                     = __pmix_mod_MOD_diff_cbt_back;
double &diff_cbt_limit                                    = __pmix_mod_MOD_diff_cbt_limit;

double &visc_cbu_back                                     = __pmix_mod_MOD_visc_cbu_back;
double &visc_cbu_limit                                    = __pmix_mod_MOD_visc_cbu_limit;

// double (*(&ric))[KMM1][JMT][IMT]                          = __pmix_mod_MOD_ric;
// double (*(&rict))[KMM1][JMT][IMT]                         = __pmix_mod_MOD_rict;
// double (*(&rict_replace))[KMM1][JMT][IMT]                 = __pmix_mod_MOD_rict_replace;

// double (*(&rict_ref))[JMT][IMT]                           = __pmix_mod_MOD_rict_ref;

// double (*(&rit))[KMM1][JMT][IMT]                          = __pmix_mod_MOD_rit;
// double (*(&riu))[KM+1][JMT][IMT]                          = __pmix_mod_MOD_riu;

double (&ricdt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]        = __pmix_mod_MOD_ricdt;
double (&ricdttms)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]     = __pmix_mod_MOD_ricdttms;

double (&ridt)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]         = __pmix_mod_MOD_ridt;

double (&s2u)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]          = __pmix_mod_MOD_s2u;
double (&s2t)[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT]          = __pmix_mod_MOD_s2t;

#ifdef SOLAR
double (&pen)[KMM1]                                       = __pmix_mod_MOD_pen;
#endif // SOLAR

#ifdef SOLARCHLORO
double (&pen_chl)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]        = __pmix_mod_MOD_pen_chl;
#endif // SOLARCHLORO
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
double ric[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
double rict[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
double rict_replace[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

double rit[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
} // namespace CppPmixMod
#endif // LICOM_ENABLE_FORTRAN
