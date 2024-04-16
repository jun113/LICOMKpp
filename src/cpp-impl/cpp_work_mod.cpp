#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_work_mod.h"

namespace CppWorkMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::KMP1;
using CppParamMod::IMT_GLOBAL;

double* work_1  = nullptr;
double* work_2  = nullptr;
double* work_3  = nullptr;
double* temp    = nullptr;
double* uk      = nullptr;
double* vk      = nullptr;
double* wkb     = nullptr;
double* wkc     = nullptr;
double* wkd     = nullptr;
double* tf      = nullptr;
double* stf     = nullptr;
double* work1_g = nullptr;
double* work2_g = nullptr;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&pxb)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_pxb_;
double (&pyb)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_pyb_;
double (&pax)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_pax_;
double (&pay)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_pay_;
double (&whx)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_whx_;
double (&why)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_why_;
double (&wgp)[MAX_BLOCKS_CLINIC][JMT][IMT]     = work_mod_mp_wgp_;
double (&wka)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = work_mod_mp_wka_;

// double (*&work_1)[KM][JMT][IMT]                = work_mod_mp_work_1_;
// double (*&work_2)[KM][JMT][IMT]                = work_mod_mp_work_2_;
// double (*&work_3)[KM+1][JMT][IMT]              = work_mod_mp_work_3_;
                                               
// double (*&temp)[KM][JMT][IMT]                  = work_mod_mp_temp_;
                                               
// double (*&uk)[KM][JMT][IMT]                    = work_mod_mp_uk_;
// double (*&vk)[KM][JMT][IMT]                    = work_mod_mp_vk_;
                                               
double (&work)[MAX_BLOCKS_CLINIC][JMT][IMT]    = work_mod_mp_work_;
                                              
double (&wkk)[KMP1]                            = work_mod_mp_wkk_;
                                               
// double (*&wkb)[KM][JMT][IMT]                   = work_mod_mp_wkb_;
// double (*&wkc)[KM][JMT][IMT]                   = work_mod_mp_wkc_;
// double (*&wkd)[KM][JMT][IMT]                   = work_mod_mp_wkd_;
                                               
// double (*&tf)[KM][JMT][IMT]                    = work_mod_mp_tf_;
// double (*&stf)[JMT][IMT]                       = work_mod_mp_stf_;
                                               
// float (*&buffer_real4)[IMT_GLOBAL]             = work_mod_mp_buffer_real4_;
                                               
// double (*&work1_g)[IMT_GLOBAL]                 = work_mod_mp_work1_g_;
// double (*&work2_g)[IMT_GLOBAL]                 = work_mod_mp_work2_g_;
// double (*&work3_g)[IMT_GLOBAL]                 = work_mod_mp_work3_g_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&pxb)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_pxb;
double (&pyb)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_pyb;
double (&pax)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_pax;
double (&pay)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_pay;
double (&whx)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_whx;
double (&why)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_why;
double (&wgp)[MAX_BLOCKS_CLINIC][JMT][IMT]      = __work_mod_MOD_wgp;
double (&wka)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __work_mod_MOD_wka;
                                               
// double (*&work_1)[KM][JMT][IMT]                 = __work_mod_MOD_work_1;
// double (*&work_2)[KM][JMT][IMT]                 = __work_mod_MOD_work_2;
// double (*&work_3)[KM+1][JMT][IMT]               = __work_mod_MOD_work_3;
                                                 
// double (*&temp)[KM][JMT][IMT]                   = __work_mod_MOD_temp;
                                                 
// double (*&uk)[KM][JMT][IMT]                     = __work_mod_MOD_uk;
// double (*&vk)[KM][JMT][IMT]                     = __work_mod_MOD_vk;
                                               
double (&work)[MAX_BLOCKS_CLINIC][JMT][IMT]     = __work_mod_MOD_work;
                                               
double (&wkk)[KMP1]                             = __work_mod_MOD_wkk;
                                               
// double (*&wkb)[KM][JMT][IMT]                    = __work_mod_MOD_wkb;
// double (*&wkc)[KM][JMT][IMT]                    = __work_mod_MOD_wkc;
// double (*&wkd)[KM][JMT][IMT]                    = __work_mod_MOD_wkd;
                                                 
// double (*&tf)[KM][JMT][IMT]                     = __work_mod_MOD_tf;
// double (*&stf)[JMT][IMT]                        = __work_mod_MOD_stf;
                                                 
// float (*&buffer_real4)[IMT_GLOBAL]              = __work_mod_MOD_buffer_real4;
                                                 
// double (*&work1_g)[IMT_GLOBAL]                  = __work_mod_MOD_work1_g;
// double (*&work2_g)[IMT_GLOBAL]                  = __work_mod_MOD_work2_g;
// double (*&work3_g)[IMT_GLOBAL]                  = __work_mod_MOD_work3_g;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppWorkMod
#endif // LICOM_ENABLOE_FORTRAN
