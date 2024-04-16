#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_hmix_del4.h"
namespace CppHmixDel4 {
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;

double* dt_nsew   = nullptr;
double* du_cnsewm = nullptr;
double* dm_cnsew  = nullptr;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&dtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dtn_;
double (&dts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dts_;
double (&dte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dte_;
double (&dtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dtw_;
double (&duc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_duc_;
double (&dun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dun_;
double (&dus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dus_;
double (&due)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_due_;
double (&duw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_duw_;
double (&dmc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dmc_;
double (&dmn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dmn_;
double (&dms)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dms_;
double (&dme)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dme_;
double (&dmw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dmw_;
double (&dum)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_dum_;
double (&ahf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_ahf_;
double (&amf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_amf_;
double (&ratio_dxy)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del4_mp_del4_ratio_dxy_;

double &ah = hmix_del4_mp_ah_;
double &am = hmix_del4_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&dtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dtn;
double (&dts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dts;
double (&dte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dte;
double (&dtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dtw;
double (&duc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_duc;
double (&dun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dun;
double (&dus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dus;
double (&due)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_due;
double (&duw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_duw;
double (&dmc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dmc;
double (&dmn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dmn;
double (&dms)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dms;
double (&dme)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dme;
double (&dmw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dmw;
double (&dum)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_dum;
double (&ahf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_ahf;
double (&amf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_amf;
double (&ratio_dxy)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del4_MOD_del4_ratio_dxy;

double &ah = __hmix_del4_MOD_ah;
double &am = __hmix_del4_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // CppHmixDel4
#endif // LICOM_ENABLE_FORTRAN
