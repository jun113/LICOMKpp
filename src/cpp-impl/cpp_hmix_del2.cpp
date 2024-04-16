#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_hmix_del2.h"
namespace CppHmixDel2 {
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (*(&dtn))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dtn_;
double (*(&dts))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dts_;
double (*(&dte))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dte_;
double (*(&dtw))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dtw_;
double (*(&duc))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_duc_;
double (*(&dun))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dun_;
double (*(&dus))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dus_;
double (*(&due))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_due_;
double (*(&duw))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_duw_;
double (*(&dmc))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmc_;
double (*(&dmn))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmn_;
double (*(&dms))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dms_;
double (*(&dme))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dme_;
double (*(&dmw))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmw_;
double (*(&dum))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dum_;
double (*(&ahf))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_ahf_;
double (*(&amf))[NY_BLOCK][NX_BLOCK] = hmix_del2_mp_amf_;

double &ah = hmix_del2_mp_ah_;
double &am = hmix_del2_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (*(&dtn))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dtn;
double (*(&dts))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dts;
double (*(&dte))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dte;
double (*(&dtw))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dtw;
double (*(&duc))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_duc;
double (*(&dun))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dun;
double (*(&dus))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dus;
double (*(&due))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_due;
double (*(&duw))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_duw;
double (*(&dmc))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmc;
double (*(&dmn))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmn;
double (*(&dms))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dms;
double (*(&dme))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dme;
double (*(&dmw))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmw;
double (*(&dum))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dum;
double (*(&ahf))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_ahf;
double (*(&amf))[NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_amf;

double &ah = __hmix_del2_MOD_ah;
double &am = __hmix_del2_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // CppHmixDel2
#endif // LICOM_ENABLE_FORTRAN
