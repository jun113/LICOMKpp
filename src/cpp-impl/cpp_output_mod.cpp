#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_output_mod.h"

namespace CppOutputMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
#ifdef DAILYACC
float (&z0daily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_z0daily_;

float (&mlddaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_mlddaily_;
float (&lthfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_lthfdaily_;
float (&sshfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_sshfdaily_;
float (&lwvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_lwvdaily_;
float (&swvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_swvdaily_;
float (&precdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_precdaily_;
float (&evapdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_evapdaily_;
float (&roffdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_roffdaily_;
float (&sudaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_sudaily_;
float (&svdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_svdaily_;
float (&tsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_tsdaily_;
float (&ssdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_ssdaily_;
float (&usdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_usdaily_;
float (&vsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_vsdaily_;
float (&wsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_wsdaily_;

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
float (&z0mon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_z0mon_;
float (&himon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_himon_;
float (&hdmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_hdmon_;
float (&qicemon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_qicemon_;
float (&lthfmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_lthfmon_;
float (&sshfmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_sshfmon_;
float (&lwvmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_lwvmon_;
float (&swvmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_swvmon_;
float (&sumon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_sumon_;
float (&svmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_svmon_;
float (&runoffmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_runoffmon_;
float (&freshmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_freshmon_;
float (&wsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_wsmon_;
float (&tsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_tsmon_;
float (&ssmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_ssmon_;
float (&usmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_usmon_;
float (&vsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = output_mod_mp_vsmon_;
float (&icmon)[MAX_BLOCKS_CLINIC][2][JMT][IMT] = output_mod_mp_icmon_;

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
float (&athkdfmon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_athkdfmon_;
#endif // ISO_TYPE_BF

#ifdef ISO
float (&axmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_axmon_iso_;
float (&aymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_aymon_iso_;
float (&azmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_azmon_iso_;
float (&dxmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_dxmon_iso_;
float (&dymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_dymon_iso_;
float (&dzmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = output_mod_mp_dzmon_iso_;
float (&wsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_wsmon_iso_;
float (&vsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_vsmon_iso_;
#ifdef ISOOUT
float (&azmon_vetisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_azmon_vetisomon_;
float (&azmon_vntisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_azmon_vntisomon_;
float (&azmon_vbtisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = output_mod_mp_azmon_vbtisomon_;
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
float (&azmon_a3mon)[KM][JMT][IMT] = output_mod_mp_azmon_a3mon_;
#endif // SMAG_OUT

double &err_norm2mon = output_mod_mp_err_norm2mon_;
#endif // LOWRES
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
#ifdef DAILYACC
float (&z0daily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_z0daily;

float (&mlddaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_mlddaily;
float (&lthfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_lthfdaily;
float (&sshfdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_sshfdaily;
float (&lwvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_lwvdaily;
float (&swvdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_swvdaily;
float (&precdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_precdaily;
float (&evapdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_evapdaily;
float (&roffdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_roffdaily;
float (&sudaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_sudaily;
float (&svdaily)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_svdaily;
float (&tsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_tsdaily;
float (&ssdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_ssdaily;
float (&usdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_usdaily;
float (&vsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_vsdaily;
float (&wsdaily)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_wsdaily;

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
float (&z0mon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_z0mon;
float (&himon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_himon;
float (&hdmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_hdmon;
float (&qicemon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_qicemon;
float (&lthfmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_lthfmon;
float (&sshfmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_sshfmon;
float (&lwvmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_lwvmon;
float (&swvmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_swvmon;
float (&sumon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_sumon;
float (&svmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_svmon;
float (&runoffmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_runoffmon;
float (&freshmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_freshmon;
float (&wsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_wsmon;
float (&tsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_tsmon;
float (&ssmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_ssmon;
float (&usmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_usmon;
float (&vsmon)[MAX_BLOCKS_CLINIC][JMT][IMT] = __output_mod_MOD_vsmon;
float (&icmon)[MAX_BLOCKS_CLINIC][2][JMT][IMT] = __output_mod_MOD_icmon;

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
float (&athkdfmon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_athkdfmon;
#endif // ISO_TYPE_BF

#ifdef ISO
float (&axmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_axmon_iso;
float (&aymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_aymon_iso;
float (&azmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_azmon_iso;
float (&dxmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_dxmon_iso;
float (&dymon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_dymon_iso;
float (&dzmon_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __output_mod_MOD_dzmon_iso;
float (&wsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_wsmon_iso;
float (&vsmon_iso)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_vsmon_iso;
#ifdef ISOOUT
float (&azmon_vetisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_azmon_vetisomon;
float (&azmon_vntisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_azmon_vntisomon;
float (&azmon_vbtisomon)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __output_mod_MOD_azmon_vbtisomon;
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
float (&azmon_a3mon)[KM][JMT][IMT] = __output_mod_MOD_azmon_a3mon;
#endif // SMAG_OUT

double &err_norm2mon = __output_mod_MOD_err_norm2mon;
#endif // LOWRES
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppOutputMod
#endif // LICOM_ENABLE_FORTRAN
