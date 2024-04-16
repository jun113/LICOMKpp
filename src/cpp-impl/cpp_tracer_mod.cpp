#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_tracer_mod.h"

namespace CppTracerMod {

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&atb)[MAX_BLOCKS_CLINIC][NTRA][KM+1][JMT][IMT]      = tracer_mod_mp_atb_;
double (&net)[MAX_BLOCKS_CLINIC][NTRA][JMT][IMT]            = tracer_mod_mp_net_;

double (&at)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_at_;
// double (*&restore_at)[NTRA][JMT][IMT]                       = tracer_mod_mp_restore_at_;

double (&pdensity)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = tracer_mod_mp_pdensity_;
double (&amld)[MAX_BLOCKS_CLINIC][JMT][IMT]                 = tracer_mod_mp_amld_;

double (&tend)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]       = tracer_mod_mp_tend_;
double (&ax)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_ax_;
double (&ay)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_ay_;
double (&az)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_az_;

double (&dx)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_dx_;
double (&dy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_dy_;
double (&dz)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = tracer_mod_mp_dz_;

double (&penetrate)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]        = tracer_mod_mp_penetrate_;

double (&dt_diff)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = tracer_mod_mp_dt_diff_;

double (&ddy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]        = tracer_mod_mp_ddy_;

double (&dt_conv)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = tracer_mod_mp_dt_conv_;
double (&dt_restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = tracer_mod_mp_dt_restore_;

#ifdef ISO
double (&aay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = tracer_mod_mp_aay_iso_;
double (&ddy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = tracer_mod_mp_ddy_iso_;

double (&ax_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_ax_iso_;
double (&ay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_ay_iso_;
double (&az_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_az_iso_;

double (&dx_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_dx_iso_;
double (&dy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_dy_iso_;
double (&dz_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = tracer_mod_mp_dz_iso_;
#endif // ISO

double (&licomqice)[MAX_BLOCKS_CLINIC][JMT][IMT]            = tracer_mod_mp_licomqice_;

double &fw_norm2                                            = tracer_mod_mp_fw_norm2_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&atb)[MAX_BLOCKS_CLINIC][NTRA][KM+1][JMT][IMT]      = __tracer_mod_MOD_atb;
double (&net)[MAX_BLOCKS_CLINIC][NTRA][JMT][IMT]            = __tracer_mod_MOD_net;

double (&at)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_at;
// double (*&restore_at)[NTRA][JMT][IMT]                       = __tracer_mod_MOD_restore_at;

double (&pdensity)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = __tracer_mod_MOD_pdensity;
double (&amld)[MAX_BLOCKS_CLINIC][JMT][IMT]                 = __tracer_mod_MOD_amld;

double (&tend)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]       = __tracer_mod_MOD_tend;
double (&ax)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_ax;
double (&ay)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_ay;
double (&az)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_az;

double (&dx)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_dx;
double (&dy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_dy;
double (&dz)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]         = __tracer_mod_MOD_dz;

double (&penetrate)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]        = __tracer_mod_MOD_penetrate;

double (&dt_diff)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = __tracer_mod_MOD_dt_diff;

double (&ddy)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]        = __tracer_mod_MOD_ddy;

double (&dt_conv)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = __tracer_mod_MOD_dt_conv;
double (&dt_restore)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT] = __tracer_mod_MOD_dt_restore;

#ifdef ISO
double (&aay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = __tracer_mod_MOD_aay_iso;
double (&ddy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = __tracer_mod_MOD_ddy_iso;

double (&ax_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_ax_iso;
double (&ay_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_ay_iso;
double (&az_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_az_iso;

double (&dx_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_dx_iso;
double (&dy_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_dy_iso;
double (&dz_iso)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]     = __tracer_mod_MOD_dz_iso;
#endif // ISO

double (&licomqice)[MAX_BLOCKS_CLINIC][JMT][IMT]            = __tracer_mod_MOD_licomqice;

double &fw_norm2                                            = __tracer_mod_MOD_fw_norm2;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppTracerMod
#endif // LICOM_ENABLE_FORTRAN
