#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_GRID_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_GRID_H_

#include "cpp_param_mod.h"
#include "cpp_precision_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppPrecisionMod::CHAR_LEN;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double grid_mp_area_u_;
extern double grid_mp_area_t_;
extern double grid_mp_volume_u_;
extern double grid_mp_volume_t_;

extern double grid_mp_volume_t_marg_;

extern double grid_mp_area_t_marg_;

extern double grid_mp_uarea_equator_;

extern double grid_mp_area_t_k_[KM];
extern double grid_mp_volume_t_k_[KM];
extern double grid_mp_volume_t_marg_k_[KM];

extern int grid_mp_sfc_layer_type_;
extern int grid_mp_kmt_kmin_;
extern int grid_mp_n_topo_smooth_;

//extern double (*grid_mp_bath_g_)[];

// extern int (*grid_mp_kmt_g_)[IMT_GLOBAL];
// extern int (*grid_mp_basin_g_)[IMT_GLOBAL];

extern double grid_mp_dxu_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dyu_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dxt_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dyt_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dxur_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dyur_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dxtr_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_dytr_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_hts_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_htw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_hun_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_hue_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_ulat_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_ulon_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_tlat_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_tlon_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_angle_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_anglet_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_fcor_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_fcort_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_uarea_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_tarea_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_uarea_r_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_tarea_r_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_ht_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_hu_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_hur_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

//extern double (*grid_mp_dzu_)[][][];
//extern double (*grid_mp_dzt_)[][][];

extern int grid_mp_basin_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern int grid_mp_kmt_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmu_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmtold_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern double grid_mp_rcalct_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_rcalcu_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern int grid_mp_kmtn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmts_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmte_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmtw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmun_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmus_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmue_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmuw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern int grid_mp_kmtee_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int grid_mp_kmtnn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

// extern int (*grid_mp_region_mask_)[NY_BLOCK][NX_BLOCK];

extern int grid_mp_num_regions_;
extern int grid_mp_num_ms_;

extern int grid_mp_nocean_u_;
extern int grid_mp_nocean_t_;
extern int grid_mp_nsurface_u_;
extern int grid_mp_nsurface_t_;

extern char grid_mp_topography_filename_[CHAR_LEN];

extern int grid_mp_jeq_;

// extern double (*grid_mp_tlat_g_)[IMT_GLOBAL];
// extern double (*grid_mp_tlon_g_)[IMT_GLOBAL];

extern double grid_mp_at0_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_atn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_ate_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_atne_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_au0_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_aus_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_auw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double grid_mp_ausw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern char grid_mp_horiz_grid_opt_[CHAR_LEN];
extern char grid_mp_vert_grid_opt_[CHAR_LEN];
extern char grid_mp_sfc_layer_opt_[CHAR_LEN];
extern char grid_mp_topography_opt_[CHAR_LEN];
extern char grid_mp_horiz_grid_file_[CHAR_LEN];
extern char grid_mp_vert_grid_file_[CHAR_LEN];
extern char grid_mp_topography_file_[CHAR_LEN];
extern char grid_mp_region_mask_file_[CHAR_LEN];
extern char grid_mp_region_info_file_[CHAR_LEN];
extern char grid_mp_bottom_cell_file_[CHAR_LEN];
extern char grid_mp_topography_outfile_[CHAR_LEN];
extern char grid_mp_basin_grid_file_[CHAR_LEN];

extern int grid_mp_stdout_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __grid_MOD_area_u;
extern double __grid_MOD_area_t;
extern double __grid_MOD_volume_u;
extern double __grid_MOD_volume_t;

extern double __grid_MOD_volume_t_marg;

extern double __grid_MOD_area_t_marg;

extern double __grid_MOD_uarea_equator;

extern double __grid_MOD_area_t_k[KM];
extern double __grid_MOD_volume_t_k[KM];
extern double __grid_MOD_volume_t_marg_k[KM];

extern int __grid_MOD_sfc_layer_type;
extern int __grid_MOD_kmt_kmin;
extern int __grid_MOD_n_topo_smooth;

//extern double (*__grid_MOD_bath_g)[];

// extern int (*__grid_MOD_kmt_g)[IMT_GLOBAL];
// extern int (*__grid_MOD_basin_g)[IMT_GLOBAL];

extern double __grid_MOD_dxu[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dyu[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dxt[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dyt[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dxur[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dyur[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dxtr[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_dytr[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_hts[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_htw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_hun[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_hue[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_ulat[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_ulon[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_tlat[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_tlon[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_angle[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_anglet[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_fcor[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_fcort[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_uarea[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_tarea[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_uarea_r[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_tarea_r[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_ht[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_hu[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_hur[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

//extern double (*__grid_MOD_dzu)[][][];
//extern double (*__grid_MOD_dzt)[][][];

extern int __grid_MOD_basin[MAX_BLOCKS_CLINIC][JMT][IMT];

extern int __grid_MOD_kmt[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmu[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmtold[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern double __grid_MOD_rcalct[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_rcalcu[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern int __grid_MOD_kmtn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmts[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmte[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmtw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmun[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmus[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmue[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmuw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern int __grid_MOD_kmtee[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int __grid_MOD_kmtnn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

// extern int (*__grid_MOD_region_mask)[NY_BLOCK][NX_BLOCK];

extern int __grid_MOD_num_regions;
extern int __grid_MOD_num_ms;

extern int __grid_MOD_nocean_u;
extern int __grid_MOD_nocean_t;
extern int __grid_MOD_nsurface_u;
extern int __grid_MOD_nsurface_t;

extern char __grid_MOD_topography_filename[CHAR_LEN];

extern int __grid_MOD_jeq;

// extern double (*__grid_MOD_tlat_g)[IMT_GLOBAL];
// extern double (*__grid_MOD_tlon_g)[IMT_GLOBAL];

extern double __grid_MOD_at0[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_atn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_ate[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_atne[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_au0[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_aus[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_auw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __grid_MOD_ausw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern char __grid_MOD_horiz_grid_opt[CHAR_LEN];
extern char __grid_MOD_vert_grid_opt[CHAR_LEN];
extern char __grid_MOD_sfc_layer_opt[CHAR_LEN];
extern char __grid_MOD_topography_opt[CHAR_LEN];
extern char __grid_MOD_horiz_grid_file[CHAR_LEN];
extern char __grid_MOD_vert_grid_file[CHAR_LEN];
extern char __grid_MOD_topography_file[CHAR_LEN];
extern char __grid_MOD_region_mask_file[CHAR_LEN];
extern char __grid_MOD_region_info_file[CHAR_LEN];
extern char __grid_MOD_bottom_cell_file[CHAR_LEN];
extern char __grid_MOD_topography_outfile[CHAR_LEN];
extern char __grid_MOD_basin_grid_file[CHAR_LEN];

extern int __grid_MOD_stdout;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_GRID_H_
