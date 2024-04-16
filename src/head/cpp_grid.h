#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_GRID_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_GRID_H_

#include "cpp_param_mod.h"
#include "cpp_precision_mod.h"

namespace CppGrid {

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppPrecisionMod::CHAR_LEN;

extern double* dxyur;
extern double* h_tu_swen;
extern int*    kmt_nsew;

extern double &area_u;
extern double &area_t;
                                                         
extern double &volume_u;
extern double &volume_t;
                                                         
extern double &volume_t_marg;
                                                         
extern double &area_t_marg;
                                                         
extern double &uarea_equator;
                                                         
extern double (&area_t_k)[KM];
                                                         
extern double (&volume_t_k)[KM];
extern double (&volume_t_marg_k)[KM];
                                                         
                                                         
extern int &sfc_layer_type;
extern int &kmt_kmin;
extern int &n_topo_smooth;

constexpr int SFC_LAYER_VARTHICK = 1;
constexpr int SFC_LAYER_RIGID    = 2;
constexpr int SFC_LAYER_OLDFREE  = 3;

// extern int (*&kmt_g)[IMT_GLOBAL];
// extern int (*&basin_g)[IMT_GLOBAL];

extern double (&dxu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dyu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dxt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dyt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dxur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dyur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dxtr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&dytr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&hts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&htw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&hun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&hue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ulat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ulon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&tlat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&tlon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&angle)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&anglet)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&fcor)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&fcort)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&uarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&tarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&uarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&tarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ht)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&hu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&hur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

//extern double (*grid_mp_dzu_)[][][];
//extern double (*grid_mp_dzt_)[][][];

extern int (&basin)[MAX_BLOCKS_CLINIC][JMT][IMT];
                                                         
extern int (&kmt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmtold)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
                                                         
extern double (&rcalct)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&rcalcu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
                                                         
extern int (&kmtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmuw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
                                                         
extern int (&kmtee)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern int (&kmtnn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
                                                         
// extern int (*&region_mask)[NY_BLOCK][NX_BLOCK];


constexpr int MAX_REGIONS = 15;
constexpr int MAX_MS      = 7;

extern int &num_regions;
extern int &num_ms;
                                                         
extern int &nocean_u;
extern int &nocean_t;
extern int &nsurface_u;
extern int &nsurface_t;
                                                         
extern char (&topography_filename)[CHAR_LEN];
                                                         
extern int &jeq;
                                                         
// extern double (*&tlat_g)[IMT_GLOBAL];
// extern double (*&tlon_g)[IMT_GLOBAL];
                                                         
extern double (&at0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&atn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ate)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&atne)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&au0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&aus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&auw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double (&ausw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
                                                         
extern char (&horiz_grid_opt)[CHAR_LEN];
extern char (&vert_grid_opt)[CHAR_LEN];
                                                         
extern char sfc_layer_opt[CHAR_LEN];
extern char topography_opt[CHAR_LEN];
extern char horiz_grid_file[CHAR_LEN];
extern char vert_grid_file[CHAR_LEN];
extern char topography_file[CHAR_LEN];
extern char region_mask_file[CHAR_LEN];
extern char region_info_file[CHAR_LEN];
extern char bottom_cell_file[CHAR_LEN];
extern char topography_outfile[CHAR_LEN];
extern char basin_grid_file[CHAR_LEN];
                                                         
extern int &stdout;

} // namespace CppGrid

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_GRID_H_
