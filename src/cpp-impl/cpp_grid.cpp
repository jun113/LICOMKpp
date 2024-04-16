#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/fortran_grid.h"

namespace CppGrid {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppPrecisionMod::CHAR_LEN;

double* dxyur     = nullptr;
double* h_tu_swen = nullptr;
int*    kmt_nsew  = nullptr;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double &area_u                                           = grid_mp_area_u_;
double &area_t                                           = grid_mp_area_t_;
                                                         
double &volume_u                                         = grid_mp_volume_u_;
double &volume_t                                         = grid_mp_volume_t_;
                                                         
// double &volume_t_marg                                    = grid_mp_volume_t_marg_;
                                                         
// double &area_t_marg                                      = grid_mp_area_t_marg_;
                                                         
// double &uarea_equator                                    = grid_mp_uarea_equator_;
                                                         
// double (&area_t_k)[KM]                                   = grid_mp_area_t_k_;
                                                         
// double (&volume_t_k)[KM]                                 = grid_mp_volume_t_k_;
// double (&volume_t_marg_k)[KM]                            = grid_mp_volume_t_marg_k_;
                                                         
                                                         
// int &sfc_layer_type                                      = grid_mp_sfc_layer_type_;
// int &kmt_kmin                                            = grid_mp_kmt_kmin_;
// int &n_topo_smooth                                       = grid_mp_n_topo_smooth_;

//extern double (*grid_mp_bath_g_)[];

// int (*(&kmt_g))[IMT_GLOBAL]                              = grid_mp_kmt_g_;
// int (*(&basin_g))[IMT_GLOBAL]                            = grid_mp_basin_g_;

double (&dxu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_dxu_;
double (&dyu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_dyu_;
double (&dxt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_dxt_;
double (&dyt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_dyt_;
double (&dxur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_dxur_;
double (&dyur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_dyur_;
// double (&dxtr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_dxtr_;
// double (&dytr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_dytr_;
double (&hts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_hts_;
double (&htw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_htw_;
double (&hun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_hun_;
double (&hue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_hue_;
double (&ulat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_ulat_;
double (&ulon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_ulon_;
double (&tlat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_tlat_;
double (&tlon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_tlon_;
double (&angle)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = grid_mp_angle_;
double (&anglet)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = grid_mp_anglet_;
double (&fcor)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = grid_mp_fcor_;
double (&fcort)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = grid_mp_fcort_;
double (&uarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = grid_mp_uarea_;
double (&tarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = grid_mp_tarea_;
double (&uarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = grid_mp_uarea_r_;
double (&tarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = grid_mp_tarea_r_;
double (&ht)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = grid_mp_ht_;
double (&hu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = grid_mp_hu_;
// double (&hur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_hur_;

//extern double (*grid_mp_dzu_)[][][];
//extern double (*grid_mp_dzt_)[][][];

int (&basin)[MAX_BLOCKS_CLINIC][JMT][IMT]                = grid_mp_basin_;
                                                         
int (&kmt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_kmt_;
int (&kmu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_kmu_;
// int (&kmtold)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = grid_mp_kmtold_;
                                                         
double (&rcalct)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = grid_mp_rcalct_;
double (&rcalcu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = grid_mp_rcalcu_;
                                                         
int (&kmtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmtn_;
int (&kmts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmts_;
int (&kmte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmte_;
int (&kmtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmtw_;
// int (&kmun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmun_;
// int (&kmus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmus_;
// int (&kmue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmue_;
// int (&kmuw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_kmuw_;
                                                         
// int (&kmtee)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = grid_mp_kmtee_;
// int (&kmtnn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = grid_mp_kmtnn_;
                                                         
// int (*(&region_mask))[NY_BLOCK][NX_BLOCK]                = grid_mp_region_mask_;
                                                         
// int &num_regions                                         = grid_mp_num_regions_;
// int &num_ms                                              = grid_mp_num_ms_;
                                                         
int &nocean_u                                            = grid_mp_nocean_u_;
int &nocean_t                                            = grid_mp_nocean_t_;
int &nsurface_u                                          = grid_mp_nsurface_u_;
int &nsurface_t                                          = grid_mp_nsurface_t_;
                                                         
// char (&topography_filename)[CHAR_LEN]                    = grid_mp_topography_filename_;
                                                         
// int &jeq                                                 = grid_mp_jeq_;
                                                         
// double (*(&tlat_g))[IMT_GLOBAL]                          = grid_mp_tlat_g_;
// double (*(&tlon_g))[IMT_GLOBAL]                          = grid_mp_tlon_g_;
                                                         
double (&at0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_at0_;
double (&atn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_atn_;
double (&ate)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_ate_;
double (&atne)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_atne_;
// double (&au0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_au0_;
// double (&aus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_aus_;
// double (&auw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = grid_mp_auw_;
// double (&ausw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = grid_mp_ausw_;
                                                         
char (&horiz_grid_opt)[CHAR_LEN]                         = grid_mp_horiz_grid_opt_;
// char (&vert_grid_opt)[CHAR_LEN]                          = grid_mp_vert_grid_opt_;
                                                         
// char (&sfc_layer_opt)[CHAR_LEN]                          = grid_mp_sfc_layer_opt_;
// char (&topography_opt)[CHAR_LEN]                         = grid_mp_topography_opt_;
char (&horiz_grid_file)[CHAR_LEN]                        = grid_mp_horiz_grid_file_;
char (&vert_grid_file)[CHAR_LEN]                         = grid_mp_vert_grid_file_;
// char (&topography_file)[CHAR_LEN]                        = grid_mp_topography_file_;
// char (&region_mask_file)[CHAR_LEN]                       = grid_mp_region_mask_file_;
// char (&region_info_file)[CHAR_LEN]                       = grid_mp_region_info_file_;
// char (&bottom_cell_file)[CHAR_LEN]                       = grid_mp_bottom_cell_file_;
// char (&topography_outfile)[CHAR_LEN]                     = grid_mp_topography_outfile_;
char (&basin_grid_file)[CHAR_LEN]                        = grid_mp_basin_grid_file_;
                                                         
int &stdout                                              = grid_mp_stdout_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double &area_u                                           = __grid_MOD_area_u;
double &area_t                                           = __grid_MOD_area_t;
                                                         
double &volume_u                                         = __grid_MOD_volume_u;
double &volume_t                                         = __grid_MOD_volume_t;
                                                         
// double &volume_t_marg                                    = __grid_MOD_volume_t_marg;
                                                         
// double &area_t_marg                                      = __grid_MOD_area_t_marg;
                                                         
// double &uarea_equator                                    = __grid_MOD_uarea_equator;
                                                         
// double (&area_t_k)[KM]                                   = __grid_MOD_area_t_k;
                                                         
// double (&volume_t_k)[KM]                                 = __grid_MOD_volume_t_k;
// double (&volume_t_marg_k)[KM]                            = __grid_MOD_volume_t_marg_k;
                                                         
                                                         
// int &sfc_layer_type                                      = __grid_MOD_sfc_layer_type;
// int &kmt_kmin                                            = __grid_MOD_kmt_kmin;
// int &n_topo_smooth                                       = __grid_MOD_n_topo_smooth;

//extern double (*__grid_MOD_bath_g_)[];

// int (*(&kmt_g))[IMT_GLOBAL]                              = __grid_MOD_kmt_g;
// int (*(&basin_g))[IMT_GLOBAL]                            = __grid_MOD_basin_g;

double (&dxu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_dxu;
double (&dyu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_dyu;
double (&dxt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_dxt;
double (&dyt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_dyt;
double (&dxur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_dxur;
double (&dyur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_dyur;
// double (&dxtr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_dxtr;
// double (&dytr)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_dytr;
double (&hts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_hts;
double (&htw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_htw;
double (&hun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_hun;
double (&hue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_hue;
double (&ulat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_ulat;
double (&ulon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_ulon;
double (&tlat)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_tlat;
double (&tlon)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_tlon;
double (&angle)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = __grid_MOD_angle;
double (&anglet)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = __grid_MOD_anglet;
double (&fcor)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]    = __grid_MOD_fcor;
double (&fcort)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = __grid_MOD_fcort;
double (&uarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = __grid_MOD_uarea;
double (&tarea)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]   = __grid_MOD_tarea;
double (&uarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __grid_MOD_uarea_r;
double (&tarea_r)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __grid_MOD_tarea_r;
double (&ht)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = __grid_MOD_ht;
double (&hu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = __grid_MOD_hu;
// double (&hur)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_hur;

//extern double (*__grid_MOD_dzu_)[][][];
//extern double (*__grid_MOD_dzt_)[][][];

int (&basin)[MAX_BLOCKS_CLINIC][JMT][IMT]                = __grid_MOD_basin;
                                                         
int (&kmt)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_kmt;
int (&kmu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_kmu;
// int (&kmtold)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]     = __grid_MOD_kmtold;
                                                         
double (&rcalct)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = __grid_MOD_rcalct;
double (&rcalcu)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]  = __grid_MOD_rcalcu;
                                                         
int (&kmtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmtn;
int (&kmts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmts;
int (&kmte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmte;
int (&kmtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmtw;
// int (&kmun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmun;
// int (&kmus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmus;
// int (&kmue)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmue;
// int (&kmuw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_kmuw;
                                                         
// int (&kmtee)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = __grid_MOD_kmtee;
// int (&kmtnn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]      = __grid_MOD_kmtnn;
                                                         
// int (*(&region_mask))[NY_BLOCK][NX_BLOCK]                = __grid_MOD_region_mask;
                                                         
// int &num_regions                                         = __grid_MOD_num_regions;
// int &num_ms                                              = __grid_MOD_num_ms;
                                                         
int &nocean_u                                            = __grid_MOD_nocean_u;
int &nocean_t                                            = __grid_MOD_nocean_t;
int &nsurface_u                                          = __grid_MOD_nsurface_u;
int &nsurface_t                                          = __grid_MOD_nsurface_t;
                                                         
// char (&topography_filename)[CHAR_LEN]                    = __grid_MOD_topography_filename;
                                                         
// int &jeq                                                 = __grid_MOD_jeq;
                                                         
// double (*(&tlat_g))[IMT_GLOBAL]                          = __grid_MOD_tlat_g;
// double (*(&tlon_g))[IMT_GLOBAL]                          = __grid_MOD_tlon_g;
                                                         
double (&at0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_at0;
double (&atn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_atn;
double (&ate)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_ate;
double (&atne)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_atne;
// double (&au0)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_au0;
// double (&aus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_aus;
// double (&auw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]        = __grid_MOD_auw;
// double (&ausw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK]       = __grid_MOD_ausw;
                                                         
char (&horiz_grid_opt)[CHAR_LEN]                         = __grid_MOD_horiz_grid_opt;
// char (&vert_grid_opt)[CHAR_LEN]                          = __grid_MOD_vert_grid_opt;
                                                         
// char (&sfc_layer_opt)[CHAR_LEN]                          = __grid_MOD_sfc_layer_opt;
// char (&topography_opt)[CHAR_LEN]                         = __grid_MOD_topography_opt;
char (&horiz_grid_file)[CHAR_LEN]                        = __grid_MOD_horiz_grid_file;
char (&vert_grid_file)[CHAR_LEN]                         = __grid_MOD_vert_grid_file;
// char (&topography_file)[CHAR_LEN]                        = __grid_MOD_topography_file;
// char (&region_mask_file)[CHAR_LEN]                       = __grid_MOD_region_mask_file;
// char (&region_info_file)[CHAR_LEN]                       = __grid_MOD_region_info_file;
// char (&bottom_cell_file)[CHAR_LEN]                       = __grid_MOD_bottom_cell_file;
// char (&topography_outfile)[CHAR_LEN]                     = __grid_MOD_topography_outfile;
char (&basin_grid_file)[CHAR_LEN]                        = __grid_MOD_basin_grid_file;
                                                         
int &stdout                                              = __grid_MOD_stdout;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppGrid

void tgrid_to_ugrid (double (&array_ugrid)[NY_BLOCK][NX_BLOCK],
    const double (&array_tgrid)[NY_BLOCK][NX_BLOCK], const int &iblock) {
  // using CppGrid::a0;
  // using CppGrid::aus;
  // using CppGrid::auw;
  // using CppGrid::ausw;
  using CppConstantMod::C0;
  using CppConstantMod::P25;

  for (int j = 0; j < NY_BLOCK-1; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      // array_ugrid[j][i] = au0 [iblock][j][i] * array_tgrid[j  ][i  ]
      //                   + aus [iblock][j][i] * array_tgrid[j+1][i  ]
      //                   + auw [iblock][j][i] * array_tgrid[j  ][i-1]
      //                   + ausw[iblock][j][i] * array_tgrid[j+1][i-1];
      array_ugrid[j][i] = P25 * array_tgrid[j  ][i  ]
                        + P25 * array_tgrid[j+1][i  ]
                        + P25 * array_tgrid[j  ][i-1]
                        + P25 * array_tgrid[j+1][i-1];
    }
  }
  for (int i = 0; i < NX_BLOCK; ++i) {
    array_ugrid[NY_BLOCK-1][i] = C0;
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    array_ugrid[j][0] = C0;
  }
  return ;
}

void ugrid_to_tgrid(
    double (&array_tgrid)[NY_BLOCK][NX_BLOCK],
    const double (&array_ugrid)[NY_BLOCK][NX_BLOCK],
    const int &iblock, const int &k) {

  using CppGrid::      at0;
  using CppGrid::      atn;
  using CppGrid::      ate;
  using CppGrid::      atne;
  using CppPconstMod:: viv;
  using CppPconstMod:: vit;
  using CppConstantMod::C0;

  const double epsln = 1.0e-25;
  for (int j = 1; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK-1; ++i) {
      array_tgrid[j][i] = vit[iblock][k][j][i] * 
         (at0 [iblock][j][i] * array_ugrid[j  ][i  ] + 
          atn [iblock][j][i] * array_ugrid[j-1][i  ] +
          ate [iblock][j][i] * array_ugrid[j  ][i+1] +
          atne[iblock][j][i] * array_ugrid[j-1][i+1]) /
         (viv[iblock][k][j  ][i  ] * at0 [iblock][j][i] +
          viv[iblock][k][j-1][i  ] * atn [iblock][j][i] +
          viv[iblock][k][j  ][i+1] * ate [iblock][j][i] +
          viv[iblock][k][j-1][i+1] * atne[iblock][j][i] + epsln);
    }
  }
  for (int i = 0; i < NX_BLOCK; ++i) {
    array_tgrid[0][i] = C0;
  } 
  for (int j = 0; j < NY_BLOCK; ++j) {
    array_tgrid[j][NX_BLOCK-1] = C0;
  } 
  return ;
}
#endif // LICOM_ENABLE_FORTRAN
