#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "../head/cpp_buf_mod.h"
#ifdef CANUTO
#include "../head/cpp_canuto_mod.h"
#endif // CANUTO
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_grid.h"
#ifndef BIHAR
#include "../head/cpp_hmix_del2.h"
#else // BIHAR
#include "../head/cpp_hmix_del4.h"
#endif // BIHAR
#ifdef ISO
#include "../head/cpp_isopyc_mod.h"
#endif // ISO
#include "../head/cpp_output_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pmix_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/kokkos_buf_mod.h"
#ifdef CANUTO
#include "../head/kokkos_canuto_mod.h"
#endif // CANUTO
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_grid.h"
#ifndef BIHAR
#include "../head/kokkos_hmix_del2.h"
#else // BIHAR
#include "../head/kokkos_hmix_del4.h"
#endif // BIHAR
#ifdef ISO
#include "../head/kokkos_isopyc_mod.h"
#endif // ISO
#include "../head/kokkos_output_mod.h"
#include "../head/kokkos_pconst_mod.h"
#include "../head/kokkos_pmix_mod.h"
#include "../head/kokkos_work_mod.h"
#include "../head/kokkos_tracer_mod.h"

#ifdef CANUTO
#include "../head/fortran_canuto_mod.h"
#endif // CANUTO

#include "../head/cpp_tmp_var.h"
#include "../head/kokkos_tmp_var.h"

#include "../head/fortran_extern_functions.h"

#include "../head/kokkos_config.hpp"

#include "Kokkos_Core.hpp"

#include <cstdio>

using CppDomain  ::nblocks_clinic;
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

using CppParamMod::S_IMT;
using CppParamMod::S_JMT;

using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;

using CppParamMod::KMP1;

using Kokkos::deep_copy;
using Kokkos::create_mirror_view;
using Kokkos::View;

static void kokkos_init_buf_mod();
#ifdef CANUTO
static void kokkos_init_canuto_mod();
#endif // CANUTO
static void kokkos_init_dyn_mod();
static void kokkos_init_forc_mod();
static void kokkos_init_grid();
#ifndef BIHAR
static void kokkos_init_hmix_del2();
#else // BIHAR
static void kokkos_init_hmix_del4();
#endif // BIHAR
#ifdef ISO
static void kokkos_init_isopyc_mod();
#endif // ISO
static void kokkos_init_output_mod();
static void kokkos_init_pconst_mod();
static void kokkos_init_pmix_mod();
static void kokkos_init_work_mod();
static void kokkos_init_tracer_mod();
static void kokkos_init_tmp_var();

void kokkos_init_view() {

  kokkos_init_buf_mod();
#ifdef CANUTO
  kokkos_init_canuto_mod();
#endif // CANUTO
  kokkos_init_dyn_mod();
  kokkos_init_forc_mod();
  kokkos_init_grid();
#ifndef BIHAR
  kokkos_init_hmix_del2();
#else // BIHAR
  kokkos_init_hmix_del4();
#endif // BIHAR
#ifdef ISO
  kokkos_init_isopyc_mod();
#endif // ISO
  kokkos_init_output_mod();
  kokkos_init_pconst_mod();
  kokkos_init_pmix_mod();
  kokkos_init_work_mod();
  kokkos_init_tracer_mod();
  kokkos_init_tmp_var();

  return ;
}

static void kokkos_init_buf_mod() {
  using namespace CppBufMod;
  using namespace KokkosBufMod;

  // auto dev = Kokkos::DefaultExecutionSpace();

  // p_v_ifrac = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);

  // new (p_v_ifrac) ViewDouble3D("pointer_view_ifrac",
  //     MAX_BLOCKS_CLINIC, JMT, IMT); 

  // ViewDouble3D::HostMirror h_v_ifrac = create_mirror_view(*p_v_ifrac);

  // for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
  //   for (int j = 0; j < JMT; ++j) {
  //     for (int i = 0; i < IMT; ++i) {
//        h_v_ifrac(iblock, j, i)  = ifrac[iblock][j][i];
  //     }
  //   }
  // }

  // deep_copy(dev, *p_v_ifrac, h_v_ifrac);

  return ;
}
#ifdef CANUTO
static void kokkos_init_canuto_mod() {
  using namespace CppCanutoMod;
  using namespace KokkosCanutoMod;

  p_v_and2on2a1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_amtaun2a1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
                                     
  p_v_sma1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_sha1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_ssa1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
                                     
  p_v_rib  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_ridb = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
                                     
  p_v_slq2b = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_smb = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_shb = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_ssb = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
                                     
  p_v_irimax = (ViewInt1D *) malloc(sizeof(ViewInt1D));
                                     
  // p_v_sisamax = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (4 * N_THETA_R_OCT + 1));
  // p_v_ra_rmax = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (4 * N_THETA_R_OCT + 1));
  // p_v_c_y_r0  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (4 * N_THETA_R_OCT + 1));
  p_v_back_ra_r = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_sm_r1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_sh_r1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_ss_r1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_slq2_r1 = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  // p_v_ria = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * NTBL);
  // p_v_slq2a = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * NTBL);
  p_v_sma = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_sha = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  auto dev = Kokkos::DefaultExecutionSpace();
  new (p_v_and2on2a1) ViewDouble1D("pointer_view_and2on2a1", 2 * MT + 1); 
  new (p_v_amtaun2a1) ViewDouble1D("pointer_view_amtaun2a1", 2 * MT + 1); 

  new (p_v_sma1) ViewDouble1D("pointer_view_sma1", 2 * MT + 1); 
  new (p_v_sha1) ViewDouble1D("pointer_view_sha1", 2 * MT + 1); 
  new (p_v_ssa1) ViewDouble1D("pointer_view_ssa1", 2 * MT + 1); 

  new (p_v_rib) ViewDouble1D("pointer_view_rib", 2 * MT + 1); 
  new (p_v_ridb) ViewDouble1D("pointer_view_ridb", 2 * MT + 1); 

  new (p_v_slq2b) ViewDouble2D("pointer_view_slq2b", 2 * MT + 1, 2 * MT + 1); 

  new (p_v_smb) ViewDouble2D("pointer_view_smb", 2 * MT + 1, 2 * MT + 1); 
  new (p_v_shb) ViewDouble2D("pointer_view_shb", 2 * MT + 1, 2 * MT + 1); 
  new (p_v_ssb) ViewDouble2D("pointer_view_ssb", 2 * MT + 1, 2 * MT + 1); 

  new (p_v_irimax) ViewInt1D("pointer_view_irimax", 2 * MT + 1); 

  new (p_v_back_ra_r) ViewDouble1D("pointer_view_back_ra_r", 4 * N_THETA_R_OCT + 1); 

  new (p_v_sm_r1) ViewDouble1D("pointer_view_sm_r1", 4 * N_THETA_R_OCT + 1); 
  new (p_v_sh_r1) ViewDouble1D("pointer_view_sh_r1", 4 * N_THETA_R_OCT + 1); 
  new (p_v_ss_r1) ViewDouble1D("pointer_view_ss_r1", 4 * N_THETA_R_OCT + 1); 
   
  new (p_v_slq2_r1) ViewDouble1D("pointer_view_slq2_r1", 4 * N_THETA_R_OCT + 1); 

  new (p_v_sma) ViewDouble1D("pointer_view_sma", NTBL); 
  new (p_v_sha) ViewDouble1D("pointer_view_sha", NTBL); 

  UnManagedViewDouble1D h_v_and2on2a1 (&and2on2a1[0], 2 * MT + 1);
  UnManagedViewDouble1D h_v_amtaun2a1 (&amtaun2a1[0], 2 * MT + 1);

  UnManagedViewDouble1D h_v_sma1 (&sma1[0], 2 * MT + 1);
  UnManagedViewDouble1D h_v_sha1 (&sha1[0], 2 * MT + 1);
  UnManagedViewDouble1D h_v_ssa1 (&ssa1[0], 2 * MT + 1);

  UnManagedViewDouble1D h_v_rib  (&rib[0], 2 * MT + 1);
  UnManagedViewDouble1D h_v_ridb (&ridb[0], 2 * MT + 1);

  UnManagedViewDouble2D h_v_slq2b (&slq2b[0][0], 2 * MT + 1, 2 * MT + 1);

  UnManagedViewDouble2D h_v_smb (&smb[0][0], 2 * MT + 1, 2 * MT + 1);
  UnManagedViewDouble2D h_v_shb (&shb[0][0], 2 * MT + 1, 2 * MT + 1);
  UnManagedViewDouble2D h_v_ssb (&ssb[0][0], 2 * MT + 1, 2 * MT + 1);

  UnManagedViewInt1D    h_v_irimax (&irimax[0], 2 * MT + 1);
  UnManagedViewDouble1D h_v_back_ra_r (&back_ra_r[0], 4 * N_THETA_R_OCT + 1);

  UnManagedViewDouble1D h_v_sm_r1 (&sm_r1[0], 4 * N_THETA_R_OCT + 1);
  UnManagedViewDouble1D h_v_sh_r1 (&sh_r1[0], 4 * N_THETA_R_OCT + 1);
  UnManagedViewDouble1D h_v_ss_r1 (&ss_r1[0], 4 * N_THETA_R_OCT + 1);

  UnManagedViewDouble1D h_v_slq2_r1 (&slq2_r1[0], 4 * N_THETA_R_OCT + 1);

  UnManagedViewDouble1D h_v_sma (&sma[0], NTBL);
  UnManagedViewDouble1D h_v_sha (&sha[0], NTBL);

  Kokkos::deep_copy(dev, *p_v_and2on2a1, h_v_and2on2a1);
  Kokkos::deep_copy(dev, *p_v_amtaun2a1, h_v_amtaun2a1);
  Kokkos::deep_copy(dev, *p_v_sma1, h_v_sma1);
  Kokkos::deep_copy(dev, *p_v_sha1, h_v_sha1);
  Kokkos::deep_copy(dev, *p_v_ssa1, h_v_ssa1);
  Kokkos::deep_copy(dev, *p_v_rib, h_v_rib);
  Kokkos::deep_copy(dev, *p_v_ridb, h_v_ridb);
  Kokkos::deep_copy(dev, *p_v_slq2b, h_v_slq2b);
  Kokkos::deep_copy(dev, *p_v_smb, h_v_smb);
  Kokkos::deep_copy(dev, *p_v_shb, h_v_shb);
  Kokkos::deep_copy(dev, *p_v_ssb, h_v_ssb);
  Kokkos::deep_copy(dev, *p_v_irimax, h_v_irimax);
  Kokkos::deep_copy(dev, *p_v_back_ra_r, h_v_back_ra_r);
  Kokkos::deep_copy(dev, *p_v_sm_r1, h_v_sm_r1);
  Kokkos::deep_copy(dev, *p_v_sh_r1, h_v_sh_r1);
  Kokkos::deep_copy(dev, *p_v_ss_r1, h_v_ss_r1);
  Kokkos::deep_copy(dev, *p_v_slq2_r1, h_v_slq2_r1);
  Kokkos::deep_copy(dev, *p_v_sma, h_v_sma);
  Kokkos::deep_copy(dev, *p_v_sha, h_v_sha);
#else

  new (p_v_and2on2a1) ViewDouble1D(&and2on2a1[0], 2 * MT + 1); 
  new (p_v_amtaun2a1) ViewDouble1D(&amtaun2a1[0], 2 * MT + 1); 

  new (p_v_sma1) ViewDouble1D(&sma1[0], 2 * MT + 1); 
  new (p_v_sha1) ViewDouble1D(&sha1[0], 2 * MT + 1); 
  new (p_v_ssa1) ViewDouble1D(&ssa1[0], 2 * MT + 1); 

  new (p_v_rib) ViewDouble1D(&rib[0], 2 * MT + 1); 
  new (p_v_ridb) ViewDouble1D(&ridb[0], 2 * MT + 1); 

  new (p_v_slq2b) ViewDouble2D(&slq2b[0][0], 2 * MT + 1, 2 * MT + 1); 

  new (p_v_smb) ViewDouble2D(&smb[0][0], 2 * MT + 1, 2 * MT + 1); 
  new (p_v_shb) ViewDouble2D(&shb[0][0], 2 * MT + 1, 2 * MT + 1); 
  new (p_v_ssb) ViewDouble2D(&ssb[0][0], 2 * MT + 1, 2 * MT + 1); 

  new (p_v_irimax) ViewInt1D(&irimax[0], 2 * MT + 1); 

  new (p_v_back_ra_r) ViewDouble1D(&back_ra_r[0], 4 * N_THETA_R_OCT + 1); 

  new (p_v_sm_r1) ViewDouble1D(&sm_r1[0], 4 * N_THETA_R_OCT + 1); 
  new (p_v_sh_r1) ViewDouble1D(&sh_r1[0], 4 * N_THETA_R_OCT + 1); 
  new (p_v_ss_r1) ViewDouble1D(&ss_r1[0], 4 * N_THETA_R_OCT + 1); 
   
  new (p_v_slq2_r1) ViewDouble1D(&slq2_r1[0], 4 * N_THETA_R_OCT + 1); 

  new (p_v_sma) ViewDouble1D(&sma[0], NTBL); 
  new (p_v_sha) ViewDouble1D(&sha[0], NTBL); 
#endif

  return ;
}
#endif // CANUTO

static void kokkos_init_dyn_mod() {
  using namespace CppDynMod;
  using namespace KokkosDynMod;

  p_v_ub = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_vb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_ubp = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_vbp = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_h0p = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_up = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_vp = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_ws = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_h0l = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_h0f = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_h0bl = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_h0bf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_utl = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_utf = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_vtl = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_vtf = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_sbcx = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_bbcx = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_sbcy = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_bbcy = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_h0 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_u = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_v = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_gg = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_dlu = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_dlv = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_dlub = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dlvb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_ub) ViewDouble3D("pointer_view_ub",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_vb) ViewDouble3D("pointer_view_vb",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_ubp) ViewDouble3D("pointer_view_ubp",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_vbp) ViewDouble3D("pointer_view_vbp",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_h0p) ViewDouble3D("pointer_view_h0p",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_up) ViewDouble4D("pointer_view_up",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vp) ViewDouble4D("pointer_view_vp",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_ws) ViewDouble4D("pointer_view_ws",
      MAX_BLOCKS_CLINIC, KMP1, JMT, IMT); 

  new (p_v_h0l) ViewDouble3D("pointer_view_h0l",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0f) ViewDouble3D("pointer_view_h0f",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0bl) ViewDouble3D("pointer_view_h0bl",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0bf) ViewDouble3D("pointer_view_h0bf",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_utl) ViewDouble4D("pointer_view_utl",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_utf) ViewDouble4D("pointer_view_utf",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vtl) ViewDouble4D("pointer_view_vtl",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vtf) ViewDouble4D("pointer_view_vtf",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_sbcx) ViewDouble3D("pointer_view_sbcx",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_bbcx) ViewDouble3D("pointer_view_bbcx",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_sbcy) ViewDouble3D("pointer_view_sbcy",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_bbcy) ViewDouble3D("pointer_view_bbcy",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_h0) ViewDouble3D("pointer_view_h0",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_u) ViewDouble4D("pointer_view_u",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_v) ViewDouble4D("pointer_view_v",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_gg) ViewDouble4D("pointer_view_gg",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_dlu) ViewDouble4D("pointer_view_dlu",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_dlv) ViewDouble4D("pointer_view_dlv",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_dlub) ViewDouble3D("pointer_view_dlub",
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_dlvb) ViewDouble3D("pointer_view_dlvb",
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble3D h_v_ub (&ub[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_vb (&vb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_ubp (&ubp[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_vbp (&vbp[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble3D h_v_h0p (&h0p[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble4D h_v_up (&up[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  UnManagedViewDouble4D h_v_vp (&vp[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  UnManagedViewDouble4D h_v_ws (&ws[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMP1, JMT, IMT); 

  UnManagedViewDouble3D h_v_h0l (&h0l[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_h0f (&h0f[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_h0bl (&h0bl[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_h0bf (&h0bf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble4D h_v_utf (&utf[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  UnManagedViewDouble4D h_v_vtf (&vtf[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  UnManagedViewDouble4D h_v_utl (&utl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  UnManagedViewDouble4D h_v_vtl (&vtl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  UnManagedViewDouble3D h_v_sbcx (&sbcx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_bbcx (&bbcx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble3D h_v_sbcy (&sbcy[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  UnManagedViewDouble3D h_v_bbcy (&bbcy[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble3D h_v_h0 (&h0[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  UnManagedViewDouble4D h_v_u (&u[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  UnManagedViewDouble4D h_v_v (&v[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  auto dev = Kokkos::DefaultExecutionSpace();

  deep_copy(dev, *p_v_ub,   h_v_ub);
  deep_copy(dev, *p_v_vb,   h_v_vb);
  deep_copy(dev, *p_v_ubp,  h_v_ubp);
  deep_copy(dev, *p_v_vbp,  h_v_vbp);
  deep_copy(dev, *p_v_h0p,  h_v_h0p);
  deep_copy(dev, *p_v_up,   h_v_up);
  deep_copy(dev, *p_v_vp,   h_v_vp);
  deep_copy(dev, *p_v_ws,   h_v_ws);
  deep_copy(dev, *p_v_h0l,  h_v_h0l);
  deep_copy(dev, *p_v_h0f,  h_v_h0f);
  deep_copy(dev, *p_v_h0bl, h_v_h0bl);
  deep_copy(dev, *p_v_h0bf, h_v_h0bf);
  deep_copy(dev, *p_v_utl,  h_v_utl);
  deep_copy(dev, *p_v_utf,  h_v_utf);
  deep_copy(dev, *p_v_vtl,  h_v_vtl);
  deep_copy(dev, *p_v_vtf,  h_v_vtf);
  deep_copy(dev, *p_v_sbcx, h_v_sbcx);
  deep_copy(dev, *p_v_bbcx, h_v_bbcx);
  deep_copy(dev, *p_v_sbcy, h_v_sbcy);
  deep_copy(dev, *p_v_bbcy, h_v_bbcy);

  deep_copy(dev, *p_v_h0,   h_v_h0);
  deep_copy(dev, *p_v_u,    h_v_u);
  deep_copy(dev, *p_v_v,    h_v_v);
#else
  new (p_v_ub) ViewDouble3D(&ub[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_vb) ViewDouble3D(&vb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_ubp) ViewDouble3D(&ubp[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_vbp) ViewDouble3D(&vbp[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_h0p) ViewDouble3D(&h0p[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_up) ViewDouble4D(&up[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vp) ViewDouble4D(&vp[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_ws) ViewDouble4D(&ws[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMP1, JMT, IMT); 

  new (p_v_h0l) ViewDouble3D(&h0l[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0f) ViewDouble3D(&h0f[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0bl) ViewDouble3D(&h0bl[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_h0bf) ViewDouble3D(&h0bf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_utl) ViewDouble4D(&utl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_utf) ViewDouble4D(&utf[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vtl) ViewDouble4D(&vtl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_vtf) ViewDouble4D(&vtf[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_sbcx) ViewDouble3D(&sbcx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_bbcx) ViewDouble3D(&bbcx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_sbcy) ViewDouble3D(&sbcy[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_bbcy) ViewDouble3D(&bbcy[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_h0) ViewDouble3D(&h0[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

  new (p_v_u) ViewDouble4D(&u[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_v) ViewDouble4D(&v[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_gg) ViewDouble4D(&gg[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_dlu) ViewDouble4D(&dlu[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 
  new (p_v_dlv) ViewDouble4D(&dlv[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT); 

  new (p_v_dlub) ViewDouble3D(&dlub[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 
  new (p_v_dlvb) ViewDouble3D(&dlvb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT); 

#endif
  return ;
}

void kokkos_init_forc_mod() {

  using namespace CppForcMod;
  using namespace KokkosForcMod;

  p_v_su       = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_sv       = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_psa      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_tsa      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_sss      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_swv      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  // p_v_uva      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  p_v_qar      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  // p_v_cld      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  // p_v_ddd      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  // p_v_qqq      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  p_v_sst      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_nswv     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dqdt     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  // p_v_chloro   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  p_v_lwv      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_seaice   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_rain     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_snow     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_fresh    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_runoff   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_lthf     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_sshf     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ustar    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_buoytur  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_buoysol  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_restore  = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_tsf      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ssf      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
#ifdef TIDEMIX
  p_v_wave_dis = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
#endif // TIDEMIX

  p_v_buffer = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  
  p_v_s_wx   = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_s_wy   = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_s_work = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_tsa3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_wspd3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_wspdu3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_wspdv3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_psa3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_qar3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_swv3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_lwv3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_rain3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_snow3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_runoff3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_seaice3 = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_windx = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_windy = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_model_sst = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_zz = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_qs = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_theta = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_core_sensible = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_core_latent = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_core_tau = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  new (p_v_su)       ViewDouble3D("pointer_view_su",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sv)       ViewDouble3D("pointer_view_sv",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_psa)      ViewDouble3D("pointer_view_psa",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_tsa)      ViewDouble3D("pointer_view_tsa",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sss)      ViewDouble3D("pointer_view_sss",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_swv)      ViewDouble3D("pointer_view_swv",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_qar)      ViewDouble3D("pointer_view_qar",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sst)      ViewDouble3D("pointer_view_sst",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_nswv)     ViewDouble3D("pointer_view_nswv",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_dqdt)     ViewDouble3D("pointer_view_dqdt",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_lwv)      ViewDouble3D("pointer_view_lwv",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_seaice)   ViewDouble3D("pointer_view_seaice",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_rain)     ViewDouble3D("pointer_view_rain",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_snow)     ViewDouble3D("pointer_view_snow",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_fresh)    ViewDouble3D("pointer_view_fresh",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_runoff)   ViewDouble3D("pointer_view_runoff",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_lthf)     ViewDouble3D("pointer_view_lthf",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sshf)     ViewDouble3D("pointer_view_sshf",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_ustar)    ViewDouble3D("pointer_view_ustar",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_buoytur)  ViewDouble3D("pointer_view_buoytur",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_buoysol)  ViewDouble3D("pointer_view_buoysol",
      MAX_BLOCKS_CLINIC, JMT, IMT);
                     
  new (p_v_restore)  ViewDouble5D("pointer_view_restore",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_tsf)      ViewDouble3D("pointer_view_tsf",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ssf)      ViewDouble3D("pointer_view_ssf",
      MAX_BLOCKS_CLINIC, JMT, IMT);

#ifdef TIDEMIX
  new (p_v_wave_dis) ViewDouble3D("pointer_view_wave_dis",
      MAX_BLOCKS_CLINIC, JMT, IMT);
#endif // TIDEMIX

  new (p_v_buffer) ViewDouble1D("pointer_view_buffer",
      std::max(S_IMT, std::max(S_JMT, IMT_GLOBAL * JMT_GLOBAL)));

  new (p_v_s_wx)   ViewDouble2D("pointer_view_s_wx", 
      S_JMT + 2, S_IMT + 2);
  new (p_v_s_wy)   ViewDouble2D("pointer_view_s_wy", 
      S_JMT + 2, S_IMT + 2);
  new (p_v_s_work) ViewDouble2D("pointer_view_s_work", 
      S_JMT + 2, S_IMT + 2);

  new (p_v_tsa3)    ViewDouble2D("pointer_view_tsa3", JMT, IMT);
  new (p_v_wspd3)   ViewDouble2D("pointer_view_wspd3", JMT, IMT);
  new (p_v_wspdu3)  ViewDouble2D("pointer_view_wspdu3", JMT, IMT);
  new (p_v_wspdv3)  ViewDouble2D("pointer_view_wspdv3", JMT, IMT);
  new (p_v_psa3)    ViewDouble2D("pointer_view_psa3", JMT, IMT);
  new (p_v_qar3)    ViewDouble2D("pointer_view_qar3", JMT, IMT);
  new (p_v_swv3)    ViewDouble2D("pointer_view_swv3", JMT, IMT);
  new (p_v_lwv3)    ViewDouble2D("pointer_view_lwv3", JMT, IMT);
  new (p_v_rain3)   ViewDouble2D("pointer_view_rain3", JMT, IMT);
  new (p_v_snow3)   ViewDouble2D("pointer_view_snow3", JMT, IMT);
  new (p_v_runoff3) ViewDouble2D("pointer_view_runoff3", JMT, IMT);
  new (p_v_seaice3) ViewDouble2D("pointer_view_seaice3", JMT, IMT);

  new (p_v_windx)         ViewDouble2D("pointer_view_windx",         JMT, IMT);
  new (p_v_windy)         ViewDouble2D("pointer_view_windy",         JMT, IMT);
  new (p_v_model_sst)     ViewDouble2D("pointer_view_model_sst",     JMT, IMT);
  new (p_v_zz)            ViewDouble2D("pointer_view_zz",            JMT, IMT);
  new (p_v_qs)            ViewDouble2D("pointer_view_qs",            JMT, IMT);
  new (p_v_theta)         ViewDouble2D("pointer_view_theta",         JMT, IMT);
  new (p_v_core_sensible) ViewDouble2D("pointer_view_core_sensible", JMT, IMT);
  new (p_v_core_latent)   ViewDouble2D("pointer_view_core_latent",   JMT, IMT);
  new (p_v_core_tau)      ViewDouble2D("pointer_view_core_tau",      JMT, IMT);

  UnManagedViewDouble3D h_v_psa (&psa[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_tsa (&tsa[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_sss (&sss[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_swv (&swv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_qar (&qar[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_sst (&sst[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_nswv (&nswv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_dqdt (&dqdt[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_lwv (&lwv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_seaice (&seaice[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_rain (&rain[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_snow (&snow[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_fresh (&fresh[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_runoff (&runoff[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_lthf (&lthf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_sshf (&sshf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_ustar (&ustar[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_buoytur (&buoytur[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_buoysol (&buoysol[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
                     
  UnManagedViewDouble5D h_v_restore (&restore[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble3D h_v_tsf (&tsf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_ssf (&ssf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

#ifdef TIDEMIX
  UnManagedViewDouble3D h_v_wave_dis (&wave_dis[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
#endif // TIDEMIX

  auto dev = Kokkos::DefaultExecutionSpace();
  deep_copy(dev, *p_v_psa,      h_v_psa);
  deep_copy(dev, *p_v_tsa,      h_v_tsa);
  deep_copy(dev, *p_v_sss,      h_v_sss);
  deep_copy(dev, *p_v_swv,      h_v_swv);
  deep_copy(dev, *p_v_qar,      h_v_qar);
  deep_copy(dev, *p_v_sst,      h_v_sst);
  deep_copy(dev, *p_v_nswv,     h_v_nswv);
  deep_copy(dev, *p_v_dqdt,     h_v_dqdt);
  deep_copy(dev, *p_v_lwv,      h_v_lwv);
  deep_copy(dev, *p_v_seaice,   h_v_seaice);
  deep_copy(dev, *p_v_rain,     h_v_rain);
  deep_copy(dev, *p_v_snow,     h_v_snow);
  deep_copy(dev, *p_v_fresh,    h_v_fresh);
  deep_copy(dev, *p_v_runoff,   h_v_runoff);
  deep_copy(dev, *p_v_lthf,     h_v_lthf);
  deep_copy(dev, *p_v_sshf,     h_v_sshf);
  deep_copy(dev, *p_v_ustar,    h_v_ustar);
  deep_copy(dev, *p_v_buoytur,  h_v_buoytur);
  deep_copy(dev, *p_v_buoysol,  h_v_buoysol);
  deep_copy(dev, *p_v_restore,  h_v_restore);
  deep_copy(dev, *p_v_tsf,      h_v_tsf);
  deep_copy(dev, *p_v_ssf,      h_v_ssf);
#ifdef TIDEMIX
  deep_copy(dev, *p_v_wave_dis, h_v_wave_dis);
#endif // TIDEMIX

#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_su)       ViewDouble3D(&su[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sv)       ViewDouble3D(&sv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_psa)      ViewDouble3D(&psa[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_tsa)      ViewDouble3D(&tsa[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sss)      ViewDouble3D(&sss[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_swv)      ViewDouble3D(&swv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_qar)      ViewDouble3D(&qar[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sst)      ViewDouble3D(&sst[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_nswv)     ViewDouble3D(&nswv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_dqdt)     ViewDouble3D(&dqdt[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_lwv)      ViewDouble3D(&lwv[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_seaice)   ViewDouble3D(&seaice[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_rain)     ViewDouble3D(&rain[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_snow)     ViewDouble3D(&snow[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_fresh)    ViewDouble3D(&fresh[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_runoff)   ViewDouble3D(&runoff[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_lthf)     ViewDouble3D(&lthf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_sshf)     ViewDouble3D(&sshf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_ustar)    ViewDouble3D(&ustar[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_buoytur)  ViewDouble3D(&buoytur[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_buoysol)  ViewDouble3D(&buoysol[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
                     
  new (p_v_restore)  ViewDouble5D(&restore[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_tsf)      ViewDouble3D(&tsf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ssf)      ViewDouble3D(&ssf[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

#ifdef TIDEMIX
  new (p_v_wave_dis) ViewDouble3D(&wave_dis[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
#endif // TIDEMIX

  CppForcMod::buffer = new double[std::max(S_IMT, std::max(S_JMT, 
      IMT_GLOBAL * JMT_GLOBAL))];
  new (p_v_buffer) ViewDouble1D(buffer,
      std::max(S_IMT, std::max(S_JMT, IMT_GLOBAL * JMT_GLOBAL)));

  CppForcMod::tsa3    = new double [JMT * IMT];
  CppForcMod::wspd3   = new double [JMT * IMT];
  CppForcMod::wspdu3  = new double [JMT * IMT];
  CppForcMod::wspdv3  = new double [JMT * IMT];
  CppForcMod::psa3    = new double [JMT * IMT];
  CppForcMod::qar3    = new double [JMT * IMT];
  CppForcMod::swv3    = new double [JMT * IMT];
  CppForcMod::lwv3    = new double [JMT * IMT];
  CppForcMod::rain3   = new double [JMT * IMT];
  CppForcMod::snow3   = new double [JMT * IMT];
  CppForcMod::runoff3 = new double [JMT * IMT];
  CppForcMod::seaice3 = new double [JMT * IMT];

  new (p_v_tsa3)    ViewDouble2D(tsa3,    JMT, IMT);
  new (p_v_wspd3)   ViewDouble2D(wspd3,   JMT, IMT);
  new (p_v_wspdu3)  ViewDouble2D(wspdu3,  JMT, IMT);
  new (p_v_wspdv3)  ViewDouble2D(wspdv3,  JMT, IMT);
  new (p_v_psa3)    ViewDouble2D(psa3,    JMT, IMT);
  new (p_v_qar3)    ViewDouble2D(qar3,    JMT, IMT);
  new (p_v_swv3)    ViewDouble2D(swv3,    JMT, IMT);
  new (p_v_lwv3)    ViewDouble2D(lwv3,    JMT, IMT);
  new (p_v_rain3)   ViewDouble2D(rain3,   JMT, IMT);
  new (p_v_snow3)   ViewDouble2D(snow3,   JMT, IMT);
  new (p_v_runoff3) ViewDouble2D(runoff3, JMT, IMT);
  new (p_v_seaice3) ViewDouble2D(seaice3, JMT, IMT);

  CppForcMod::s_wx   = new double [(S_JMT+2) * (S_IMT+2)];
  CppForcMod::s_wy   = new double [(S_JMT+2) * (S_IMT+2)];
  CppForcMod::s_work = new double [(S_JMT+2) * (S_IMT+2)];

  new (p_v_s_wx)   ViewDouble2D(s_wx,   S_JMT + 2, S_IMT + 2);
  new (p_v_s_wy)   ViewDouble2D(s_wy,   S_JMT + 2, S_IMT + 2);
  new (p_v_s_work) ViewDouble2D(s_work, S_JMT + 2, S_IMT + 2);

  CppForcMod::windx         = new double [JMT * IMT];
  CppForcMod::windy         = new double [JMT * IMT];
  CppForcMod::model_sst     = new double [JMT * IMT];
  CppForcMod::zz            = new double [JMT * IMT];
  CppForcMod::qs            = new double [JMT * IMT];
  CppForcMod::theta         = new double [JMT * IMT];
  CppForcMod::core_sensible = new double [JMT * IMT];
  CppForcMod::core_latent   = new double [JMT * IMT];
  CppForcMod::core_tau      = new double [JMT * IMT];

  new (p_v_windx)         ViewDouble2D(windx,         JMT, IMT);
  new (p_v_windy)         ViewDouble2D(windy,         JMT, IMT);
  new (p_v_model_sst)     ViewDouble2D(model_sst,     JMT, IMT);
  new (p_v_zz)            ViewDouble2D(zz,            JMT, IMT);
  new (p_v_qs)            ViewDouble2D(qs,            JMT, IMT);
  new (p_v_theta)         ViewDouble2D(theta,         JMT, IMT);
  new (p_v_core_sensible) ViewDouble2D(core_sensible, JMT, IMT);
  new (p_v_core_latent)   ViewDouble2D(core_latent,   JMT, IMT);
  new (p_v_core_tau)      ViewDouble2D(core_tau,      JMT, IMT);

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  // Free memory in Fortran
  free_forc_mod_();
  return ;

}
static void kokkos_init_grid() {
  using namespace CppGrid;
  using namespace KokkosGrid;

  // p_v_area_t_k        = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KM);
  // p_v_volume_t_k      = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KM);
  // p_v_volume_t_marg_k = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KM);

  // p_v_kmt_g   = (ViewInt2D *) malloc(sizeof(ViewInt2D) *
  //     JMT_GLOBAL * IMT_GLOBAL);
  // p_v_basin_g = (ViewInt2D *) malloc(sizeof(ViewInt2D) *
  //     JMT_GLOBAL * IMT_GLOBAL);

  p_v_dxu     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dyu     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  // p_v_dxt     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_dyt     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  p_v_dxur    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dyur    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_dxyur   = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // p_v_dxtr    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_dytr    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  p_v_hts     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_htw     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_hun     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_hue     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_h_tu_swen = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
//   p_v_ulat    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_ulon    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  p_v_tlat    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
//   p_v_tlon    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_angle   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_anglet  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  p_v_fcor    = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_fcort   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_uarea   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_tarea   = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_uarea_r = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_tarea_r = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
//   p_v_ht      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_hu      = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_hur     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);

  // p_v_basin = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);

  p_v_kmt = (ViewInt3D *) malloc(sizeof(ViewInt3D));

  p_v_kmu = (ViewInt3D *) malloc(sizeof(ViewInt3D));
  // p_v_kmtold = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
                                                         
  // p_v_rcalct = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_rcalcu = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
                                                         
  p_v_kmtn = (ViewInt3D *) malloc(sizeof(ViewInt3D));
  p_v_kmts = (ViewInt3D *) malloc(sizeof(ViewInt3D));
  p_v_kmte = (ViewInt3D *) malloc(sizeof(ViewInt3D));
  p_v_kmtw = (ViewInt3D *) malloc(sizeof(ViewInt3D));
//   p_v_kmun = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_kmus = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_kmue = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_kmuw = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);

  p_v_kmt_nsew = (ViewInt4D *) malloc(sizeof(ViewInt4D));
                                                         
  // p_v_kmtee = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
  // p_v_kmtnn = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
  //     MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
                                                         
  // p_v_region_mask = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
      // MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);

  // p_v_tlat_g = (ViewDouble2D *) malloc(sizeof(ViewDouble2D) *
  //     JMT_GLOBAL * IMT_GLOBAL);
  // p_v_tlon_g = (ViewDouble2D *) malloc(sizeof(ViewDouble2D) *
  //     JMT_GLOBAL * IMT_GLOBAL);
                                                         
  p_v_at0  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_atn  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ate  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_atne = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
//   p_v_au0  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_aus  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_auw  = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);
//   p_v_ausw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
//       MAX_BLOCKS_CLINIC * NY_BLOCK * NX_BLOCK);

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  auto dev = Kokkos::DefaultExecutionSpace();

  new (p_v_dxu)     ViewDouble3D("pointer_view_dxu",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dyu)     ViewDouble3D("pointer_view_dyu",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dxur)    ViewDouble3D("pointer_view_dxur",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dyur)    ViewDouble3D("pointer_view_dyur",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dxyur)   ViewDouble4D("pointer_view_dxyur",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 2);
  new (p_v_hts)     ViewDouble3D("pointer_view_hts",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_htw)     ViewDouble3D("pointer_view_htw",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_hun)     ViewDouble3D("pointer_view_hun",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_hue)     ViewDouble3D("pointer_view_hue",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_h_tu_swen) ViewDouble5D("pointer_view_h_tu_swen",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 2, 2);
  new (p_v_tlat)    ViewDouble3D("pointer_view_tlat",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_fcor)    ViewDouble3D("pointer_view_fcor",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_fcort)   ViewDouble3D("pointer_view_fcort",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_uarea)   ViewDouble3D("pointer_view_uarea",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_tarea)   ViewDouble3D("pointer_view_tarea",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_uarea_r) ViewDouble3D("pointer_view_uarea_r",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_tarea_r) ViewDouble3D("pointer_view_uarea_r",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmt)    ViewInt3D("pointer_view_kmt",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmu)    ViewInt3D("pointer_view_kmu",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmtn) ViewInt3D("pointer_view_kmtn",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmts) ViewInt3D("pointer_view_kmts",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmte) ViewInt3D("pointer_view_kmte",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmtw) ViewInt3D("pointer_view_kmtw",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmt_nsew) ViewInt4D("pointer_view_kmt_nsew",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 4);
  new (p_v_at0)  ViewDouble3D("pointer_view_at0",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_atn)  ViewDouble3D("pointer_view_atn",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_ate)  ViewDouble3D("pointer_view_ate",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_atne) ViewDouble3D("pointer_view_atne",
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

 UnManagedViewDouble3D h_v_dxu (&dxu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_dyu (&dyu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_dxur (&dxur[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_dyur (&dyur[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_hts (&hts[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_htw (&htw[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_hun (&hun[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_hue (&hue[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_tlat (&tlat[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_fcor (&fcor[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_fcort (&fcort[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_uarea (&uarea[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_tarea (&tarea[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_uarea_r (&uarea_r[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_tarea_r (&uarea_r[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_basin (&basin[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);
 UnManagedViewInt3D h_v_kmt (&kmt[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_kmu (&kmu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_kmtn (&kmtn[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_kmts (&kmts[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_kmte (&kmte[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewInt3D h_v_kmtw (&kmtw[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_at0 (&at0[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_atn (&atn[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_ate (&ate[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
 UnManagedViewDouble3D h_v_atne (&atne[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

  ViewDouble4D::HostMirror h_v_dxyur     = create_mirror_view(*p_v_dxyur);
  ViewDouble5D::HostMirror h_v_h_tu_swen = create_mirror_view(*p_v_h_tu_swen);
  ViewInt4D::HostMirror    h_v_kmt_nsew  = create_mirror_view(*p_v_kmt_nsew);

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for(int j = 0; j < NY_BLOCK; ++j) {
      for(int i = 0; i < NX_BLOCK; ++i) {
        h_v_dxyur(iblock, j, i, 0) = dxur[iblock][j][i];
        h_v_dxyur(iblock, j, i, 1) = dyur[iblock][j][i];

        // htw hun hts hue
        h_v_h_tu_swen(iblock, j, i, 0, 0) = htw[iblock][j][i];
        h_v_h_tu_swen(iblock, j, i, 0, 1) = hun[iblock][j][i];
        h_v_h_tu_swen(iblock, j, i, 1, 0) = hts[iblock][j][i];
        h_v_h_tu_swen(iblock, j, i, 1, 1) = hue[iblock][j][i];
      }
    }
  }//end for

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        h_v_kmt_nsew(iblock, j, i, 0) = kmtn[iblock][j][i];
        h_v_kmt_nsew(iblock, j, i, 1) = kmts[iblock][j][i];
        h_v_kmt_nsew(iblock, j, i, 2) = kmte[iblock][j][i];
        h_v_kmt_nsew(iblock, j, i, 3) = kmtw[iblock][j][i];
      }
    }
  }// end for

  deep_copy(dev, *p_v_dxu,       h_v_dxu);
  deep_copy(dev, *p_v_dyu,       h_v_dyu);
  deep_copy(dev, *p_v_dxur,      h_v_dxur);
  deep_copy(dev, *p_v_dyur,      h_v_dyur);
  deep_copy(dev, *p_v_dxyur,     h_v_dxyur);
  deep_copy(dev, *p_v_hts,       h_v_hts);
  deep_copy(dev, *p_v_htw,       h_v_htw);
  deep_copy(dev, *p_v_hun,       h_v_hun);
  deep_copy(dev, *p_v_hue,       h_v_hue);
  deep_copy(dev, *p_v_h_tu_swen, h_v_h_tu_swen);
  deep_copy(dev, *p_v_tlat,      h_v_tlat);
  deep_copy(dev, *p_v_fcor,      h_v_fcor);
  deep_copy(dev, *p_v_fcort,     h_v_fcort);
  deep_copy(dev, *p_v_uarea,     h_v_uarea);
  deep_copy(dev, *p_v_tarea,     h_v_tarea);
  deep_copy(dev, *p_v_uarea_r,   h_v_uarea_r);
  deep_copy(dev, *p_v_tarea_r,   h_v_tarea_r);
  deep_copy(dev, *p_v_kmt,       h_v_kmt);
  deep_copy(dev, *p_v_kmu,       h_v_kmu);
  deep_copy(dev, *p_v_kmtn,      h_v_kmtn);
  deep_copy(dev, *p_v_kmts,      h_v_kmts);
  deep_copy(dev, *p_v_kmte,      h_v_kmte);
  deep_copy(dev, *p_v_kmtw,      h_v_kmtw);
  deep_copy(dev, *p_v_kmt_nsew,  h_v_kmt_nsew);
  deep_copy(dev, *p_v_at0,       h_v_at0);
  deep_copy(dev, *p_v_atn,       h_v_atn); 
  deep_copy(dev, *p_v_ate,       h_v_ate);
  deep_copy(dev, *p_v_atne,      h_v_atne);

#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_dxu)     ViewDouble3D(&dxu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dyu)     ViewDouble3D(&dyu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dxur)    ViewDouble3D(&dxur[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_dyur)    ViewDouble3D(&dyur[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

  dxyur = new double [MAX_BLOCKS_CLINIC * 2 * NY_BLOCK * NX_BLOCK];
  new (p_v_dxyur)   ViewDouble4D(dxyur,
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 2);

  new (p_v_hts)     ViewDouble3D(&hts[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_htw)     ViewDouble3D(&htw[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_hun)     ViewDouble3D(&hun[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_hue)     ViewDouble3D(&hue[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

  h_tu_swen = new double [MAX_BLOCKS_CLINIC * 2 * 2 * NY_BLOCK * NX_BLOCK];
  new (p_v_h_tu_swen) ViewDouble5D(h_tu_swen,
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 2, 2);

  new (p_v_tlat)    ViewDouble3D(&tlat[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_fcor)    ViewDouble3D(&fcor[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_fcort)   ViewDouble3D(&fcort[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_uarea)   ViewDouble3D(&uarea[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_tarea)   ViewDouble3D(&tarea[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_uarea_r) ViewDouble3D(&uarea_r[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_tarea_r) ViewDouble3D(&tarea_r[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmt)    ViewInt3D(&kmt[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmu)    ViewInt3D(&kmu[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmtn) ViewInt3D(&kmtn[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmts) ViewInt3D(&kmts[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmte) ViewInt3D(&kmte[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_kmtw) ViewInt3D(&kmtw[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

  kmt_nsew = new int [MAX_BLOCKS_CLINIC * 4 * NY_BLOCK * NX_BLOCK];
  new (p_v_kmt_nsew) ViewInt4D(kmt_nsew,
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK, 4);

  new (p_v_at0)  ViewDouble3D(&at0[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_atn)  ViewDouble3D(&atn[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_ate)  ViewDouble3D(&ate[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);
  new (p_v_atne) ViewDouble3D(&atne[0][0][0],
      MAX_BLOCKS_CLINIC, NY_BLOCK, NX_BLOCK);

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for(int j = 0; j < NY_BLOCK; ++j) {
      for(int i = 0; i < NX_BLOCK; ++i) {
        (*p_v_dxyur)(iblock, j, i, 0) = dxur[iblock][j][i];
        (*p_v_dxyur)(iblock, j, i, 1) = dyur[iblock][j][i];

        // htw hun hts hue
        (*p_v_h_tu_swen)(iblock, j, i, 0, 0) = htw[iblock][j][i];
        (*p_v_h_tu_swen)(iblock, j, i, 0, 1) = hun[iblock][j][i];
        (*p_v_h_tu_swen)(iblock, j, i, 1, 0) = hts[iblock][j][i];
        (*p_v_h_tu_swen)(iblock, j, i, 1, 1) = hue[iblock][j][i];
      }
    }
  }//end for

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        (*p_v_kmt_nsew)(iblock, j, i, 0) = kmtn[iblock][j][i];
        (*p_v_kmt_nsew)(iblock, j, i, 1) = kmts[iblock][j][i];
        (*p_v_kmt_nsew)(iblock, j, i, 2) = kmte[iblock][j][i];
        (*p_v_kmt_nsew)(iblock, j, i, 3) = kmtw[iblock][j][i];
      }
    }
  }// end for

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  free_grid_();

  return ;
}

#ifndef BIHAR
static void kokkos_init_hmix_del2() {
  using namespace CppHmixDel2;
  using namespace KokkosHmixDel2;

  auto dev = Kokkos::DefaultExecutionSpace();
  p_v_dtn = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dts = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dte = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dtw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_duc = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dun = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dus = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_due = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_duw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmc = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmn = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dms = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dme = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dum = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ahf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_amf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  new (p_v_dtn) ViewDouble3D("pointer_view_dtn", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dts) ViewDouble3D("pointer_view_dts", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dte) ViewDouble3D("pointer_view_dte", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dtw) ViewDouble3D("pointer_view_dtw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duc) ViewDouble3D("pointer_view_duc", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dun) ViewDouble3D("pointer_view_dun", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dus) ViewDouble3D("pointer_view_dus", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_due) ViewDouble3D("pointer_view_due", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duw) ViewDouble3D("pointer_view_duw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmc) ViewDouble3D("pointer_view_dmc", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmn) ViewDouble3D("pointer_view_dmn", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dms) ViewDouble3D("pointer_view_dms", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dme) ViewDouble3D("pointer_view_dme", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmw) ViewDouble3D("pointer_view_dmw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dum) ViewDouble3D("pointer_view_dum", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_ahf) ViewDouble3D("pointer_view_ahf", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_amf) ViewDouble3D("pointer_view_amf", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 

  ViewDouble3D::HostMirror h_v_dtn = create_mirror_view(*p_v_dtn);
  ViewDouble3D::HostMirror h_v_dts = create_mirror_view(*p_v_dts);
  ViewDouble3D::HostMirror h_v_dte = create_mirror_view(*p_v_dte);
  ViewDouble3D::HostMirror h_v_dtw = create_mirror_view(*p_v_dtw);
  ViewDouble3D::HostMirror h_v_duc = create_mirror_view(*p_v_duc);
  ViewDouble3D::HostMirror h_v_dun = create_mirror_view(*p_v_dun);
  ViewDouble3D::HostMirror h_v_dus = create_mirror_view(*p_v_dus);
  ViewDouble3D::HostMirror h_v_due = create_mirror_view(*p_v_due);
  ViewDouble3D::HostMirror h_v_duw = create_mirror_view(*p_v_duw);
  ViewDouble3D::HostMirror h_v_dmc = create_mirror_view(*p_v_dmc);
  ViewDouble3D::HostMirror h_v_dmn = create_mirror_view(*p_v_dmn);
  ViewDouble3D::HostMirror h_v_dms = create_mirror_view(*p_v_dms);
  ViewDouble3D::HostMirror h_v_dme = create_mirror_view(*p_v_dme);
  ViewDouble3D::HostMirror h_v_dmw = create_mirror_view(*p_v_dmw);
  ViewDouble3D::HostMirror h_v_dum = create_mirror_view(*p_v_dum);
  ViewDouble3D::HostMirror h_v_ahf = create_mirror_view(*p_v_ahf);
  ViewDouble3D::HostMirror h_v_amf = create_mirror_view(*p_v_amf);

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        h_v_dtn(iblock, j, i) = dtn[iblock][j][i];
        h_v_dts(iblock, j, i) = dts[iblock][j][i];
        h_v_dte(iblock, j, i) = dte[iblock][j][i];
        h_v_dtw(iblock, j, i) = dtw[iblock][j][i];
        h_v_duc(iblock, j, i) = duc[iblock][j][i];
        h_v_dun(iblock, j, i) = dun[iblock][j][i];
        h_v_dus(iblock, j, i) = dus[iblock][j][i];
        h_v_due(iblock, j, i) = due[iblock][j][i];
        h_v_duw(iblock, j, i) = duw[iblock][j][i];
        h_v_dmc(iblock, j, i) = dmc[iblock][j][i];
        h_v_dmn(iblock, j, i) = dmn[iblock][j][i];
        h_v_dms(iblock, j, i) = dms[iblock][j][i];
        h_v_dme(iblock, j, i) = dme[iblock][j][i];
        h_v_dmw(iblock, j, i) = dmw[iblock][j][i];
        h_v_dum(iblock, j, i) = dum[iblock][j][i];
        h_v_ahf(iblock, j, i) = ahf[iblock][j][i];
        h_v_amf(iblock, j, i) = amf[iblock][j][i];
      }
    }
  }

  deep_copy(dev, *p_v_dtn, h_v_dtn);
  deep_copy(dev, *p_v_dts, h_v_dts);
  deep_copy(dev, *p_v_dte, h_v_dte);
  deep_copy(dev, *p_v_dtw, h_v_dtw);
  deep_copy(dev, *p_v_duc, h_v_duc);
  deep_copy(dev, *p_v_dun, h_v_dun);
  deep_copy(dev, *p_v_dus, h_v_dus);
  deep_copy(dev, *p_v_due, h_v_due);
  deep_copy(dev, *p_v_duw, h_v_duw);
  deep_copy(dev, *p_v_dmc, h_v_dmc);
  deep_copy(dev, *p_v_dmn, h_v_dmn);
  deep_copy(dev, *p_v_dms, h_v_dms);
  deep_copy(dev, *p_v_dme, h_v_dme);
  deep_copy(dev, *p_v_dmw, h_v_dmw);
  deep_copy(dev, *p_v_dum, h_v_dum);
  deep_copy(dev, *p_v_ahf, h_v_ahf);
  deep_copy(dev, *p_v_amf, h_v_amf);
  return ;
}
#else // BIHAR
static void kokkos_init_hmix_del4() {
  using namespace CppHmixDel4;
  using namespace KokkosHmixDel4;

  p_v_dtn = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dts = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dte = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dtw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_duc = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dun = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dus = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_due = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_duw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmc = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmn = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dms = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dme = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dmw = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dum = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ahf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_amf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  // p_v_ratio_dxy = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     NX_BLOCK * NY_BLOCK);

  p_v_dt_nsew = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_du_cnsewm = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_dm_cnsew  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_dtn) ViewDouble3D("pointer_view_dtn", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dts) ViewDouble3D("pointer_view_dts", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dte) ViewDouble3D("pointer_view_dte", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dtw) ViewDouble3D("pointer_view_dtw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duc) ViewDouble3D("pointer_view_duc", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dun) ViewDouble3D("pointer_view_dun", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dus) ViewDouble3D("pointer_view_dus", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_due) ViewDouble3D("pointer_view_due", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duw) ViewDouble3D("pointer_view_duw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmc) ViewDouble3D("pointer_view_dmc", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmn) ViewDouble3D("pointer_view_dmn", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dms) ViewDouble3D("pointer_view_dms", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dme) ViewDouble3D("pointer_view_dme", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmw) ViewDouble3D("pointer_view_dmw", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dum) ViewDouble3D("pointer_view_dum", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_ahf) ViewDouble3D("pointer_view_ahf", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_amf) ViewDouble3D("pointer_view_amf", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 

  new (p_v_dt_nsew) ViewDouble4D("pointer_view_dt_nsew", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 4); 
  new (p_v_du_cnsewm) ViewDouble4D("pointer_view_du_cnsewm", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 6); 
  new (p_v_dm_cnsew) ViewDouble4D("pointer_view_dm_cnsew", 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 5); 

  UnManagedViewDouble3D h_v_dtn (&dtn[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dts (&dts[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dte (&dte[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dtw (&dtw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_duc (&duc[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dun (&dun[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dus (&dus[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_due (&due[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_duw (&duw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dmc (&dmc[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dmn (&dmn[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dms (&dms[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dme (&dme[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dmw (&dmw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_dum (&dum[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_ahf (&ahf[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  UnManagedViewDouble3D h_v_amf (&amf[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 

  ViewDouble4D::HostMirror h_v_dt_nsew   = create_mirror_view(*p_v_dt_nsew);
  ViewDouble4D::HostMirror h_v_du_cnsewm = create_mirror_view(*p_v_du_cnsewm);
  ViewDouble4D::HostMirror h_v_dm_cnsew  = create_mirror_view(*p_v_dm_cnsew);
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {

        h_v_dt_nsew(iblock, j, i, 0) = dtn[iblock][j][i];
        h_v_dt_nsew(iblock, j, i, 1) = dts[iblock][j][i];
        h_v_dt_nsew(iblock, j, i, 2) = dte[iblock][j][i];
        h_v_dt_nsew(iblock, j, i, 3) = dtw[iblock][j][i];

        h_v_du_cnsewm(iblock, j, i, 0) = duc[iblock][j][i];
        h_v_du_cnsewm(iblock, j, i, 1) = dun[iblock][j][i];
        h_v_du_cnsewm(iblock, j, i, 2) = dus[iblock][j][i];
        h_v_du_cnsewm(iblock, j, i, 3) = due[iblock][j][i];
        h_v_du_cnsewm(iblock, j, i, 4) = duw[iblock][j][i];
        h_v_du_cnsewm(iblock, j, i, 5) = dum[iblock][j][i];

        h_v_dm_cnsew(iblock, j, i, 0) = dmc[iblock][j][i];
        h_v_dm_cnsew(iblock, j, i, 1) = dmn[iblock][j][i];
        h_v_dm_cnsew(iblock, j, i, 2) = dms[iblock][j][i];
        h_v_dm_cnsew(iblock, j, i, 3) = dme[iblock][j][i];
        h_v_dm_cnsew(iblock, j, i, 4) = dmw[iblock][j][i];
      }
    }
  }
  auto dev = Kokkos::DefaultExecutionSpace();

  deep_copy(dev, *p_v_dtn, h_v_dtn);
  deep_copy(dev, *p_v_dts, h_v_dts);
  deep_copy(dev, *p_v_dte, h_v_dte);
  deep_copy(dev, *p_v_dtw, h_v_dtw);
  deep_copy(dev, *p_v_duc, h_v_duc);
  deep_copy(dev, *p_v_dun, h_v_dun);
  deep_copy(dev, *p_v_dus, h_v_dus);
  deep_copy(dev, *p_v_due, h_v_due);
  deep_copy(dev, *p_v_duw, h_v_duw);
  deep_copy(dev, *p_v_dmc, h_v_dmc);
  deep_copy(dev, *p_v_dmn, h_v_dmn);
  deep_copy(dev, *p_v_dms, h_v_dms);
  deep_copy(dev, *p_v_dme, h_v_dme);
  deep_copy(dev, *p_v_dmw, h_v_dmw);
  deep_copy(dev, *p_v_dum, h_v_dum);
  deep_copy(dev, *p_v_ahf, h_v_ahf);
  deep_copy(dev, *p_v_amf, h_v_amf);

  deep_copy(dev, *p_v_dt_nsew,   h_v_dt_nsew);
  deep_copy(dev, *p_v_du_cnsewm, h_v_du_cnsewm);
  deep_copy(dev, *p_v_dm_cnsew,  h_v_dm_cnsew);
#else  // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_dtn) ViewDouble3D(&dtn[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dts) ViewDouble3D(&dts[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dte) ViewDouble3D(&dte[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dtw) ViewDouble3D(&dtw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duc) ViewDouble3D(&duc[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dun) ViewDouble3D(&dun[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dus) ViewDouble3D(&dus[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_due) ViewDouble3D(&due[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_duw) ViewDouble3D(&duw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmc) ViewDouble3D(&dmc[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmn) ViewDouble3D(&dmn[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dms) ViewDouble3D(&dms[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dme) ViewDouble3D(&dme[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dmw) ViewDouble3D(&dmw[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_dum) ViewDouble3D(&dum[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_ahf) ViewDouble3D(&ahf[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 
  new (p_v_amf) ViewDouble3D(&amf[0][0][0], 
      nblocks_clinic, NY_BLOCK, NX_BLOCK); 

  dt_nsew   = new double[nblocks_clinic * 4 * NY_BLOCK * NX_BLOCK];
  du_cnsewm = new double[nblocks_clinic * 6 * NY_BLOCK * NX_BLOCK];
  dm_cnsew  = new double[nblocks_clinic * 5 * NY_BLOCK * NX_BLOCK];

  new (p_v_dt_nsew) ViewDouble4D(dt_nsew, 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 4); 
  new (p_v_du_cnsewm) ViewDouble4D(du_cnsewm, 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 6); 
  new (p_v_dm_cnsew) ViewDouble4D(dm_cnsew, 
      nblocks_clinic, NY_BLOCK, NX_BLOCK, 5); 

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        (*p_v_dt_nsew)(iblock, j, i, 0) = dtn[iblock][j][i];
        (*p_v_dt_nsew)(iblock, j, i, 1) = dts[iblock][j][i];
        (*p_v_dt_nsew)(iblock, j, i, 2) = dte[iblock][j][i];
        (*p_v_dt_nsew)(iblock, j, i, 3) = dtw[iblock][j][i];

        (*p_v_du_cnsewm)(iblock, j, i, 0) = duc[iblock][j][i];
        (*p_v_du_cnsewm)(iblock, j, i, 1) = dun[iblock][j][i];
        (*p_v_du_cnsewm)(iblock, j, i, 2) = dus[iblock][j][i];
        (*p_v_du_cnsewm)(iblock, j, i, 3) = due[iblock][j][i];
        (*p_v_du_cnsewm)(iblock, j, i, 4) = duw[iblock][j][i];
        (*p_v_du_cnsewm)(iblock, j, i, 5) = dum[iblock][j][i];

        (*p_v_dm_cnsew)(iblock, j, i, 0) = dmc[iblock][j][i];
        (*p_v_dm_cnsew)(iblock, j, i, 1) = dmn[iblock][j][i];
        (*p_v_dm_cnsew)(iblock, j, i, 2) = dms[iblock][j][i];
        (*p_v_dm_cnsew)(iblock, j, i, 3) = dme[iblock][j][i];
        (*p_v_dm_cnsew)(iblock, j, i, 4) = dmw[iblock][j][i];
      }
    }
  }
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  return ;
}
#endif // BIHAR

// #ifdef ISO
// static void kokkos_init_isopyc_mod() {
//   using namespace CppIsopycMod;
//   using namespace KokkosIsopycMod;
//   auto dev = Kokkos::DefaultExecutionSpace();
//   p_v_dptlim = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (NRPL+1));

//   p_v_fzisop = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KM);

//   p_v_ahisop = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       IMT * JMT * MAX_BLOCKS_CLINIC);
//   p_v_athkdf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       IMT * JMT * MAX_BLOCKS_CLINIC);

//   p_v_e = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) * 
//       IMT * KMP1 * JMT * 3 * MAX_BLOCKS_CLINIC);

//   p_v_rhoi = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) * 
//       IMT * (KM+1) * JMT * NRPL * MAX_BLOCKS_CLINIC);

//   p_v_k1 = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) * 
//       IMT * (KM+1) * JMT * 1 * MAX_BLOCKS_CLINIC);
//   p_v_k2 = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) * 
//       IMT * (KM+1) * JMT * 1 * MAX_BLOCKS_CLINIC);
//   p_v_k3 = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) * 
//       IMT * (KM+1) * JMT * 3 * MAX_BLOCKS_CLINIC);

//   p_v_adv_vetiso = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) * 
//       IMT * KM * JMT * MAX_BLOCKS_CLINIC);
//   p_v_adv_vbtiso = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) * 
//       IMT * (KM+1) * JMT * MAX_BLOCKS_CLINIC);
//   p_v_adv_vntiso = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) * 
//       IMT * KM * JMT * MAX_BLOCKS_CLINIC);

// #ifdef isopycmixspatialvar
//   p_v_dciso1 = (ViewFloat4D *) malloc(sizeof(ViewFloat4D) * 
//       IMT * KM * JMT * MAX_BLOCKS_CLINIC);
//   p_v_dciso2 = (ViewFloat4D *) malloc(sizeof(ViewFloat4D) * 
//       IMT * KM * JMT * MAX_BLOCKS_CLINIC);
// #endif // isopycmixspatialvar

//   p_v_kisrpl = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KM);
//   p_v_krplin = (ViewInt1D *) malloc(sizeof(ViewInt1D) * NRPL);

//   p_v_zt   = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) *  KM);
//   p_v_dzw  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (KM+1));
//   p_v_dzwr = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * (KM+1));
//   p_v_dzr  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) *  KM);

//   p_v_tmask = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) * 
//       IMT * KM * JMT * MAX_BLOCKS_CLINIC);

//   p_v_f3 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
//       IMT * JMT * MAX_BLOCKS_CLINIC);

//   new (p_v_dptlim) ViewDouble1D("pointer_view_dptlim", NRPL+1); 

//   new (p_v_fzisop) ViewDouble1D("pointer_view_fzisop", KM); 

//   new (p_v_ahisop) ViewDouble3D("pointer_view_ahisop", MAX_BLOCKS_CLINIC, JMT, IMT); 
//   new (p_v_athkdf) ViewDouble3D("pointer_view_athkdf", MAX_BLOCKS_CLINIC, JMT, IMT); 

//   new (p_v_e) ViewDouble5D("pointer_view_e", MAX_BLOCKS_CLINIC, 3, JMT, KMP1, IMT); 
//   new (p_v_rhoi) ViewDouble5D("pointer_view_rhoi", MAX_BLOCKS_CLINIC, NRPL, JMT, KM+1, IMT); 

//   new (p_v_k1) ViewDouble5D("pointer_view_k1", MAX_BLOCKS_CLINIC, 1, JMT, KM+1, IMT); 
//   new (p_v_k2) ViewDouble5D("pointer_view_k2", MAX_BLOCKS_CLINIC, 1, JMT, KM+1, IMT); 
//   new (p_v_k3) ViewDouble5D("pointer_view_k3", MAX_BLOCKS_CLINIC, 3, JMT, KM+1, IMT); 

//   new (p_v_adv_vetiso) ViewDouble4D("pointer_view_adv_vetiso", MAX_BLOCKS_CLINIC, JMT, KM,   IMT); 
//   new (p_v_adv_vbtiso) ViewDouble4D("pointer_view_adv_vbtiso", MAX_BLOCKS_CLINIC, JMT, KM+1, IMT); 
//   new (p_v_adv_vntiso) ViewDouble4D("pointer_view_adv_vntiso", MAX_BLOCKS_CLINIC, JMT, KM,   IMT); 

// #ifdef isopycmixspatialvar
//   new (p_v_dciso1) ViewFloat4D("pointer_view_dciso1", MAX_BLOCKS_CLINIC, JMT, KM, IMT); 
//   new (p_v_dciso2) ViewFloat4D("pointer_view_dciso2", MAX_BLOCKS_CLINIC, JMT, KM, IMT); 
// #endif // isopycmixspatialvar

//   new (p_v_kisrpl) ViewDouble1D("pointer_view_kisrpl", KM); 

//   new (p_v_krplin) ViewInt1D("pointer_view_krplin", NRPL); 

//   new (p_v_zt)   ViewDouble1D("pointer_view_zt",   KM); 
//   new (p_v_dzw)  ViewDouble1D("pointer_view_dzw",  KM+1); 
//   new (p_v_dzwr) ViewDouble1D("pointer_view_dzwr", KM+1); 
//   new (p_v_dzr)  ViewDouble1D("pointer_view_dzr",  KM); 

//   new (p_v_tmask) ViewDouble4D("pointer_view_tmask", MAX_BLOCKS_CLINIC, JMT, KM, IMT); 
//   new (p_v_f3) ViewDouble3D("pointer_view_f3", MAX_BLOCKS_CLINIC, JMT, IMT); 

//   ViewDouble1D::HostMirror h_v_dptlim = create_mirror_view(*p_v_dptlim);

//   ViewDouble1D::HostMirror h_v_fzisop = create_mirror_view(*p_v_fzisop);
//   ViewDouble3D::HostMirror h_v_ahisop = create_mirror_view(*p_v_ahisop);
//   ViewDouble3D::HostMirror h_v_athkdf = create_mirror_view(*p_v_athkdf);

//   ViewDouble5D::HostMirror h_v_e = create_mirror_view(*p_v_e);

//   ViewDouble5D::HostMirror h_v_rhoi = create_mirror_view(*p_v_rhoi);

//   ViewDouble5D::HostMirror h_v_k1 = create_mirror_view(*p_v_k1);
//   ViewDouble5D::HostMirror h_v_k2 = create_mirror_view(*p_v_k2);
//   ViewDouble5D::HostMirror h_v_k3 = create_mirror_view(*p_v_k3);

//   ViewDouble4D::HostMirror h_v_adv_vetiso = create_mirror_view(*p_v_adv_vetiso);
//   ViewDouble4D::HostMirror h_v_adv_vbtiso = create_mirror_view(*p_v_adv_vbtiso);
//   ViewDouble4D::HostMirror h_v_adv_vntiso = create_mirror_view(*p_v_adv_vntiso);

// #ifdef isopycmixspatialvar
//   ViewFloat4D::HostMirror h_v_dciso1 = create_mirror_view(*p_v_dciso1);
//   ViewFloat4D::HostMirror h_v_dciso2 = create_mirror_view(*p_v_dciso2);
// #endif // isopycmixspatialvar

//   ViewDouble1D::HostMirror h_v_kisrpl = create_mirror_view(*p_v_kisrpl);

//   ViewInt1D::HostMirror h_v_krplin = create_mirror_view(*p_v_krplin);

//   ViewDouble1D::HostMirror h_v_zt = create_mirror_view(*p_v_zt);
//   ViewDouble1D::HostMirror h_v_dzw = create_mirror_view(*p_v_dzw);
//   ViewDouble1D::HostMirror h_v_dzwr = create_mirror_view(*p_v_dzwr);
//   ViewDouble1D::HostMirror h_v_dzr = create_mirror_view(*p_v_dzr);

//   ViewDouble4D::HostMirror h_v_tmask = create_mirror_view(*p_v_tmask);
//   ViewDouble3D::HostMirror h_v_f3 = create_mirror_view(*p_v_f3);

//   for (int n = 0; n < CppIsopycMod::NRPL + 1; ++n) {
//     h_v_dptlim(n) = dptlim[n];
//   }
//   for (int k = 0; k < KM; ++k) {
//     h_v_fzisop(k) = fzisop[k];
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int i = 0; i < IMT; ++i) {
//         h_v_ahisop(iblock, j, i) = ahisop[iblock][j][i];
//         h_v_athkdf(iblock, j, i) = athkdf[iblock][j][i];
//       }
//     }
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int m = 0; m < 3; ++m) {
//       for (int j = 0; j < JMT; ++j) {
//         for (int k = 0; k < KMP1; ++k) {
//           for (int i = 0; i < IMT; ++i) {
//             //h_v_e(i, k, j, m, iblock) = e[iblock][m][j][k][i];
//           }
//         }
//       }
//     }
//   }

//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int m = 0; m < NRPL; ++m) {
//       for (int j = 0; j < JMT; ++j) {
//         for (int k = 0; k < KM+1; ++k) {
//           for (int i = 0; i < IMT; ++i) {
//             //h_v_rhoi(i, k, j, m, iblock) = rhoi[iblock][m][j][k][i];
//           }
//         }
//       }
//     }
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int k = 0; k < KM+1; ++k) {
//         for (int i = 0; i < IMT; ++i) {
//           //h_v_k1(iblock, 0, j, k, i) = k1[iblock][0][j][k][i];
//         }
//       }
//     }
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int m = 0; m < 3; ++m) {
//       for (int j = 0; j < JMT; ++j) {
//         for (int k = 0; k < KM+1; ++k) {
//           for (int i = 0; i < IMT; ++i) {
//             //h_v_k3(i, k, j, m, iblock) = k3[iblock][m][j][k][i];
//           }
//         }
//       }
//     }
//   }

//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int k = 0; k < KM+1; ++k) {
//         for (int i = 0; i < IMT; ++i) {
// //          h_v_adv_vbtiso(iblock, j, k, i) = adv_vbtiso[iblock][j][k][i];
//         }
//       }
//     }
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int k = 0; k < KM; ++k) {
//         for (int i = 0; i < IMT; ++i) {
// //         h_v_adv_vetiso(iblock, j, k, i) = adv_vetiso[iblock][j][k][i];
// //          h_v_adv_vntiso(iblock, j, k, i) = adv_vntiso[iblock][j][k][i];
//         }
//       }
//     }
//   }

// #ifdef isopycmixspatialvar
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int k = 0; k < KM; ++k) {
//         for (int i = 0; i < IMT; ++i) {
//           h_v_dciso1(iblock, j, k, i) = adv_dciso1[iblock][j][k][i];
//           h_v_dciso2(iblock, j, k, i) = adv_dciso2[iblock][j][k][i];
//         }
//       }
//     }
//   }
// #endif // isopycmixspatialvar
//   for (int k = 0; k < KM; ++k) {
//     h_v_kisrpl(k) = kisrpl[k];
//   }
//   for (int n = 0; n < CppIsopycMod::NRPL; ++n) {
//     h_v_krplin(n) = krplin[n];
//   }
//   for (int k = 0; k < KM; ++k) {
//     h_v_zt(k) = zt[k];
//   }
//   for (int k = 0; k < KM+1; ++k) {
//     h_v_dzw(k) = dzw[k];
//     h_v_dzwr(k) = dzwr[k];
//   }
//   for (int k = 0; k < KM; ++k) {
//     h_v_dzr(k) = dzr[k];
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int k = 0; k < KM; ++k) {
//         for (int i = 0; i < IMT; ++i) {
//           h_v_tmask(iblock, j, k, i) = tmask[iblock][j][k][i];
//         }
//       }
//     }
//   }
//   for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
//     for (int j = 0; j < JMT; ++j) {
//       for (int i = 0; i < IMT; ++i) {
//         h_v_f3(iblock, j, i) = f3[iblock][j][i];
//       }
//     }
//   }

//   deep_copy(dev, *p_v_dptlim, h_v_dptlim);

//   deep_copy(dev, *p_v_fzisop, h_v_fzisop);

//   deep_copy(dev, *p_v_ahisop, h_v_ahisop);
//   deep_copy(dev, *p_v_athkdf, h_v_athkdf);

//   deep_copy(dev, *p_v_e, h_v_e);

//   deep_copy(dev, *p_v_rhoi, h_v_rhoi);

//   deep_copy(dev, *p_v_k1, h_v_k1);
//   deep_copy(dev, *p_v_k2, h_v_k2);
//   deep_copy(dev, *p_v_k3, h_v_k3);

//   deep_copy(dev, *p_v_adv_vetiso, h_v_adv_vetiso);
//   deep_copy(dev, *p_v_adv_vbtiso, h_v_adv_vbtiso);
//   deep_copy(dev, *p_v_adv_vntiso, h_v_adv_vntiso);

// #ifdef isopycmixspatialvar
//   deep_copy(dev, *p_v_dciso1, h_v_dciso1);
//   deep_copy(dev, *p_v_dciso2, h_v_dciso2);
// #endif // isopycmixspatialvar

//   deep_copy(dev, *p_v_kisrpl, h_v_kisrpl);

//   deep_copy(dev, *p_v_krplin, h_v_krplin);

//   deep_copy(dev, *p_v_zt, h_v_zt);

//   deep_copy(dev, *p_v_dzw, h_v_dzw);
//   deep_copy(dev, *p_v_dzwr, h_v_dzwr);

//   deep_copy(dev, *p_v_dzr, h_v_dzr);

//   deep_copy(dev, *p_v_tmask, h_v_tmask);
//   deep_copy(dev, *p_v_f3, h_v_f3);
//   return ;
// }
// #endif // ISO
static void kokkos_init_output_mod() {
  using namespace CppOutputMod;
  using namespace KokkosOutputMod;

  auto dev = Kokkos::DefaultExecutionSpace();
#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES

  p_v_icmon = (ViewFloat4D *) malloc(sizeof(ViewFloat4D)
      * IMT * JMT * 2 * MAX_BLOCKS_CLINIC);

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
#endif // ISO_TYPE_BF

#ifdef ISO
#ifdef ISOOUT
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
#endif // SMAG_OUT

#endif // LOWRES

#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

#ifdef LOWRES

  new (p_v_icmon) ViewFloat4D("pointer_view_icmon",
      MAX_BLOCKS_CLINIC, 2, JMT, IMT);

  UnManagedViewFloat4D h_v_icmon (&icmon[0][0][0][0],
      MAX_BLOCKS_CLINIC, 2, JMT, IMT);

  deep_copy(dev, *p_v_icmon, h_v_icmon);

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
#endif // ISO_TYPE_BF

#ifdef ISO
#ifdef ISOOUT
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
#endif // SMAG_OUT

#endif // LOWRES

#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC
#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE

#ifdef LOWRES

  new (p_v_icmon) ViewFloat4D(&icmon[0][0][0][0],
      MAX_BLOCKS_CLINIC, 2, JMT, IMT);

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
#endif // ISO_TYPE_BF

#ifdef ISO
#ifdef ISOOUT
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
#endif // SMAG_OUT

#endif // LOWRES

#ifdef DAILYACC

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  return ;
}

static void kokkos_init_pconst_mod() {

  using namespace CppPconstMod;
  using namespace KokkosPconstMod;

  auto dev = Kokkos::DefaultExecutionSpace();
  // p_v_j_global = (ViewInt1D *) malloc(sizeof(ViewInt1D) * JMT);
  // p_v_i_global = (ViewInt1D *) malloc(sizeof(ViewInt1D) * IMT);

  p_v_vit = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_viv = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // p_v_ahv_back = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * JMT_GLOBAL);

  // p_v_na = (ViewInt3D *) malloc(sizeof(ViewInt3D) *
      // MAX_BLOCKS_CLINIC * JMT * IMT);

#if (defined NETCDF) || (defined ALL)
  p_v_lon = (ViewFloat1D *) malloc(sizeof(ViewFloat1D));
  p_v_lat = (ViewFloat1D *) malloc(sizeof(ViewFloat1D));

  p_v_lon_o = (ViewFloat2D *) malloc(sizeof(ViewFloat2D));
  p_v_lat_o = (ViewFloat2D *) malloc(sizeof(ViewFloat2D));

  // p_v_ulon_o = (ViewFloat2D *) malloc(sizeof(ViewFloat2D) * 
  //     IMT_GLOBAL * JMT_GLOBAL);
  // p_v_ulat_o = (ViewFloat2D *) malloc(sizeof(ViewFloat2D) * 
  //     IMT_GLOBAL * JMT_GLOBAL);

  // p_v_lev  = (ViewFloat1D *) malloc(sizeof(ViewFloat1D) *  KM);
  // p_v_lev1 = (ViewFloat1D *) malloc(sizeof(ViewFloat1D) * (KM + 1));
#endif // NETCDF || ALL

  p_v_s_lon = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_s_lat = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));


  p_v_zkt  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_dzp  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_odzp = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_odzt = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_odz_pt = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_zkp  = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_ebea = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ebeb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  // p_v_ebla = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);
  // p_v_eblb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) *
  //     MAX_BLOCKS_CLINIC * JMT * IMT);

  p_v_epea = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_epeb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_epla = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_eplb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_rrd1 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_rrd2 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_ohbt = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_ohbu = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_dzph = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_hbx = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_hby = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_snlat = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_to = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
  p_v_so = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_c = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_po = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_akmu = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_akmt = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_akt = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  // p_v_am = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * JMT);
  // p_v_ah = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * JMT);

  // p_v_am3 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KM * JMT * IMT);
  // p_v_ah3 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KM * JMT * IMT);

  // p_v_amx = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KM * JMT * IMT);
  // p_v_amy = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KM * JMT * IMT);

  // p_v_nmonth  = (ViewInt1D *) malloc(sizeof(ViewInt1D) * 12);
  // p_v_nnmonth = (ViewInt1D *) malloc(sizeof(ViewInt1D) * 12);

#ifdef TIDEMIX
  p_v_fz_tide = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_ak_tide = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_fztidal    = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_richardson = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp3_tidal  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  // p_v_ak_tide1   = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KM * JMT * IMT);
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
  p_v_wp1_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp2_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp3_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp4_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp5_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp6_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp7_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp8_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp10_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp11_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp12_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp13_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_wk1_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wk2_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wk3_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wk4_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_fcor_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_fcort_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_alpha_canuto = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_beta_canuto  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
#endif // CANUTOMIXOUT

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_vit) ViewDouble4D("pointer_view_vit",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_viv) ViewDouble4D("pointer_view_viv",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

#if (defined NETCDF) || (defined ALL)
  new (p_v_lon) ViewFloat1D("pointer_view_lon", IMT_GLOBAL);
  new (p_v_lat) ViewFloat1D("pointer_view_lat", JMT_GLOBAL);

  new (p_v_lon_o) ViewFloat2D("pointer_view_lon_o", JMT_GLOBAL, IMT_GLOBAL);
  new (p_v_lat_o) ViewFloat2D("pointer_view_lat_o", JMT_GLOBAL, IMT_GLOBAL);

#endif // NETCDF || ALL
  new (p_v_s_lon) ViewDouble2D("pointer_view_s_lon", S_JMT, S_IMT);
  new (p_v_s_lat) ViewDouble2D("pointer_view_s_lat", S_JMT, S_IMT);


  new (p_v_zkt)  ViewDouble1D("pointer_view_zkt",  KM);
  new (p_v_dzp)  ViewDouble1D("pointer_view_dzp",  KM);
  new (p_v_odzp) ViewDouble1D("pointer_view_odzp", KM);
  new (p_v_odzt) ViewDouble1D("pointer_view_odzt", KM);

  new (p_v_odz_pt) ViewDouble2D("pointer_view_odz_pt", KM, 2);

  new (p_v_zkp)  ViewDouble1D("pointer_view_zkp",  KMP1);

  new (p_v_ebea) ViewDouble3D("pointer_view_ebea",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ebeb) ViewDouble3D("pointer_view_ebeb",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_epea) ViewDouble3D("pointer_view_epea",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_epeb) ViewDouble3D("pointer_view_epeb",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_epla) ViewDouble3D("pointer_view_epla",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_eplb) ViewDouble3D("pointer_view_eplb",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_rrd1) ViewDouble3D("pointer_view_rrd1",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_rrd2) ViewDouble3D("pointer_view_rrd2",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_ohbt) ViewDouble3D("pointer_view_ohbt",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ohbu) ViewDouble3D("pointer_view_ohbu",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_dzph) ViewDouble3D("pointer_view_dzph",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_hbx) ViewDouble3D("pointer_view_hbx",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_hby) ViewDouble3D("pointer_view_hby",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_snlat) ViewDouble3D("pointer_view_snlat",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_to) ViewDouble1D("pointer_view_to", KM);
  new (p_v_so) ViewDouble1D("pointer_view_so", KM);

  new (p_v_c) ViewDouble2D("pointer_view_c", 9, KM);

  new (p_v_po) ViewDouble1D("pointer_view_po", KM);

  new (p_v_akmu) ViewDouble4D("pointer_view_akmu",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_akmt) ViewDouble4D("pointer_view_akmt",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_akt) ViewDouble5D("pointer_view_akt",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef TIDEMIX
  new (p_v_fz_tide) ViewDouble4D("pointer_view_fz_tide",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ak_tide) ViewDouble4D("pointer_view_ak_tide",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_fztidal)    ViewDouble4D("pointer_view_fztidal",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_richardson) ViewDouble4D("pointer_view_richardson",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp3_tidal)  ViewDouble4D("pointer_view_wp3_tidal",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
  new (p_v_wp1_canuto)  ViewDouble4D("pointer_view_wp1_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp2_canuto)  ViewDouble4D("pointer_view_wp2_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp3_canuto)  ViewDouble4D("pointer_view_wp3_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp4_canuto)  ViewDouble4D("pointer_view_wp4_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp5_canuto)  ViewDouble4D("pointer_view_wp5_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp6_canuto)  ViewDouble4D("pointer_view_wp6_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp7_canuto)  ViewDouble4D("pointer_view_wp7_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp8_canuto)  ViewDouble4D("pointer_view_wp8_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp10_canuto) ViewDouble4D("pointer_view_wp10_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp11_canuto) ViewDouble4D("pointer_view_wp11_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp12_canuto) ViewDouble4D("pointer_view_wp12_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp13_canuto) ViewDouble4D("pointer_view_wp13_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_wk1_canuto) ViewDouble4D("pointer_view_wk1_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk2_canuto) ViewDouble4D("pointer_view_wk2_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk3_canuto) ViewDouble4D("pointer_view_wk3_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk4_canuto) ViewDouble4D("pointer_view_wk4_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_fcor_canuto)  ViewDouble4D("pointer_view_fcor_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_fcort_canuto) ViewDouble4D("pointer_view_fcort_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_alpha_canuto) ViewDouble4D("pointer_view_alpha_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_beta_canuto)  ViewDouble4D("pointer_view_beta_canuto",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // CANUTOMIXOUT

  UnManagedViewDouble4D h_v_vit (&vit[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_viv (&viv[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

#if (defined NETCDF) || (defined ALL)
  UnManagedViewFloat1D h_v_lon (&lon[0], IMT_GLOBAL);
  UnManagedViewFloat1D h_v_lat (&lat[0], JMT_GLOBAL);

  UnManagedViewFloat2D h_v_lon_o (&lon_o[0][0], JMT_GLOBAL, IMT_GLOBAL);
  UnManagedViewFloat2D h_v_lat_o (&lat_o[0][0], JMT_GLOBAL, IMT_GLOBAL);

#endif // NETCDF || ALL

  UnManagedViewDouble1D h_v_zkt  (&zkt[0], KM);
  UnManagedViewDouble1D h_v_dzp  (&dzp[0], KM);
  UnManagedViewDouble1D h_v_odzp (&odzp[0], KM);
  UnManagedViewDouble1D h_v_odzt (&odzt[0], KM);

  UnManagedViewDouble2D::HostMirror h_v_odz_pt = create_mirror_view(*p_v_odz_pt);

  UnManagedViewDouble1D h_v_zkp (&zkp[0], KMP1);

  UnManagedViewDouble3D h_v_ebea (&ebea[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_ebeb (&ebeb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_epea (&epea[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_epeb (&epeb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_epla (&epla[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_eplb (&eplb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_rrd1 (&rrd1[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_rrd2 (&rrd2[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_ohbt (&ohbt[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_ohbu (&ohbu[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_dzph (&dzph[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_hbx (&hbx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_hby (&hby[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_snlat (&snlat[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble1D h_v_to (&to[0], KM);
  UnManagedViewDouble1D h_v_so (&so[0], KM);

  UnManagedViewDouble2D h_v_c (&c[0][0], 9, KM);

  UnManagedViewDouble1D h_v_po (&po[0], KM);

  UnManagedViewDouble4D h_v_akmu (&akmu[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_akmt (&akmt[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_akt (&akt[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef TIDEMIX
  UnManagedViewDouble4D h_v_fz_tide (&fz_tide[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_ak_tide (&ak_tide[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble4D h_v_fztidal (&fztidal[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_richardson (&richardson[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp3_tidal (&wp3_tidal[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
  UnManagedViewDouble4D h_v_wp1_canuto (&wp1_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp2_canuto (&wp2_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp3_canuto (&wp3_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp4_canuto (&wp4_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp5_canuto (&wp5_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp6_canuto (&wp6_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp7_canuto (&wp7_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp8_canuto (&wp8_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp10_canuto (&wp10_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp11_canuto (&wp11_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp12_canuto (&wp12_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wp13_canuto (&wp13_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble4D h_v_wk1_canuto (&wk1_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wk2_canuto (&wk2_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wk3_canuto (&wk3_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_wk4_canuto (&wk4_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble4D h_v_fcor_canuto (&fcor_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_fcort_canuto (&fcort_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_alpha_canuto (&alpha_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble4D h_v_beta_canuto (&beta_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // CANUTOMIXOUT

  for (int k = 0; k < KM; ++k) {
    h_v_odz_pt(k, 0) = odzp[k];
    h_v_odz_pt(k, 1) = odzt[k];
  }

  deep_copy(dev, *p_v_vit,          h_v_vit);
  deep_copy(dev, *p_v_viv,          h_v_viv);
#if (defined NETCDF) || (defined ALL)
  deep_copy(dev, *p_v_lon,          h_v_lon);
  deep_copy(dev, *p_v_lat,          h_v_lat);
  deep_copy(dev, *p_v_lon_o,        h_v_lon_o);
  deep_copy(dev, *p_v_lat_o,        h_v_lat_o);
#endif // NETCDF || ALL
  deep_copy(dev, *p_v_zkt,          h_v_zkt);
  deep_copy(dev, *p_v_dzp,          h_v_dzp);
  deep_copy(dev, *p_v_odzp,         h_v_odzp);
  deep_copy(dev, *p_v_odzt,         h_v_odzt);
  deep_copy(dev, *p_v_odz_pt,       h_v_odz_pt);
  deep_copy(dev, *p_v_zkp,          h_v_zkp);
  deep_copy(dev, *p_v_ebea,         h_v_ebea);
  deep_copy(dev, *p_v_ebeb,         h_v_ebeb);
  deep_copy(dev, *p_v_epea,         h_v_epea);
  deep_copy(dev, *p_v_epeb,         h_v_epeb);
  deep_copy(dev, *p_v_epla,         h_v_epla);
  deep_copy(dev, *p_v_eplb,         h_v_eplb);
  deep_copy(dev, *p_v_rrd1,         h_v_rrd1);
  deep_copy(dev, *p_v_rrd2,         h_v_rrd2);
  deep_copy(dev, *p_v_ohbt,         h_v_ohbt);
  deep_copy(dev, *p_v_ohbu,         h_v_ohbu);
  deep_copy(dev, *p_v_dzph,         h_v_dzph);
  deep_copy(dev, *p_v_hbx,          h_v_hbx);
  deep_copy(dev, *p_v_hby,          h_v_hby);
  deep_copy(dev, *p_v_snlat,        h_v_snlat);
  deep_copy(dev, *p_v_to,           h_v_to);
  deep_copy(dev, *p_v_so,           h_v_so);
  deep_copy(dev, *p_v_c,            h_v_c);
  deep_copy(dev, *p_v_po,           h_v_po);
  deep_copy(dev, *p_v_akmu,         h_v_akmu);
  deep_copy(dev, *p_v_akmt,         h_v_akmt);
  deep_copy(dev, *p_v_akt,          h_v_akt);
#ifdef TIDEMIX
  deep_copy(dev, *p_v_fz_tide,      h_v_fz_tide);
  deep_copy(dev, *p_v_ak_tide,      h_v_ak_tide);
  deep_copy(dev, *p_v_fztidal,      h_v_fztidal);
  deep_copy(dev, *p_v_richardson,   h_v_richardson);
  deep_copy(dev, *p_v_wp3_tidal,    h_v_wp3_tidal);
#endif // TIDEMIX
#ifdef CANUTOMIXOUT
  deep_copy(dev, *p_v_wp1_canuto,   h_v_wp1_canuto);
  deep_copy(dev, *p_v_wp2_canuto,   h_v_wp2_canuto);
  deep_copy(dev, *p_v_wp3_canuto,   h_v_wp3_canuto);
  deep_copy(dev, *p_v_wp4_canuto,   h_v_wp4_canuto);
  deep_copy(dev, *p_v_wp5_canuto,   h_v_wp5_canuto);
  deep_copy(dev, *p_v_wp6_canuto,   h_v_wp6_canuto);
  deep_copy(dev, *p_v_wp7_canuto,   h_v_wp7_canuto);
  deep_copy(dev, *p_v_wp8_canuto,   h_v_wp8_canuto);
  deep_copy(dev, *p_v_wp10_canuto,  h_v_wp10_canuto);
  deep_copy(dev, *p_v_wp11_canuto,  h_v_wp11_canuto);
  deep_copy(dev, *p_v_wp12_canuto,  h_v_wp12_canuto);
  deep_copy(dev, *p_v_wp13_canuto,  h_v_wp13_canuto);

  deep_copy(dev, *p_v_wk1_canuto,   h_v_wk1_canuto);
  deep_copy(dev, *p_v_wk2_canuto,   h_v_wk2_canuto);
  deep_copy(dev, *p_v_wk3_canuto,   h_v_wk3_canuto);
  deep_copy(dev, *p_v_wk4_canuto,   h_v_wk4_canuto);

  deep_copy(dev, *p_v_fcor_canuto,  h_v_fcor_canuto);
  deep_copy(dev, *p_v_fcort_canuto, h_v_fcort_canuto);

  deep_copy(dev, *p_v_alpha_canuto, h_v_alpha_canuto);
  deep_copy(dev, *p_v_beta_canuto,  h_v_beta_canuto);
#endif // CANUTOMIXOUT

#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_vit) ViewDouble4D(&vit[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_viv) ViewDouble4D(&viv[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

#if (defined NETCDF) || (defined ALL)
  new (p_v_lon) ViewFloat1D(&lon[0], IMT_GLOBAL);
  new (p_v_lat) ViewFloat1D(&lat[0], JMT_GLOBAL);

  new (p_v_lon_o) ViewFloat2D(&lon_o[0][0], JMT_GLOBAL, IMT_GLOBAL);
  new (p_v_lat_o) ViewFloat2D(&lat_o[0][0], JMT_GLOBAL, IMT_GLOBAL);

#endif // NETCDF || ALL
	s_lon  = new double[S_IMT * S_JMT];
	s_lat  = new double[S_IMT * S_JMT];
  new (p_v_s_lon) ViewDouble2D(s_lon, S_JMT, S_IMT);
  new (p_v_s_lat) ViewDouble2D(s_lat, S_JMT, S_IMT);


  new (p_v_zkt)  ViewDouble1D(&zkt[0],  KM);
  new (p_v_dzp)  ViewDouble1D(&dzp[0],  KM);
  new (p_v_odzp) ViewDouble1D(&odzp[0], KM);
  new (p_v_odzt) ViewDouble1D(&odzt[0], KM);

  odz_pt = new double [2 * KM];
  new (p_v_odz_pt) ViewDouble2D(odz_pt, KM, 2);

  new (p_v_zkp)  ViewDouble1D(&zkp[0],  KMP1);

  new (p_v_ebea) ViewDouble3D(&ebea[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ebeb) ViewDouble3D(&ebeb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_epea) ViewDouble3D(&epea[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_epeb) ViewDouble3D(&epeb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_epla) ViewDouble3D(&epla[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_eplb) ViewDouble3D(&eplb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_rrd1) ViewDouble3D(&rrd1[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_rrd2) ViewDouble3D(&rrd2[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_ohbt) ViewDouble3D(&ohbt[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_ohbu) ViewDouble3D(&ohbu[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_dzph) ViewDouble3D(&dzph[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_hbx) ViewDouble3D(&hbx[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_hby) ViewDouble3D(&hby[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_snlat) ViewDouble3D(&snlat[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_to) ViewDouble1D(&to[0], KM);
  new (p_v_so) ViewDouble1D(&so[0], KM);

  new (p_v_c) ViewDouble2D(&c[0][0], 9, KM);

  new (p_v_po) ViewDouble1D(&po[0], KM);

  new (p_v_akmu) ViewDouble4D(&akmu[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_akmt) ViewDouble4D(&akmt[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_akt) ViewDouble5D(&akt[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef TIDEMIX
  new (p_v_fz_tide) ViewDouble4D(&fz_tide[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ak_tide) ViewDouble4D(&ak_tide[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_fztidal)    ViewDouble4D(&fztidal[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_richardson) ViewDouble4D(&richardson[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp3_tidal)  ViewDouble4D(&wp3_tidal[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // TIDEMIX

#ifdef CANUTOMIXOUT
  new (p_v_wp1_canuto)  ViewDouble4D(&wp1_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp2_canuto)  ViewDouble4D(&wp2_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp3_canuto)  ViewDouble4D(&wp3_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp4_canuto)  ViewDouble4D(&wp4_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp5_canuto)  ViewDouble4D(&wp5_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp6_canuto)  ViewDouble4D(&wp6_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp7_canuto)  ViewDouble4D(&wp7_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp8_canuto)  ViewDouble4D(&wp8_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp10_canuto) ViewDouble4D(&wp10_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp11_canuto) ViewDouble4D(&wp11_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp12_canuto) ViewDouble4D(&wp12_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp13_canuto) ViewDouble4D(&wp13_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_wk1_canuto) ViewDouble4D(&wk1_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk2_canuto) ViewDouble4D(&wk2_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk3_canuto) ViewDouble4D(&wk3_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wk4_canuto) ViewDouble4D(&wk4_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_fcor_canuto)  ViewDouble4D(&fcor_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_fcort_canuto) ViewDouble4D(&fcort_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_alpha_canuto) ViewDouble4D(&alpha_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_beta_canuto)  ViewDouble4D(&beta_canuto[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // CANUTOMIXOUT

  for (int k = 0; k < KM; ++k) {
    (*p_v_odz_pt)(k, 0) = odzp[k];
    (*p_v_odz_pt)(k, 1) = odzt[k];
  }

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  return ;
}

static void kokkos_init_pmix_mod() {

  using namespace CppPmixMod;
  using namespace KokkosPmixMod;

  auto dev = Kokkos::DefaultExecutionSpace();
  // p_v_ric  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KMM1 * JMT * IMT);
  p_v_rict = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // p_v_rict_replace = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KMM1 * JMT * IMT);

  // p_v_rict_ref = (ViewDouble3D *) malloc(sizeof(ViewDouble3D) * 
  //     MAX_BLOCKS_CLINIC * JMT * IMT);

  p_v_rit = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // p_v_riu = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * (KM+1) * JMT * IMT);

  p_v_ricdt    = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_ricdttms = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // p_v_ridt = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KMM1 * JMT * IMT);

  // p_v_s2u = (ViewDouble4D *) malloc(sizeof(ViewDouble4D) *
  //     MAX_BLOCKS_CLINIC * KMM1 * JMT * IMT);
  p_v_s2t = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

#ifdef SOLAR
  p_v_pen = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));
#endif // SOLAR

#ifdef SOLARCHLORO
  p_v_pen_chl = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
#endif // SOLARCHLORO

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  new (p_v_rict) ViewDouble4D("pointer_view_rict",
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

  new (p_v_rit) ViewDouble4D("pointer_view_rit",
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

  new (p_v_ricdt) ViewDouble4D("pointer_view_ricdt",
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);
  new (p_v_ricdttms) ViewDouble4D("pointer_view_ricdttms",
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);
  new (p_v_s2t) ViewDouble4D("pointer_view_s2t",
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

#ifdef SOLAR
  new (p_v_pen) ViewDouble1D("pointer_view_pen", KMM1);
#endif // SOLAR

#ifdef SOLARCHLORO
  new (p_v_pen_chl) ViewDouble4D("pointer_view_pen_chl",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // SOLARCHLORO

  UnManagedViewDouble4D h_v_s2t (&s2t[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

#ifdef SOLAR
  UnManagedViewDouble1D h_v_pen (&pen[0], KMM1);
#endif // SOLAR

#ifdef SOLARCHLORO
  UnManagedViewDouble4D h_v_pen_chl (&pen_chl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // SOLARCHLORO

  deep_copy(dev, *p_v_s2t, h_v_s2t);
#ifdef SOLAR
  deep_copy(dev, *p_v_pen, h_v_pen);
#endif // SOLAR

#ifdef SOLARCHLORO
  deep_copy(dev, *p_v_pen_chl, h_v_pen_chl);
#endif // SOLARCHLORO
#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_rict) ViewDouble4D(&rict[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

  new (p_v_rit) ViewDouble4D(&rit[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

  new (p_v_ricdt) ViewDouble4D(&ricdt[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);
  new (p_v_ricdttms) ViewDouble4D(&ricdttms[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);
  new (p_v_s2t) ViewDouble4D(&s2t[0][0][0][0],
      MAX_BLOCKS_CLINIC, KMM1, JMT, IMT);

#ifdef SOLAR
  new (p_v_pen) ViewDouble1D (&pen[0], KMM1);
#endif // SOLAR

#ifdef SOLARCHLORO
  new (p_v_pen_chl) ViewDouble4D (&pen_chl[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // SOLARCHLORO

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
  return ;
}

static void kokkos_init_work_mod() {

  using namespace CppWorkMod;
  using namespace KokkosWorkMod;

  p_v_pxb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_pyb = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_pax = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_pay = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_whx = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_why = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_wgp = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_wka = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_work_1 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_work_2 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_work_3 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_temp = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_uk = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_vk = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_work = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  // p_v_wkk = (ViewDouble1D *) malloc(sizeof(ViewDouble1D) * KMP1);

  p_v_wkb = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wkc = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wkd = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_tf = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_stf = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  // p_v_buffer_real4 = (ViewFloat2D *) malloc(sizeof(ViewFloat2D) *
  //     IMT_GLOBAL * JMT_GLOBAL);

  // TODO can opt
  p_v_work1_g = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_work2_g = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  // p_v_work3_g = (ViewDouble2D *) malloc(sizeof(ViewDouble2D) *
  //     IMT_GLOBAL * JMT_GLOBAL);

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  new (p_v_pxb) ViewDouble3D("pointer_view_pxb",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_pyb) ViewDouble3D("pointer_view_pyb",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_pax) ViewDouble3D("pointer_view_pax",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_pay) ViewDouble3D("pointer_view_pay", 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_whx) ViewDouble3D("pointer_view_whx", 
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_why) ViewDouble3D("pointer_view_why", 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_wgp) ViewDouble3D("pointer_view_wgp", 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_wka) ViewDouble4D("pointer_view_wka",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_work_1) ViewDouble4D("pointer_view_work_1",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_work_2) ViewDouble4D("pointer_view_work_2",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_work_3) ViewDouble4D("pointer_view_work_2",
      MAX_BLOCKS_CLINIC, KM+1, JMT, IMT);

  new (p_v_temp) ViewDouble4D("pointer_view_temp",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_uk) ViewDouble4D("pointer_view_uk",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_vk) ViewDouble4D("pointer_view_vk",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_work) ViewDouble3D("pointer_view_work",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_wkb) ViewDouble4D("pointer_view_wkb",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wkc) ViewDouble4D("pointer_view_wkc",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wkd) ViewDouble4D("pointer_view_wkd",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_tf) ViewDouble4D("pointer_view_tf",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_stf) ViewDouble3D("pointer_view_stf",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_work1_g) ViewDouble2D("pointer_view_work1_g",
      JMT_GLOBAL, IMT_GLOBAL); 
  new (p_v_work2_g) ViewDouble2D("pointer_view_work2_g",
      JMT_GLOBAL, IMT_GLOBAL); 

  UnManagedViewDouble3D h_v_pxb (&pxb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_pyb (&pyb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_pax (&pax[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_pay (&pay[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_whx (&whx[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);
  UnManagedViewDouble3D h_v_why (&why[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble3D h_v_wgp (&wgp[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble4D h_v_wka (&wka[0][0][0][0], 
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble3D h_v_work (&work[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  auto dev = Kokkos::DefaultExecutionSpace();
  deep_copy(dev, *p_v_pxb,  h_v_pxb);
  deep_copy(dev, *p_v_pyb,  h_v_pyb);
  deep_copy(dev, *p_v_pax,  h_v_pax);
  deep_copy(dev, *p_v_pay,  h_v_pay);
  deep_copy(dev, *p_v_whx,  h_v_whx);
  deep_copy(dev, *p_v_why,  h_v_why);
  deep_copy(dev, *p_v_wgp,  h_v_wgp);
  deep_copy(dev, *p_v_wka,  h_v_wka);
  deep_copy(dev, *p_v_work, h_v_work);
#else // KOKKOS_ENABLE_DEVICE_MEM_SPACE
  new (p_v_pxb) ViewDouble3D(&pxb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_pyb) ViewDouble3D(&pyb[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_pax) ViewDouble3D(&pax[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_pay) ViewDouble3D(&pay[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_whx) ViewDouble3D(&whx[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_why) ViewDouble3D(&why[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_wgp) ViewDouble3D(&wgp[0][0][0], 
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_wka) ViewDouble4D(&wka[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);


  work_1 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  work_2 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  work_3 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_work_1) ViewDouble4D(work_1,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_work_2) ViewDouble4D(work_2,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_work_3) ViewDouble4D(work_3,
      MAX_BLOCKS_CLINIC, KM+1, JMT, IMT);

  temp = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];

  new (p_v_temp) ViewDouble4D(temp,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  uk = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  vk = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];

  new (p_v_uk) ViewDouble4D(uk,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_vk) ViewDouble4D(vk,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_work) ViewDouble3D(&work[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  wkb = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  wkc = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  wkd = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_wkb) ViewDouble4D(wkb,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wkc) ViewDouble4D(wkc,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wkd) ViewDouble4D(wkd,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  tf  = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  stf = new double[MAX_BLOCKS_CLINIC * JMT * IMT];

  new (p_v_tf) ViewDouble4D(tf,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_stf) ViewDouble3D(stf,
      MAX_BLOCKS_CLINIC, JMT, IMT);

  work1_g = new double[JMT_GLOBAL * IMT_GLOBAL];
  work2_g = new double[JMT_GLOBAL * IMT_GLOBAL];
  new (p_v_work1_g) ViewDouble2D(work1_g,
      JMT_GLOBAL, IMT_GLOBAL); 
  new (p_v_work2_g) ViewDouble2D(work2_g,
      JMT_GLOBAL, IMT_GLOBAL); 
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
  free_work_mod_();
  return ;
}

static void kokkos_init_tracer_mod() {

  using namespace CppTracerMod;
  using namespace KokkosTracerMod;

  p_v_atb = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_net = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_at = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
//extern double (*(&restore_at))[NTRA][JMT][IMT];

  p_v_pdensity = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_amld     = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_tend = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_ax = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_ay = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_az = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_dx = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_dy = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_dz = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_penetrate = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_dt_diff    = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  // p_v_ddy        = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) *
  //     IMT * JMT * KM * NTRA * MAX_BLOCKS_CLINIC);

  p_v_dt_conv    = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  // p_v_dt_restore = (ViewDouble5D *) malloc(sizeof(ViewDouble5D) *
  //     IMT * JMT * KM * NTRA * MAX_BLOCKS_CLINIC);

#ifdef ISO
  p_v_aay_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_ddy_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_ax_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_ay_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_az_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));

  p_v_dx_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_dy_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
  p_v_dz_iso = (ViewDouble5D *) malloc(sizeof(ViewDouble5D));
#endif // ISO

  p_v_licomqice = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_atb) ViewDouble5D("pointer_view_atb",
      MAX_BLOCKS_CLINIC, NTRA, KM+1, JMT, IMT);

  new (p_v_net) ViewDouble4D("pointer_view_net",
      MAX_BLOCKS_CLINIC, NTRA, JMT, IMT);

  new (p_v_at) ViewDouble5D("pointer_view_at",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_pdensity) ViewDouble4D("pointer_view_pdensity",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_amld) ViewDouble3D("pointer_view_amld",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_tend) ViewDouble5D("pointer_view_tend",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ax) ViewDouble5D("pointer_view_ax",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ay) ViewDouble5D("pointer_view_ay",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_az) ViewDouble5D("pointer_view_az",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dx) ViewDouble5D("pointer_view_dx",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dy) ViewDouble5D("pointer_view_dy",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dz) ViewDouble5D("pointer_view_dz",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_penetrate) ViewDouble4D("pointer_view_penetrate",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_dt_diff) ViewDouble5D("pointer_view_dt_diff",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dt_conv) ViewDouble5D("pointer_view_dt_conv",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef ISO
  new (p_v_aay_iso) ViewDouble5D("pointer_view_aay_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ddy_iso) ViewDouble5D("pointer_view_ddy_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_ax_iso) ViewDouble5D("pointer_view_ax_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ay_iso) ViewDouble5D("pointer_view_ay_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_az_iso) ViewDouble5D("pointer_view_az_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dx_iso) ViewDouble5D("pointer_view_dx_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dy_iso) ViewDouble5D("pointer_view_dy_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dz_iso) ViewDouble5D("pointer_view_dz_iso",
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
#endif // ISO

  new (p_v_licomqice) ViewDouble3D("pointer_view_licomqice",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble5D h_v_atb (&atb[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM+1, JMT, IMT);

  UnManagedViewDouble4D h_v_net (&net[0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, JMT, IMT);

  UnManagedViewDouble5D h_v_at (&at[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble4D h_v_pdensity (&pdensity[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  UnManagedViewDouble3D h_v_amld (&amld[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  UnManagedViewDouble5D h_v_tend (&tend[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_ax (&ax[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_ay (&ay[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_az (&az[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_dx (&dx[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_dy (&dy[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_dz (&dz[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble4D h_v_penetrate (&penetrate[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_dt_diff (&dt_diff[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_dt_conv (&dt_conv[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef ISO
  UnManagedViewDouble5D h_v_aay_iso (&aay_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_ddy_iso (&ddy_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_ax_iso (&ax_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_ay_iso (&ay_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_az_iso (&az_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  UnManagedViewDouble5D h_v_dx_iso (&dx_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_dy_iso (&dy_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  UnManagedViewDouble5D h_v_dz_iso (&dz_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
#endif // ISO

  UnManagedViewDouble3D h_v_licomqice (&licomqice[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  auto dev = Kokkos::DefaultExecutionSpace();

  deep_copy(dev, *p_v_atb, h_v_atb);
  deep_copy(dev, *p_v_net, h_v_net);

  deep_copy(dev, *p_v_at,  h_v_at);

  deep_copy(dev, *p_v_pdensity, h_v_pdensity);

  deep_copy(dev, *p_v_amld, h_v_amld);
  deep_copy(dev, *p_v_tend, h_v_tend);

  deep_copy(dev, *p_v_ax, h_v_ax);
  deep_copy(dev, *p_v_ay, h_v_ay);
  deep_copy(dev, *p_v_az, h_v_az);

  deep_copy(dev, *p_v_dx, h_v_dx);
  deep_copy(dev, *p_v_dy, h_v_dy);
  deep_copy(dev, *p_v_dz, h_v_dz);

  deep_copy(dev, *p_v_penetrate,h_v_penetrate);

  deep_copy(dev, *p_v_dt_diff, h_v_dt_diff);

  deep_copy(dev, *p_v_dt_conv, h_v_dt_conv);

#ifdef ISO
  deep_copy(dev, *p_v_aay_iso, h_v_aay_iso);
  deep_copy(dev, *p_v_ddy_iso, h_v_ddy_iso);

  deep_copy(dev, *p_v_ax_iso, h_v_ax_iso);
  deep_copy(dev, *p_v_ay_iso, h_v_ay_iso);
  deep_copy(dev, *p_v_az_iso, h_v_az_iso);
  deep_copy(dev, *p_v_dx_iso, h_v_dx_iso);
  deep_copy(dev, *p_v_dy_iso, h_v_dy_iso);
  deep_copy(dev, *p_v_dz_iso, h_v_dz_iso);
#endif // ISO
  deep_copy(dev, *p_v_licomqice, h_v_licomqice);
#else  // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  new (p_v_atb) ViewDouble5D(&atb[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM+1, JMT, IMT);

  new (p_v_net) ViewDouble4D(&net[0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, JMT, IMT);

  new (p_v_at) ViewDouble5D(&at[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_pdensity) ViewDouble4D(&pdensity[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_amld) ViewDouble3D(&amld[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_tend) ViewDouble5D(&tend[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ax) ViewDouble5D(&ax[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ay) ViewDouble5D(&ay[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_az) ViewDouble5D(&az[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dx) ViewDouble5D(&dx[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dy) ViewDouble5D(&dy[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dz) ViewDouble5D(&dz[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_penetrate) ViewDouble4D(&penetrate[0][0][0][0],
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_dt_diff) ViewDouble5D(&dt_diff[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dt_conv) ViewDouble5D(&dt_conv[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

#ifdef ISO
  new (p_v_aay_iso) ViewDouble5D(&aay_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ddy_iso) ViewDouble5D(&ddy_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_ax_iso) ViewDouble5D(&ax_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_ay_iso) ViewDouble5D(&ay_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_az_iso) ViewDouble5D(&az_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);

  new (p_v_dx_iso) ViewDouble5D(&dx_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dy_iso) ViewDouble5D(&dy_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
  new (p_v_dz_iso) ViewDouble5D(&dz_iso[0][0][0][0][0],
      MAX_BLOCKS_CLINIC, NTRA, KM, JMT, IMT);
#endif // ISO

  new (p_v_licomqice) ViewDouble3D(&licomqice[0][0][0],
      MAX_BLOCKS_CLINIC, JMT, IMT);

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
  return ;
}
static void kokkos_init_tmp_var() {
  using namespace CppTmpVar;
  using namespace KokkosTmpVar;

  // READYT
  p_v_work1 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_work2 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));

  p_v_pp  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_ppa = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_ppb = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_ppc = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_alpha = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_beta  = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  // READYC
  p_v_wp12 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_wp13 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_riv1 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));
  p_v_riv2 = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

  p_v_wk1 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_wk2 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_wk3 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_wp3 = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
#ifdef BCKMEX
  ViewDouble3D *p_v_diff_back    = nullptr;
  ViewDouble3D *p_v_diff_back_sh = nullptr;
  ViewDouble3D *p_v_diff_back_nn = nullptr;
#endif // BCKMEX

  p_v_uv_ws_face = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

#ifndef SMAG
  p_v_div_out = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

#ifdef BIHAR
  p_v_curl = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_cc = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

  p_v_d2uk = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
  p_v_d2vk = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
#endif // BIHAR
#endif // SMAG

  // TRACER
  p_v_vtl_ori = (ViewDouble4D *) malloc(sizeof(ViewDouble4D)
      * IMT * JMT * KM * MAX_BLOCKS_CLINIC);

  p_v_adv_tt = (ViewDouble3D *) malloc(sizeof(ViewDouble3D)
      * IMT * JMT * KM);

  p_v_at_00_max_min = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

#ifdef BIHAR
  p_v_dt2k = (ViewDouble3D *) malloc(sizeof(ViewDouble3D));
#endif // BIHAR

  p_v_nn = (ViewInt1D *) malloc(sizeof(ViewInt1D));
  p_v_xs = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

  p_v_c_cnsew = (ViewDouble4D *) malloc(sizeof(ViewDouble4D));

// BAROTR
  p_v_gradx = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));
  p_v_grady = (ViewDouble2D *) malloc(sizeof(ViewDouble2D));

// POP Halo Update
  p_v_halo_buffer = (ViewDouble1D *) malloc(sizeof(ViewDouble1D));

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  new (p_v_work1) ViewDouble3D("pointer_view_work1",
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_work2) ViewDouble3D("pointer_view_work2",
      MAX_BLOCKS_CLINIC, JMT, IMT);

  new (p_v_pp)  ViewDouble4D("pointer_view_pp",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppa) ViewDouble4D("pointer_view_ppa",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppb) ViewDouble4D("pointer_view_ppb",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppc) ViewDouble4D("pointer_view_ppc",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_alpha) ViewDouble4D("pointer_view_alpha",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_beta)  ViewDouble4D("pointer_view_beta",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_wp12) ViewDouble4D("pointer_view_wp12",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp13) ViewDouble4D("pointer_view_wp13",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_riv1) ViewDouble4D("pointer_view_riv1",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_riv2) ViewDouble4D("pointer_view_riv2",
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_wk1) ViewDouble3D("pointer_view_wk1",
      KM, JMT, IMT);
  new (p_v_wk2) ViewDouble3D("pointer_view_wk2",
      KM, JMT, IMT);
  new (p_v_wk3) ViewDouble3D("pointer_view_wk3",
      KM, JMT, IMT);
  new (p_v_wp3) ViewDouble3D("pointer_view_wp3",
      KM, JMT, IMT);
#ifdef BCKMEX
  ViewDouble3D *p_v_diff_back    = nullptr;
  ViewDouble3D *p_v_diff_back_sh = nullptr;
  ViewDouble3D *p_v_diff_back_nn = nullptr;
#endif // BCKMEX

  new (p_v_uv_ws_face) ViewDouble4D("pointer_view_uv_ws_face", KM, JMT, IMT, 2);
#ifndef SMAG

  new (p_v_div_out)   ViewDouble2D("pointer_view_div_out", JMT, IMT);

#ifdef BIHAR
  new (p_v_curl) ViewDouble2D("pointer_view_curl", JMT, IMT);

  new (p_v_cc) ViewDouble2D("pointer_view_cc", JMT, IMT);
  new (p_v_d2uk) ViewDouble3D("pointer_view_d2uk", KM, JMT, IMT);
  new (p_v_d2vk) ViewDouble3D("pointer_view_d2vk", KM, JMT, IMT);
#endif // BIHAR
#endif // SMAG

  new (p_v_vtl_ori) ViewDouble4D("pointer_view_vtl_ori", 
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  new (p_v_adv_tt) ViewDouble3D("pointer_view_adv_tt", 
      KM, JMT, IMT);

  new (p_v_at_00_max_min) ViewDouble4D("pointer_view_at_00_max_min", 
      KM, JMT, IMT, 3);

#ifdef BIHAR
  new (p_v_dt2k) ViewDouble3D("pointer_view_dt2k", 
      KM, JMT, IMT);
#endif // BIHAR

  new (p_v_nn) ViewInt1D("pointer_view_nn", JMT);
  new (p_v_xs) ViewDouble1D("pointer_view_xs", IMT);

  new (p_v_c_cnsew) ViewDouble4D("pointer_view_c_cnsew", 
      KM, NY_BLOCK, NX_BLOCK, 5);
  new (p_v_gradx) ViewDouble2D("pointer_view_gradx", JMT, IMT);
  new (p_v_grady) ViewDouble2D("pointer_view_grady", JMT, IMT);
#ifdef ISO
#ifdef LDD97
  ViewDouble4D v_f1("view_f1", MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  ViewDouble4D v_f2("view_f2", MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // LDD97
#endif // ISO

  new (p_v_halo_buffer) ViewDouble1D("pointer_view_halo_buffer", 
      8 * KM * (IMT + JMT));
#else  // KOKKOS_ENABLE_DEVICE_MEM_SPACE

  // READYT
  work1 = new double[MAX_BLOCKS_CLINIC * JMT * IMT];
  work2 = new double[MAX_BLOCKS_CLINIC * JMT * IMT];
  new (p_v_work1) ViewDouble3D(work1,
      MAX_BLOCKS_CLINIC, JMT, IMT);
  new (p_v_work2) ViewDouble3D(work2,
      MAX_BLOCKS_CLINIC, JMT, IMT);

  pp  = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  ppa = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  ppb = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  ppc = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_pp)  ViewDouble4D(pp,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppa) ViewDouble4D(ppa,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppb) ViewDouble4D(ppb,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_ppc) ViewDouble4D(ppc,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  alpha = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  beta  = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_alpha) ViewDouble4D(alpha,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_beta)  ViewDouble4D(beta,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  wp12 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  wp13 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_wp12) ViewDouble4D(wp12,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_wp13) ViewDouble4D(wp13,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  riv1 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  riv2 = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_riv1) ViewDouble4D(riv1,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  new (p_v_riv2) ViewDouble4D(riv2,
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  double* wk1 = new double[KM * JMT * IMT];
  double* wk2 = new double[KM * JMT * IMT];
  double* wk3 = new double[KM * JMT * IMT];
  double* wp3 = new double[KM * JMT * IMT];
  new (p_v_wk1) ViewDouble3D(wk1, KM, JMT, IMT);
  new (p_v_wk2) ViewDouble3D(wk2, KM, JMT, IMT);
  new (p_v_wk3) ViewDouble3D(wk3, KM, JMT, IMT);
  new (p_v_wp3) ViewDouble3D(wp3, KM, JMT, IMT);

#ifdef BCKMEX
  ViewDouble3D *p_v_diff_back    = nullptr;
  ViewDouble3D *p_v_diff_back_sh = nullptr;
  ViewDouble3D *p_v_diff_back_nn = nullptr;
#endif // BCKMEX

  uv_ws_face = new double[MAX_BLOCKS_CLINIC * 2 * KM * JMT * IMT];
  new (p_v_uv_ws_face) ViewDouble4D(uv_ws_face, KM, JMT, IMT, 2);
#ifndef SMAG

  div_out = new double[JMT * IMT];
  new (p_v_div_out) ViewDouble2D(div_out, JMT, IMT);

#ifdef BIHAR
  curl   = new double[JMT * IMT];
  new (p_v_curl)   ViewDouble2D(curl,   JMT, IMT);

  cc = new double[JMT * IMT];
  new (p_v_cc) ViewDouble2D(cc, JMT, IMT);
  d2uk = new double[KM * JMT * IMT];
  d2vk = new double[KM * JMT * IMT];
  new (p_v_d2uk) ViewDouble3D(d2uk, KM, JMT, IMT);
  new (p_v_d2vk) ViewDouble3D(d2vk, KM, JMT, IMT);

#endif // BIHAR
#endif // SMAG

  vtl_ori = new double[MAX_BLOCKS_CLINIC * KM * JMT * IMT];
  new (p_v_vtl_ori) ViewDouble4D(vtl_ori, 
      MAX_BLOCKS_CLINIC, KM, JMT, IMT);

  adv_tt = new double[KM * JMT * IMT];
  new (p_v_adv_tt) ViewDouble3D(adv_tt, KM, JMT, IMT);

  at_00_max_min = new double[KM * JMT * IMT * 3];
  new (p_v_at_00_max_min) ViewDouble4D(at_00_max_min, 
      KM, JMT, IMT, 3);

#ifdef BIHAR
  dt2k = new double[KM * JMT * IMT];
  new (p_v_dt2k) ViewDouble2D(dt2k, KM, JMT, IMT);
#endif // BIHAR

  nn = new int[JMT];
  xs = new double[IMT];
  new (p_v_nn) ViewInt1D(nn, JMT);
  new (p_v_xs) ViewDouble1D(xs, IMT);

  c_cnsew = new double[5 * KM * NY_BLOCK * NX_BLOCK];
  new (p_v_c_cnsew) ViewDouble4D(c_cnsew, KM, NY_BLOCK, NX_BLOCK, 5);
#ifdef ISO
#ifdef LDD97
  ViewDouble4D v_f1(view_f1, MAX_BLOCKS_CLINIC, KM, JMT, IMT);
  ViewDouble4D v_f2(view_f2, MAX_BLOCKS_CLINIC, KM, JMT, IMT);
#endif // LDD97
#endif // ISO

  gradx = new double[JMT * IMT];
  grady = new double[JMT * IMT];
  new (p_v_gradx) ViewDouble2D(gradx, JMT, IMT);
  new (p_v_grady) ViewDouble2D(grady, JMT, IMT);

  halo_buffer = new double[8 * KM * (IMT + JMT)];
  new (p_v_halo_buffer) ViewDouble1D(halo_buffer, 
      8 * KM * (IMT + JMT));
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE
  return ;
}

#endif // LICOM_ENABLE_KOKKOS
