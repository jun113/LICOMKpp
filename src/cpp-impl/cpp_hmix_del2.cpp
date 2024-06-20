#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/fortran_hmix_del2.h"
namespace CppHmixDel2 {
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&dtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dtn_;
double (&dts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dts_;
double (&dte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dte_;
double (&dtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dtw_;
double (&duc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_duc_;
double (&dun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dun_;
double (&dus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dus_;
double (&due)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_due_;
double (&duw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_duw_;
double (&dmc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmc_;
double (&dmn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmn_;
double (&dms)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dms_;
double (&dme)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dme_;
double (&dmw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dmw_;
double (&dum)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_dum_;
double (&ahf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_ahf_;
double (&amf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = hmix_del2_mp_amf_;

double &ah = hmix_del2_mp_ah_;
double &am = hmix_del2_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&dtn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dtn;
double (&dts)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dts;
double (&dte)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dte;
double (&dtw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dtw;
double (&duc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_duc;
double (&dun)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dun;
double (&dus)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dus;
double (&due)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_due;
double (&duw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_duw;
double (&dmc)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmc;
double (&dmn)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmn;
double (&dms)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dms;
double (&dme)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dme;
double (&dmw)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dmw;
double (&dum)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_dum;
double (&ahf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_ahf;
double (&amf)[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK] = __hmix_del2_MOD_amf;

double &ah = __hmix_del2_MOD_ah;
double &am = __hmix_del2_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

void hdiffu_del2 (const int &k,
    double (&hduk)[NY_BLOCK][NX_BLOCK],
    double (&hdvk)[NY_BLOCK][NX_BLOCK],
    const double (&umixk)[NY_BLOCK][NX_BLOCK],
    const double (&vmixk)[NY_BLOCK][NX_BLOCK]) {

  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;

  using CppConstantMod::C0;

  using CppGrid::kmu;
  using CppPconstMod::viv;

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hduk[j][i] = C0;
      hdvk[j][i] = C0;
    }
  }
  const int bid = 0;
  // const int ib = this_block.ib;
  // const int ie = this_block.ie;
  // const int jb = this_block.jb;
  // const int je = this_block.je;
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      const double cc = duc[bid][j][i] + dum[bid][j][i];
    
      hduk[j][i] = am * ((cc * umixk[j  ][i  ] +
              dun[bid][j][i] * umixk[j-1][i  ] +
              dus[bid][j][i] * umixk[j+1][i  ] +
              due[bid][j][i] * umixk[j  ][i+1] +
              duw[bid][j][i] * umixk[j  ][i-1]) +
             (dmc[bid][j][i] * vmixk[j  ][i  ] +
              dmn[bid][j][i] * vmixk[j-1][i  ] +
              dms[bid][j][i] * vmixk[j+1][i  ] +
              dme[bid][j][i] * vmixk[j  ][i+1] +
              dmw[bid][j][i] * vmixk[j  ][i-1])) *
                  viv[bid][k][j][i];

      hdvk[j][i] = am * ((cc * vmixk[j  ][i  ] +
              dun[bid][j][i] * vmixk[j-1][i  ] +
              dus[bid][j][i] * vmixk[j+1][i  ] +
              due[bid][j][i] * vmixk[j  ][i+1] +
              duw[bid][j][i] * vmixk[j  ][i-1]) -
             (dmc[bid][j][i] * umixk[j  ][i  ] +
              dmn[bid][j][i] * umixk[j-1][i  ] +
              dms[bid][j][i] * umixk[j+1][i  ] +
              dme[bid][j][i] * umixk[j  ][i+1] +
              dmw[bid][j][i] * umixk[j  ][i-1])) *
                  viv[bid][k][j][i];
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k > kmu[bid][j][i] - 1) {
        hduk[j][i] = C0;
        hdvk[j][i] = C0;
      }
    }
  }
  return ;
}

void hdifft_del2 (const int &k, 
    double (&hdtk)[NY_BLOCK][NX_BLOCK],
    const double (&tmix)[NY_BLOCK][NX_BLOCK]) {

  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;

  using CppGrid::kmt;
  using CppGrid::kmtn;
  using CppGrid::kmts;
  using CppGrid::kmte;
  using CppGrid::kmtw;

  using CppConstantMod::C0;

  const int bid = 0;
  const int kk = k + 1;

  double cc[NY_BLOCK][NX_BLOCK];
  double cn[NY_BLOCK][NX_BLOCK];
  double cs[NY_BLOCK][NX_BLOCK];
  double ce[NY_BLOCK][NX_BLOCK];
  double cw[NY_BLOCK][NX_BLOCK];

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cn[j][i] = 
          (kk <= kmtn[bid][j][i] && kk <= kmt[bid][j][i]) 
              ? dtn[bid][j][i] : C0;
      cs[j][i] = 
          (kk <= kmts[bid][j][i] && kk <= kmt[bid][j][i]) 
              ? dts[bid][j][i] : C0;
      ce[j][i] = 
          (kk <= kmte[bid][j][i] && kk <= kmt[bid][j][i]) 
              ? dte[bid][j][i] : C0;
      cw[j][i] = 
          (kk <= kmtw[bid][j][i] && kk <= kmt[bid][j][i]) 
              ? dtw[bid][j][i] : C0;

    }
  }

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cc[j][i] = - (cn[j][i] + cs[j][i] + ce[j][i] + cw[j][i]);
    }
  }

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hdtk[j][i] = C0;
    }
  }

  // const int ib = this_block.ib;
  // const int ie = this_block.ie;
  // const int jb = this_block.jb;
  // const int je = this_block.je;

  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdtk[j][i] = ah * (cc[j][i] * tmix[j  ][i  ]
                       + cn[j][i] * tmix[j-1][i  ]
                       + cs[j][i] * tmix[j+1][i  ]
                       + ce[j][i] * tmix[j  ][i+1]
                       + cw[j][i] * tmix[j  ][i-1]);
    }
  }
  return ;
}

} // CppHmixDel2
#endif // LICOM_ENABLE_FORTRAN
