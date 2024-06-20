#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_extern_functions.h"

#include "../head/fortran_hmix_del4.h"

#include <vector>
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

void hdiffu_del4 (const int &k,
    double (&hduk)[NY_BLOCK][NX_BLOCK],
    double (&hdvk)[NY_BLOCK][NX_BLOCK],
    const double (&umixk)[NY_BLOCK][NX_BLOCK],
    const double (&vmixk)[NY_BLOCK][NX_BLOCK]) {

  using CppBlocks::ib;
  using CppBlocks::ie;
  using CppBlocks::jb;
  using CppBlocks::je;
  
  using CppConstantMod::C0;
  using CppDomain::nblocks_clinic;
  using CppGrid::kmu;
  using CppGrid::uarea;

  //const int bid = this_block.local_id - 1;
  const int bid = 0;
  // TODO this_block

  std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
      am_factor(nblocks_clinic);

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        am_factor[iblock][j][i] = 1.0;
      }
    }
  }
  double div_out[NY_BLOCK][NX_BLOCK];
  div(k, div_out, umixk, vmixk);  

  double curl[NY_BLOCK][NX_BLOCK];
  zcurl(k, curl, umixk, vmixk);

  double gradx1[NY_BLOCK][NX_BLOCK];
  double grady1[NY_BLOCK][NX_BLOCK];
  grad(k, gradx1, grady1, curl);

  double gradx2[NY_BLOCK][NX_BLOCK];
  double grady2[NY_BLOCK][NX_BLOCK];
  grad(k, gradx2, grady2, div_out);

  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      double dxdy = std::pow(std::sqrt(uarea[bid][j][i]), 5) * 45.0;
      am_factor[bid][j][i] = std::sqrt(
          std::pow(gradx1[j][i], 2) + std::pow(gradx2[j][i], 2) + 
          std::pow(grady1[j][i], 2) + std::pow(grady2[j][i], 2)) *
          dxdy / std::abs(am * amf[bid][j][i]);
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK; ++j) {
      for (int i = 0; i < NX_BLOCK; ++i) {
        am_factor[iblock][j][i] = std::min(40.0, am_factor[iblock][j][i]);
        am_factor[iblock][j][i] = std::max(1.0,  am_factor[iblock][j][i]);
      }
    }
  }
  double cc[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      cc[j][i] = duc[bid][j][i] + dum[bid][j][i];
    }
  }
  double d2uk[NY_BLOCK][NX_BLOCK];
  double d2vk[NY_BLOCK][NX_BLOCK];
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      d2uk[j][i] = C0;
      d2vk[j][i] = C0;
    }
  }
  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2uk[j][i] = (cc[j][i] * umixk[j  ][i  ] +
              dun[bid][j][i] * umixk[j-1][i  ] +
              dus[bid][j][i] * umixk[j+1][i  ] +
              due[bid][j][i] * umixk[j  ][i+1] +
              duw[bid][j][i] * umixk[j  ][i-1]) +
             (dmc[bid][j][i] * vmixk[j  ][i  ] +
              dmn[bid][j][i] * vmixk[j-1][i  ] +
              dms[bid][j][i] * vmixk[j+1][i  ] +
              dme[bid][j][i] * vmixk[j  ][i+1] +
              dmw[bid][j][i] * vmixk[j  ][i-1]);
    }
  }
  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2vk[j][i] = (cc[j][i] * vmixk[j  ][i  ] +
              dun[bid][j][i] * vmixk[j-1][i  ] +
              dus[bid][j][i] * vmixk[j+1][i  ] +
              due[bid][j][i] * vmixk[j  ][i+1] +
              duw[bid][j][i] * vmixk[j  ][i-1]) -
             (dmc[bid][j][i] * umixk[j  ][i  ] +
              dmn[bid][j][i] * umixk[j-1][i  ] +
              dms[bid][j][i] * umixk[j+1][i  ] +
              dme[bid][j][i] * umixk[j  ][i+1] +
              dmw[bid][j][i] * umixk[j  ][i-1]);
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k <= kmu[bid][j][i]-1) {
        d2uk[j][i] = am_factor[bid][j][i] *
            amf[bid][j][i] * d2uk[j][i];
        d2vk[j][i] = am_factor[bid][j][i] *
            amf[bid][j][i] * d2vk[j][i];
      } else {
        d2uk[j][i] = C0;
        d2vk[j][i] = C0;
      }
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hduk[j][i] = C0;
      hdvk[j][i] = C0;
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hduk[j][i] = am * ((cc[j][i] * d2uk[j  ][i  ] +
                    dun[bid][j][i] * d2uk[j-1][i  ] +
                    dus[bid][j][i] * d2uk[j+1][i  ] +
                    due[bid][j][i] * d2uk[j  ][i+1] +
                    duw[bid][j][i] * d2uk[j  ][i-1]) +
                   (dmc[bid][j][i] * d2vk[j  ][i  ] +
                    dmn[bid][j][i] * d2vk[j-1][i  ] +
                    dms[bid][j][i] * d2vk[j+1][i  ] +
                    dme[bid][j][i] * d2vk[j  ][i+1] +
                    dmw[bid][j][i] * d2vk[j  ][i-1]));
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdvk[j][i] = am * ((cc[j][i] * d2vk[j  ][i  ] +
                    dun[bid][j][i] * d2vk[j-1][i  ] +
                    dus[bid][j][i] * d2vk[j+1][i  ] +
                    due[bid][j][i] * d2vk[j  ][i+1] +
                    duw[bid][j][i] * d2vk[j  ][i-1]) +
                   (dmc[bid][j][i] * d2uk[j  ][i  ] +
                    dmn[bid][j][i] * d2uk[j-1][i  ] +
                    dms[bid][j][i] * d2uk[j+1][i  ] +
                    dme[bid][j][i] * d2uk[j  ][i+1] +
                    dmw[bid][j][i] * d2uk[j  ][i-1]));
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      if (k > kmu[bid][j][i]-1) {
        hduk[j][i] = C0;
        hdvk[j][i] = C0;
      }
    }
  }
  return ;
}

// void hdiffu_del4(
//     const int &k,
//     double (&hduk)[NY_BLOCK][NX_BLOCK],
//     double (&hdvk)[NY_BLOCK][NX_BLOCK],
//     const double (&umixk)[NY_BLOCK][NX_BLOCK],
//     const double (&vmixk)[NY_BLOCK][NX_BLOCK],
//     const block &this_block) {
  
//   using CppConstantMod::C0;
//   using CppConstantMod::P5;
//   using CppGrid::dxu;
//   using CppGrid::dyu;
//   using CppGrid::dxur;
//   using CppGrid::dyur;
//   using CppGrid::htw;
//   using CppGrid::hts;
//   using CppGrid::kmt;
//   using CppGrid::kmu;
//   using CppGrid::uarea;
//   using CppGrid::tarea_r;

//   using CppHmixDel4::am;
//   using CppHmixDel4::amf;
//   using CppHmixDel4::duc;
//   using CppHmixDel4::due;
//   using CppHmixDel4::dum;
//   using CppHmixDel4::dun;
//   using CppHmixDel4::dus;
//   using CppHmixDel4::duw;
//   using CppHmixDel4::dmc;
//   using CppHmixDel4::dme;
//   using CppHmixDel4::dmn;
//   using CppHmixDel4::dms;
//   using CppHmixDel4::dmw;
//   //const int bid = this_block.local_id - 1;
//   const int bid = 0;
//   // TODO this_block

//   double div_out[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       div_out[j][i] = C0;
//     }
//   }
//   for (int j = 1; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK-1; ++i) {
//       if (k <= kmt[bid][j][i]-1) {
//         div_out[j][i] = P5 * (
//             (umixk[j  ][i+1] + umixk[j-1][i+1] * htw[bid][j  ][i+1]) -
//             (umixk[j  ][i  ] + umixk[j-1][i  ] * htw[bid][j  ][i  ]) +
//             (vmixk[j  ][i+1] + vmixk[j  ][i  ] * hts[bid][j  ][i  ]) -
//             (vmixk[j-1][i+1] + vmixk[j-1][i  ] * hts[bid][j-1][i  ])) *
//                 tarea_r[bid][j][i];
//       }
//     }
//   }
//   double curl[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       curl[j][i] = C0;
//     }
//   }
//   for (int j = 1; j < NY_BLOCK; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       if (k <= kmt[bid][j][i]-1) {
//         curl[j][i] = P5 * (
//             vmixk[j  ][i  ] * dyu[bid][j  ][i  ] +
//             vmixk[j-1][i  ] * dyu[bid][j-1][i  ] -
//             vmixk[j  ][i-1] * dyu[bid][j  ][i-1] -
//             vmixk[j-1][i-1] * dyu[bid][j-1][i-1] -
//             umixk[j  ][i  ] * dxu[bid][j  ][i  ] -
//             umixk[j  ][i-1] * dxu[bid][j  ][i-1] +
//             umixk[j-1][i  ] * dxu[bid][j-1][i  ] +
//             umixk[j-1][i-1] * dxu[bid][j-1][i-1]);
//       }
//     }
//   }
//   double gradx1[NY_BLOCK][NX_BLOCK];
//   double grady1[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       gradx1[j][i] = C0;
//       grady1[j][i] = C0;
//     }
//   }
//   for (int j = 0; j < NY_BLOCK-1; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       gradx1[j][i] = dxur[bid][j][i] * P5 * (
//           curl[j+1][i  ] - curl[j  ][i-1] -
//           curl[j+1][i-1] + curl[j  ][i  ]);
//       grady1[j][i] = dyur[bid][j][i] * P5 * (
//           curl[j+1][i  ] - curl[j  ][i-1] +
//           curl[j+1][i-1] - curl[j  ][i  ]);
//     }
//   }
//   double gradx2[NY_BLOCK][NX_BLOCK];
//   double grady2[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       gradx2[j][i] = C0;
//       grady2[j][i] = C0;
//     }
//   }
//   for (int j = 0; j < NY_BLOCK-1; ++j) {
//     for (int i = 1; i < NX_BLOCK; ++i) {
//       gradx2[j][i] = dxur[bid][j][i] * P5 * (
//           div_out[j+1][i  ] - div_out[j  ][i-1] -
//           div_out[j+1][i-1] + div_out[j  ][i  ]);
//       grady2[j][i] = dyur[bid][j][i] * P5 * (
//           div_out[j+1][i  ] - div_out[j  ][i-1] +
//           div_out[j+1][i-1] - div_out[j  ][i  ]);
//     }
//   }
//   const int ib = this_block.ib;
//   const int ie = this_block.ie;
//   const int jb = this_block.jb;
//   const int je = this_block.je;

//   std::vector<std::array<std::array<double, NX_BLOCK>, NY_BLOCK>> 
//       am_factor(nblocks_clinic);

//   for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
//     for (int j = 0; j < NY_BLOCK; ++j) {
//       for (int i = 0; i < NX_BLOCK; ++i) {
//         am_factor[iblock][j][i] = 1.0;
//       }
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       double dxdy = pow(sqrt(uarea[bid][j][i]), 5) * 45.0;
//       am_factor[bid][j][i] = sqrt(
//           pow(gradx1[j][i], 2) + pow(gradx2[j][i], 2) + 
//           pow(grady1[j][i], 2) + pow(grady2[j][i], 2)) *
//           dxdy / fabs(am * amf[bid][j][i]);
//     }
//   }
//   for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
//     for (int j = 0; j < NY_BLOCK; ++j) {
//       for (int i = 0; i < NX_BLOCK; ++i) {
//         am_factor[iblock][j][i] = fmin(40.0, am_factor[iblock][j][i]);
//         am_factor[iblock][j][i] = fmax(1.0,  am_factor[iblock][j][i]);
//       }
//     }
//   }
//   double cc[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       cc[j][i] = duc[bid][j][i] + dum[bid][j][i];
//     }
//   }
//   double d2uk[NY_BLOCK][NX_BLOCK];
//   double d2vk[NY_BLOCK][NX_BLOCK];
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       d2uk[j][i] = C0;
//       d2vk[j][i] = C0;
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       d2uk[j][i] = (cc[j][i] * umixk[j  ][i  ] +
//               dun[bid][j][i] * umixk[j-1][i  ] +
//               dus[bid][j][i] * umixk[j+1][i  ] +
//               due[bid][j][i] * umixk[j  ][i+1] +
//               duw[bid][j][i] * umixk[j  ][i-1]) +
//              (dmc[bid][j][i] * vmixk[j  ][i  ] +
//               dmn[bid][j][i] * vmixk[j-1][i  ] +
//               dms[bid][j][i] * vmixk[j+1][i  ] +
//               dme[bid][j][i] * vmixk[j  ][i+1] +
//               dmw[bid][j][i] * vmixk[j  ][i-1]);
//     }
//   }
//   for (int j = jb-2; j < je+1; ++j) {
//     for (int i = ib-2; i < ie+1; ++i) {
//       d2vk[j][i] = (cc[j][i] * vmixk[j  ][i  ] +
//               dun[bid][j][i] * vmixk[j-1][i  ] +
//               dus[bid][j][i] * vmixk[j+1][i  ] +
//               due[bid][j][i] * vmixk[j  ][i+1] +
//               duw[bid][j][i] * vmixk[j  ][i-1]) +
//              (dmc[bid][j][i] * umixk[j  ][i  ] +
//               dmn[bid][j][i] * umixk[j-1][i  ] +
//               dms[bid][j][i] * umixk[j+1][i  ] +
//               dme[bid][j][i] * umixk[j  ][i+1] +
//               dmw[bid][j][i] * umixk[j  ][i-1]);
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       if (k <= kmu[bid][j][i]-1) {
//         d2uk[j][i] = am_factor[bid][j][i] *
//             amf[bid][j][i] * d2uk[j][i];
//         d2vk[j][i] = am_factor[bid][j][i] *
//             amf[bid][j][i] * d2vk[j][i];
//       } else {
//         d2uk[j][i] = C0;
//         d2vk[j][i] = C0;
//       }
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       hduk[j][i] = C0;
//       hdvk[j][i] = C0;
//     }
//   }
//   for (int j = jb-1; j < je; ++j) {
//     for (int i = ib-1; i < ie; ++i) {
//       hduk[j][i] = am * ((cc[j][i] * d2uk[j  ][i  ] +
//                     dun[bid][j][i] * d2uk[j-1][i  ] +
//                     dus[bid][j][i] * d2uk[j+1][i  ] +
//                     due[bid][j][i] * d2uk[j  ][i+1] +
//                     duw[bid][j][i] * d2uk[j  ][i-1]) +
//                    (dmc[bid][j][i] * d2vk[j  ][i  ] +
//                     dmn[bid][j][i] * d2vk[j-1][i  ] +
//                     dms[bid][j][i] * d2vk[j+1][i  ] +
//                     dme[bid][j][i] * d2vk[j  ][i+1] +
//                     dmw[bid][j][i] * d2vk[j  ][i-1]));
//     }
//   }
//   for (int j = jb-1; j < je; ++j) {
//     for (int i = ib-1; i < ie; ++i) {
//       hdvk[j][i] = am * ((cc[j][i] * d2vk[j  ][i  ] +
//                     dun[bid][j][i] * d2vk[j-1][i  ] +
//                     dus[bid][j][i] * d2vk[j+1][i  ] +
//                     due[bid][j][i] * d2vk[j  ][i+1] +
//                     duw[bid][j][i] * d2vk[j  ][i-1]) +
//                    (dmc[bid][j][i] * d2uk[j  ][i  ] +
//                     dmn[bid][j][i] * d2uk[j-1][i  ] +
//                     dms[bid][j][i] * d2uk[j+1][i  ] +
//                     dme[bid][j][i] * d2uk[j  ][i+1] +
//                     dmw[bid][j][i] * d2uk[j  ][i-1]));
//     }
//   }
//   for (int j = 0; j < NY_BLOCK; ++j) {
//     for (int i = 0; i < NX_BLOCK; ++i) {
//       if (k > kmu[bid][j][i]-1) {
//         hduk[j][i] = C0;
//         hdvk[j][i] = C0;
//       }
//     }
//   }
//   return ;
// }
void hdifft_del4(const int &k, 
    double (&d2tk)[NY_BLOCK][NX_BLOCK],
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
  double cc[NY_BLOCK][NX_BLOCK];
  double cn[NY_BLOCK][NX_BLOCK];
  double cs[NY_BLOCK][NX_BLOCK];
  double ce[NY_BLOCK][NX_BLOCK];
  double cw[NY_BLOCK][NX_BLOCK];

  const int kk = k + 1;
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
      cc[j][i] = -(cn[j][i] + cs[j][i] + ce[j][i] + cw[j][i]);
    }
  }

  for (int j = jb-2; j < je+1; ++j) {
    for (int i = ib-2; i < ie+1; ++i) {
      d2tk[j][i] = ahf[bid][j][i] * (
          cc[j][i] * tmix[j  ][i  ]
        + cn[j][i] * tmix[j-1][i  ]
        + cs[j][i] * tmix[j+1][i  ]
        + ce[j][i] * tmix[j  ][i+1]
        + cw[j][i] * tmix[j  ][i-1]);
    }
  }
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      hdtk[j][i] = C0;
    }
  }
  for (int j = jb-1; j < je; ++j) {
    for (int i = ib-1; i < ie; ++i) {
      hdtk[j][i] = ah * (cc[j][i] * d2tk[j  ][i  ]
                       + cn[j][i] * d2tk[j-1][i  ]
                       + cs[j][i] * d2tk[j+1][i  ]
                       + ce[j][i] * d2tk[j  ][i+1]
                       + cw[j][i] * d2tk[j  ][i-1]);
    }
  }
  return ;
}

} // CppHmixDel4
#endif // LICOM_ENABLE_FORTRAN
