#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_work_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
void upwell (
    const double (&uwk)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&vwk)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&h0wk)[MAX_BLOCKS_CLINIC][JMT][IMT]){

  using CppConstantMod::C0;
  using CppConstantMod::P5;
	using CppDomain::nblocks_clinic;
  using CppDynMod::ws;
  // using CppGrid::au0;
  // using CppGrid::aus;
  // using CppGrid::auw;
  // using CppGrid::ausw;
  using CppGrid::kmt;
  using CppGrid::hts;
  using CppGrid::htw;
  using CppGrid::tarea_r;
	using CppParamMod::NX_BLOCK;
	using CppParamMod::NY_BLOCK;
  using CppPconstMod::dzp;
  using CppPconstMod::vit;
  using CppPconstMod::ohbu;
  using CppPconstMod::ohbt;
  using CppWorkMod::wka;
  using CppWorkMod::work;

	double uk[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
	double vk[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wka[iblock][k][j][i] = 0.0;
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK-1; ++j) {
      for (int i = 1; i < NX_BLOCK; ++i) {
        // work[iblock][j][i] = 
        //     au0 [iblock][j][i] * h0wk[iblock][j  ][i  ] +
        //     aus [iblock][j][i] * h0wk[iblock][j+1][i  ] +
        //     auw [iblock][j][i] * h0wk[iblock][j  ][i-1] +
        //     ausw[iblock][j][i] * h0wk[iblock][j+1][i-1];
        work[iblock][j][i] = 
            CppConstantMod::P25 * h0wk[iblock][j  ][i  ] +
            CppConstantMod::P25 * h0wk[iblock][j+1][i  ] +
            CppConstantMod::P25 * h0wk[iblock][j  ][i-1] +
            CppConstantMod::P25 * h0wk[iblock][j+1][i-1];
      }
    }
    for (int i = 0; i < NX_BLOCK; ++i) {
      work[iblock][NY_BLOCK-1][i] = C0;
    }
    for (int j = 0; j < NY_BLOCK; ++j) {
      work[iblock][j][0] = C0;
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          uk[iblock][k][j][i] = (1.0 + work[iblock][j][i] * 
              ohbu[iblock][j][i]) * uwk[iblock][k][j][i];
          vk[iblock][k][j][i] = (1.0 + work[iblock][j][i] * 
              ohbu[iblock][j][i]) * vwk[iblock][k][j][i];
        }
      }
    }
  }
  const int bid = 0;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < NY_BLOCK; ++j) {
        for (int i = 0; i < NX_BLOCK; ++i) {
          wka[iblock][k][j][i] = C0;
        }
      }
      for (int j = 1; j < NY_BLOCK; ++j) {
        for (int i = 0; i < NX_BLOCK-1; ++i) {
          if (k <= kmt[bid][j][i]-1) {
            wka[iblock][k][j][i] = P5 * (
                (uk[iblock][k][j  ][i+1] + uk[iblock][k][j-1][i+1]) * htw[bid][j  ][i+1] -
                (uk[iblock][k][j  ][i  ] + uk[iblock][k][j-1][i  ]) * htw[bid][j  ][i  ] +
                (vk[iblock][k][j  ][i+1] + vk[iblock][k][j  ][i  ]) * hts[bid][j  ][i  ] -
                (vk[iblock][k][j-1][i+1] + vk[iblock][k][j-1][i  ]) * hts[bid][j-1][i  ]) *
              tarea_r[bid][j][i];
          }
        }
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work[iblock][j][i] -= dzp[k] * 
              wka[iblock][k][j][i] * vit[iblock][k][j][i];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ws[iblock][k][j][i] = vit[iblock][k][j][i] *
              (ws[iblock][k-1][j][i] + dzp[k-1] * 
              (work[iblock][j][i] * ohbt[iblock][j][i] + 
                   wka[iblock][k-1][j][i]));
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        work[iblock][j][i] = 1.0 / (1.0 + h0wk[iblock][j][i] *
            ohbt[iblock][j][i]);
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ws[iblock][k][j][i] *= work[iblock][j][i];
        }
      }
    }
  }
  return ;
}
void upwell(
    const double (*uwk)[KM][JMT][IMT], 
    const double (*vwk)[KM][JMT][IMT], 
    const double (*h0wk)[JMT][IMT]){

  using CppConstantMod::C0;
  using CppConstantMod::P5;
	using CppDomain::nblocks_clinic;
  using CppDynMod::ws;
  // using CppGrid::au0;
  // using CppGrid::aus;
  // using CppGrid::auw;
  // using CppGrid::ausw;
  using CppGrid::kmt;
  using CppGrid::hts;
  using CppGrid::htw;
  using CppGrid::tarea_r;
	using CppParamMod::NX_BLOCK;
	using CppParamMod::NY_BLOCK;
  using CppPconstMod::dzp;
  using CppPconstMod::vit;
  using CppPconstMod::ohbu;
  using CppPconstMod::ohbt;
  using CppWorkMod::wka;
  using CppWorkMod::work;

	double uk[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
	double vk[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wka[iblock][k][j][i] = 0.0;
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < NY_BLOCK-1; ++j) {
      for (int i = 1; i < NX_BLOCK; ++i) {
        // work[iblock][j][i] = 
        //     au0 [iblock][j][i] * h0wk[iblock][j  ][i  ] +
        //     aus [iblock][j][i] * h0wk[iblock][j+1][i  ] +
        //     auw [iblock][j][i] * h0wk[iblock][j  ][i-1] +
        //     ausw[iblock][j][i] * h0wk[iblock][j+1][i-1];
        work[iblock][j][i] = 
            CppConstantMod::P25 * h0wk[iblock][j  ][i  ] +
            CppConstantMod::P25 * h0wk[iblock][j+1][i  ] +
            CppConstantMod::P25 * h0wk[iblock][j  ][i-1] +
            CppConstantMod::P25 * h0wk[iblock][j+1][i-1];
      }
    }
    for (int i = 0; i < NX_BLOCK; ++i) {
      work[iblock][NY_BLOCK-1][i] = C0;
    }
    for (int j = 0; j < NY_BLOCK; ++j) {
      work[iblock][j][0] = C0;
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          uk[iblock][k][j][i] = (1.0 + work[iblock][j][i] * 
              ohbu[iblock][j][i]) * uwk[iblock][k][j][i];
          vk[iblock][k][j][i] = (1.0 + work[iblock][j][i] * 
              ohbu[iblock][j][i]) * vwk[iblock][k][j][i];
        }
      }
    }
  }
  const int bid = 0;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < NY_BLOCK; ++j) {
        for (int i = 0; i < NX_BLOCK; ++i) {
          wka[iblock][k][j][i] = C0;
        }
      }
      for (int j = 1; j < NY_BLOCK; ++j) {
        for (int i = 0; i < NX_BLOCK-1; ++i) {
          if (k <= kmt[bid][j][i]-1) {
            wka[iblock][k][j][i] = P5 * (
                (uk[iblock][k][j  ][i+1] + uk[iblock][k][j-1][i+1]) * htw[bid][j  ][i+1] -
                (uk[iblock][k][j  ][i  ] + uk[iblock][k][j-1][i  ]) * htw[bid][j  ][i  ] +
                (vk[iblock][k][j  ][i+1] + vk[iblock][k][j  ][i  ]) * hts[bid][j  ][i  ] -
                (vk[iblock][k][j-1][i+1] + vk[iblock][k][j-1][i  ]) * hts[bid][j-1][i  ]) *
              tarea_r[bid][j][i];
          }
        }
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work[iblock][j][i] -= dzp[k] * 
              wka[iblock][k][j][i] * vit[iblock][k][j][i];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ws[iblock][k][j][i] = vit[iblock][k][j][i] *
              (ws[iblock][k-1][j][i] + dzp[k-1] * 
              (work[iblock][j][i] * ohbt[iblock][j][i] + 
                   wka[iblock][k-1][j][i]));
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        work[iblock][j][i] = 1.0 / (1.0 + h0wk[iblock][j][i] *
            ohbt[iblock][j][i]);
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ws[iblock][k][j][i] *= work[iblock][j][i];
        }
      }
    }
  }
  return ;
}
#endif // LICOM_ENABLE_FORTRAN