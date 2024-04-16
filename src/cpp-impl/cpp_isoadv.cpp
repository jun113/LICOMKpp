#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#ifdef ISO

#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_isopyc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

#include <cmath>

using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::KM;
using CppParamMod::MAX_BLOCKS_CLINIC;

void isoadv (const double (&k1)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
		const double (&k2)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  			double (&adv_vetiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]) {
  
	using CppDomain::nblocks_clinic;
  using CppGrid::hts;
  using CppGrid::htw;
  using CppGrid::kmt;
  using CppGrid::tarea_r;
  // using CppIsopycMod::k1;
  // using CppIsopycMod::k2;
  using CppIsopycMod::dzr;
  using CppIsopycMod::tmask;
  using CppIsopycMod::athkdf;
  using CppIsopycMod::adv_vntiso;
  using CppIsopycMod::adv_vbtiso;
  // using CppIsopycMod::adv_vetiso;
  using CppPconstMod::dzp;

  const double c0 = 0.0;
  const double p5 = 0.5;

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < IMT; ++i) {
          adv_vntiso[iblock][j][k][i] = c0;
          adv_vetiso[iblock][j][k][i] = c0;
        }
      }
    }
  }
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int k = 0; k < KM+1; ++k) {
        for (int i = 0; i < IMT; ++i) {
          adv_vbtiso[iblock][j][k][i] = c0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 1; k < KM-1; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double fxa = - p5 * dzr[k] * athkdf[iblock][j][i];
          adv_vntiso[iblock][j][k][i] = fxa 
              * tmask[iblock][j][k][i] * tmask[iblock][j+1][k][i] 
              * (k2[iblock][0][j][k][i] - k2[iblock][0][j][k+2][i]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const double fxa = - p5 * dzr[0] * athkdf[iblock][j][i];
        adv_vntiso[iblock][j][0][i] = - fxa 
            * tmask[iblock][j][0][i] * tmask[iblock][j+1][0][i]
                * (k2[iblock][0][j][1][i] + k2[iblock][0][j][2][i]);
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j+1][i]) - 1;
        if (k != -1) {
          adv_vntiso[iblock][j][k][i] = - p5 * dzr[k] * athkdf[iblock][j][i]
              * tmask[iblock][j][k][i] * tmask[iblock][j+1][k][i]
                  * (k2[iblock][0][j][k+1][i] + k2[iblock][0][j][k][i]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 1; k < KM-1; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double fxa = - p5 * dzr[k] * athkdf[iblock][j][i];
          adv_vetiso[iblock][j][k][i] = fxa 
              * tmask[iblock][j][k][i] * tmask[iblock][j][k][i+1]
                  * (k1[iblock][0][j][k][i] - k1[iblock][0][j][k+2][i]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const double fxa = - p5 * dzr[0] * athkdf[iblock][j][i];
        adv_vetiso[iblock][j][0][i] = - fxa 
            * tmask[iblock][j][0][i] * tmask[iblock][j][0][i+1]
                * (k1[iblock][0][j][1][i] + k1[iblock][0][j][2][i]); 
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j][i+1]) - 1;
        if (k != -1) {
          adv_vetiso[iblock][j][k][i] = - p5 * dzr[k] * athkdf[iblock][j][i]
              * tmask[iblock][j][k][i] * tmask[iblock][j][k][i+1]
                  * (k1[iblock][0][j][k+1][i] + k1[iblock][0][j][k][i]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 2; j < JMT-2; ++j) {
      for (int k = 0; k < KM-1; ++k) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_vbtiso[iblock][j][k+1][i] = dzp[k]
              * ((adv_vetiso[iblock][j  ][k][i  ] * htw[iblock][j  ][i+1]
                - adv_vetiso[iblock][j  ][k][i-1] * htw[iblock][j  ][i  ])
                    * tarea_r[iblock][j][i]
               + (adv_vntiso[iblock][j  ][k][i  ] * hts[iblock][j  ][i  ]
                - adv_vntiso[iblock][j-1][k][i  ] * hts[iblock][j-1][i  ])
                    * tarea_r[iblock][j][i]);
        }
      }
    }
  }

  // Note: adv_vbtiso(imt, 0:km, jmt, max_blocks_clinic)
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 2; j < JMT-2; ++j) {
      for (int k = 1; k <= KM-1; ++k) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_vbtiso[iblock][j][k][i] += adv_vbtiso[iblock][j][k-1][i];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {
        adv_vbtiso[iblock][j][kmt[iblock][j][i]][i] = c0;
      }
    }
  }

  return ;
}

#endif // ISO

#endif // LICOM_ENABLE_FORTRAN