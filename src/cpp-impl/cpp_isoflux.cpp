#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#ifdef ISO

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_isopyc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"

using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::KM;
using CppParamMod::MAX_BLOCKS_CLINIC;

void isoflux(const int &mtrace, 
    const double (&k1)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
    const double (&k2)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
    const double (&k3)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT],
    const double (&adv_vetiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT],
        double (&tf)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]) {

	using CppConstantMod::C0;
	using CppConstantMod::P5;
	using CppConstantMod::P25;
	using CppDomain::nblocks_clinic;
  using CppGrid::hue;
  using CppGrid::hun;
  using CppGrid::hts;
  using CppGrid::htw;
  using CppGrid::kmt;
  using CppGrid::tarea_r;
  using CppPconstMod::vit;
  // using CppIsopycMod::k1;
  // using CppIsopycMod::k2;
  // using CppIsopycMod::k3;
  using CppIsopycMod::dzr;
  using CppIsopycMod::dzw;
  using CppIsopycMod::dzwr;
  using CppIsopycMod::ahisop;
  using CppIsopycMod::adv_vbtiso;
  // using CppIsopycMod::adv_vetiso;
  using CppIsopycMod::adv_vntiso;
  // using CppWorkMod::tf;
  // using CppWorkMod::temp;
  // using CppWorkMod::work_1;
  // using CppWorkMod::work_2;
  // using CppWorkMod::work_3;
  using CppTracerMod::atb;
  using CppTracerMod::ax_iso;
  using CppTracerMod::ay_iso;
  using CppTracerMod::az_iso;
  using CppTracerMod::dx_iso;
  using CppTracerMod::dy_iso;
  using CppTracerMod::dz_iso;
  using CppTracerMod::ddy_iso;

	double work_1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
	double work_2[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
	double work_3[MAX_BLOCKS_CLINIC][KM+1][JMT][IMT];
	double temp[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  const int m = mtrace;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM-1; ++k) {
      for (int j = 1; j < JMT-2; ++j) {
        for (int i = 1; i < IMT-2; ++i) {
          temp[iblock][k][j][i] = P25 * dzr[k] 
              * ((atb[iblock][m][k  ][j+1][i] - atb[iblock][m][k+2][j+1][i])
               + (atb[iblock][m][k  ][j  ][i] - atb[iblock][m][k+2][j  ][i]));
        }
      }
    }
  }
  // k = 0
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-2; ++j) {
      for (int i = 1; i < IMT-2; ++i) {
        temp[iblock][0][j][i] = 0.25 * dzr[0] 
            * ((atb[iblock][m][1][j+1][i] - atb[iblock][m][2][j+1][i])
             + (atb[iblock][m][1][j  ][i] - atb[iblock][m][2][j  ][i]));
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        temp[iblock][KM-1][j][i] = C0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j+1][i]) - 1;
        if (k != -1) {
          const double fxe = dzw[k] + dzw[k+1];
          const double fxa = 0.5 
              * (atb[iblock][m][k  ][j+1][i] + atb[iblock][m][k  ][j][i]);
          const double fxb = 0.5 
              * (atb[iblock][m][k+1][j+1][i] + atb[iblock][m][k+1][j][i]);
          const double fxc = dzwr[k] * (fxb * fxe - fxa * dzw[k+1]);

          temp[iblock][k][j][i] = dzr[k] * (0.5 * (fxa + fxb) - fxc);
        }
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int i = 0; i < IMT; ++i) {
        work_1[iblock][k][0][i] = C0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_1[iblock][k][j][i] = hts[iblock][j][i] * (ahisop[iblock][j][i]
              * (atb[iblock][m][k+1][j+1][i] - atb[iblock][m][k+1][j][i])
                  / hue[iblock][j][i] + ahisop[iblock][j][i] * k2[iblock][0][j][k+1][i] 
                      * temp[iblock][k][j][i]) * vit[iblock][k][j][i] 
                          * vit[iblock][k][j+1][i]; 
        
          ddy_iso[iblock][m][k][j][i] = work_1[iblock][k][j][i];
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM-1; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          temp[iblock][k][j][i] = P25 * dzr[k]
              * (atb[iblock][m][k][j][i+1] - atb[iblock][m][k+2][j][i+1]
               + atb[iblock][m][k][j][i  ] - atb[iblock][m][k+2][j][i  ]);
        }
      }
    }
  }

  // k = 0
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        temp[iblock][0][j][i] = P25 * dzr[0]
            * (atb[iblock][m][1][j][i+1] - atb[iblock][m][2][j][i+1]
             + atb[iblock][m][1][j][i  ] - atb[iblock][m][2][j][i  ]);
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        temp[iblock][KM-1][j][i] = C0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j][i+1]) - 1;
        if (k != -1) {
          const double fxe = dzw[k] + dzw[k+1];
          const double fxa = 0.5 * 
              (atb[iblock][m][k  ][j][i] + atb[iblock][m][k  ][j][i+1]);
          const double fxb = 0.5 * 
              (atb[iblock][m][k+1][j][i] + atb[iblock][m][k+1][j][i+1]);
          const double fxc = dzwr[k] * (fxb * fxe - fxa * dzw[k+1]);

          temp[iblock][k][j][i] = dzr[k] * (P5 * (fxa + fxb) - fxc);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_2[iblock][k][j][i] = htw[iblock][j][i+1] * (ahisop[iblock][j][i]
              * (atb[iblock][m][k+1][j][i+1] - atb[iblock][m][k+1][j][i])
                  / hun[iblock][j][i+1] + ahisop[iblock][j][i] 
                      * k1[iblock][0][j][k+1][i] * temp[iblock][k][j][i])
                          * vit[iblock][k][j][i+1] * vit[iblock][k][j][i];
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_3[iblock][k][j][i] = ahisop[iblock][j][i] * P25 * vit[iblock][k][j][i] 
              * (k3[iblock][0][j][k][i] 
                  * (vit[iblock][k  ][j  ][i-1] 
                      * (atb[iblock][m][k+1][j  ][i  ] - atb[iblock][m][k+1][j  ][i-1]) 
                          / hun[iblock][j][i  ]
                   + vit[iblock][k-1][j][i-1]
                      * (atb[iblock][m][k  ][j  ][i  ] - atb[iblock][m][k  ][j  ][i-1]) 
                          / hun[iblock][j  ][i  ]
                   + vit[iblock][k  ][j  ][i+1]
                      * (atb[iblock][m][k+1][j  ][i+1] - atb[iblock][m][k+1][j  ][i  ]) 
                          / hun[iblock][j  ][i+1]
                   + vit[iblock][k-1][j  ][i+1]
                      * (atb[iblock][m][k  ][j  ][i+1] - atb[iblock][m][k  ][j  ][i  ]) 
                          / hun[iblock][j  ][i+1])
              + k3[iblock][1][j][k][i]
                  * (vit[iblock][k  ][j-1][i  ]
                      * (atb[iblock][m][k+1][j  ][i  ] - atb[iblock][m][k+1][j-1][i  ])
                          / hue[iblock][j-1][i  ]
                   + vit[iblock][k-1][j-1][i  ]
                      * (atb[iblock][m][k  ][j  ][i  ] - atb[iblock][m][k  ][j-1][i  ])
                          / hue[iblock][j-1][i  ]
                   + vit[iblock][k  ][j+1][i  ]
                      * (atb[iblock][m][k+1][j+1][i  ] - atb[iblock][m][k+1][j  ][i  ])
                          / hue[iblock][j  ][i  ]
                   + vit[iblock][k-1][j+1][i  ]
                      * (atb[iblock][m][k  ][j+1][i  ] - atb[iblock][m][k  ][j  ][i  ])
                          / hue[iblock][j  ][i  ]));
        }
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work_3[iblock][0][j][i]  = C0;
        work_3[iblock][KM][j][i] = C0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          tf[iblock][k][j][i] = tf[iblock][k][j][i] 
              + tarea_r[iblock][j][i]
                  * (work_1[iblock][k][j][i] - work_1[iblock][k  ][j-1][i  ])
              + tarea_r[iblock][j][i]
                  * (work_2[iblock][k][j][i] - work_2[iblock][k  ][j  ][i-1])
              + dzr[k]
                  * (work_3[iblock][k][j][i] - work_3[iblock][k+1][j  ][i  ]);

          dx_iso[iblock][m][k][j][i] = tarea_r[iblock][j][i]
              * (work_2[iblock][k][j][i] - work_2[iblock][k  ][j  ][i-1]);
          dy_iso[iblock][m][k][j][i] = tarea_r[iblock][j][i]
              * (work_1[iblock][k][j][i] - work_1[iblock][k  ][j-1][i  ]);
          dz_iso[iblock][m][k][j][i] = dzr[k]
              * (work_3[iblock][k][j][i] - work_3[iblock][k+1][j  ][i  ]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_1[iblock][k][j][i] = adv_vntiso[iblock][j][k][i]
              * (atb[iblock][m][k+1][j+1][i] + atb[iblock][m][k+1][j][i])
                  * hts[iblock][j][i];
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_2[iblock][k][j][i] = adv_vetiso[iblock][j][k][i]
              * (atb[iblock][m][k+1][j][i+1] + atb[iblock][m][k+1][j][i])
                  * htw[iblock][j][i+1];
        }
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work_3[iblock][0][j][i]  = C0;
        work_3[iblock][KM][j][i] = C0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          work_3[iblock][k][j][i] = adv_vbtiso[iblock][j][k][i]
              * (atb[iblock][m][k+1][j][i] + atb[iblock][m][k][j][i]);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          tf[iblock][k][j][i] = tf[iblock][k][j][i]
              - P5 * (work_1[iblock][k][j  ][i  ] - work_1[iblock][k][j-1][i  ]
                    + work_2[iblock][k][j  ][i  ] - work_2[iblock][k][j  ][i-1])
                        * tarea_r[iblock][j][i]
              - P5 * dzr[k] 
                   * (work_3[iblock][k][j][i] - work_3[iblock][k+1][j][i]);

          ax_iso[iblock][m][k][j][i] = - P5
              * (work_2[iblock][k  ][j  ][i  ] - work_2[iblock][k  ][j  ][i-1])
                  * tarea_r[iblock][j][i];
          ay_iso[iblock][m][k][j][i] = - P5
              * (work_1[iblock][k  ][j  ][i  ] - work_1[iblock][k  ][j-1][i  ])
                  * tarea_r[iblock][j][i];
          az_iso[iblock][m][k][j][i] = - P5 * dzr[k]
              * (work_3[iblock][k  ][j  ][i  ] - work_3[iblock][k+1][j  ][i  ]);
        }
      }
    }
  }

  return ;
}

#endif // ISO

#endif // LICOM_ENABLE_FORTRAN