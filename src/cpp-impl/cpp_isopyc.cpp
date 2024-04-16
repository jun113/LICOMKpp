#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#ifdef ISO

#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_isopyc_mod.h"
#include "../head/cpp_tracer_mod.h"

#include "../head/cpp_extern_functions.h"

#include <cmath>

void isopyc(double (&k1)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  					double (&k2)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  					double (&k3)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT],
						double (&adv_vetiso)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]) {

	using CppDomain::nblocks_clinic;
  using CppPconstMod::to;
  using CppPconstMod::so;
  // using CppIsopycMod::rhoi;
  using CppIsopycMod::krplin;
  using CppTracerMod::atb;

	double e[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT];
	double rhoi[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int m = 0; m < NRPL; ++m) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          rhoi[iblock][m][j][0][i] = 0.0;
        }
      }
    }
  }
   
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int m = 0; m < NRPL; ++m) {
      double tref0 = to[krplin[m] - 1];
      double sref0 = so[krplin[m] - 1];
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            double tq = atb[iblock][0][k+1][j][i] - tref0;
            double sq = atb[iblock][1][k+1][j][i] - sref0;
            int    kq = krplin[m];
            rhoi[iblock][m][j][k+1][i] = dens(tq, sq, kq-1);
          }
        }
      }
    }
  }

  k2_3 (rhoi, e, k2);
  k1_3 (rhoi, e, k1);
  k3_123 (rhoi, e, k3);
  isoadv (k1, k2, adv_vetiso);

  return ;
}

void k2_3 (const double (&rhoi)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&e)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&k2)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT]) {

	using CppDomain::nblocks_clinic;
  using CppGrid::dxu;
  using CppGrid::hue;
  using CppGrid::kmt;
  // using CppIsopycMod::e;
  // using CppIsopycMod::k2;
  using CppIsopycMod::dzr;
  using CppIsopycMod::dzw;
  using CppIsopycMod::dzwr;
  // using CppIsopycMod::rhoi;
  using CppIsopycMod::tmask;
  using CppIsopycMod::slmxr;
  using CppIsopycMod::fzisop;
  using CppIsopycMod::kisrpl;
  int m;
  double fxa, fxb, fxc, fxd, fxe;

  const double c1e10 = 1.0e10;
  const double eps   = 1.0e-25;
  const double p5    = 0.5;
  const double c0    = 0.0;
  const double c1    = 1.0;
  const double p25   = 0.25;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 1; k < KM-1; ++k) {
        m = static_cast<int>(kisrpl[k]) - 1;
        fxd = c1e10 * p25 * dzr[k];

        for (int i = 1; i < IMT-1; ++i) {
          e[iblock][2][j][k][i] = fxd 
              * (rhoi[iblock][m][j  ][k][i] - rhoi[iblock][m][j  ][k+2][i]
               + rhoi[iblock][m][j+1][k][i] - rhoi[iblock][m][j+1][k+2][i]);
        }
      }
    }
  }
  fxd = c1e10 * dzr[0];
  fxe = dzw[0] + dzw[1];
  m = static_cast<int>(kisrpl[0]) - 1;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        fxa = p5 * (rhoi[iblock][m][j][2][i] + rhoi[iblock][m][j+1][2][i]);
        fxb = p5 * (rhoi[iblock][m][j][1][i] + rhoi[iblock][m][j+1][1][i]);
        fxc = dzwr[1] * (fxb * fxe - fxa * dzw[0]);
        e[iblock][2][j][0][i] = fxd * (fxc - p5 * (fxa + fxb));
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        e[iblock][2][j][KM-1][i] = c0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j+1][i]);
        if (k != 0) {
          fxe = dzw[k-1] + dzw[k];
          m = static_cast<int>(kisrpl[k-1]) - 1;
          fxa = p5 * (rhoi[iblock][m][j][k-1][i] + rhoi[iblock][m][j+1][k-1][i]);
          fxb = p5 * (rhoi[iblock][m][j][k  ][i] + rhoi[iblock][m][j+1][k  ][i]);
          fxc = dzwr[k-1] * (fxb * fxe - fxa * dzw[k]);
          e[iblock][2][j][k-1][i] = dzr[k-1] * c1e10 * (p5 * (fxa + fxb) - fxc);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        m = static_cast<int>(kisrpl[k]) - 1;
        for (int i = 1; i < IMT-1; ++i) {
          e[iblock][0][j][k][i] = p5 * c1e10
              / (dxu[iblock][j][i] + dxu[iblock][j][i+1])
              * (rhoi[iblock][m][j+1][k+1][i+1] - rhoi[iblock][m][j+1][k+1][i-1]
               + rhoi[iblock][m][j  ][k+1][i+1] - rhoi[iblock][m][j  ][k+1][i-1]);

          e[iblock][1][j][k][i] = tmask[iblock][j][k][i] 
              * tmask[iblock][j+1][k][i] / hue[iblock][j][i] * c1e10 
                  * (rhoi[iblock][m][j+1][k+1][i] - rhoi[iblock][m][j][k+1][i]); 

        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double olmask = 
              tmask[iblock][j+1][k][i-1] * tmask[iblock][j+1][k][i+1]
            * tmask[iblock][j  ][k][i-1] * tmask[iblock][j  ][k][i+1];
          if (olmask < c1) {
            e[iblock][0][j][k][i] = c0;
          }
        }
      }
    }
  }

#ifdef LDD97

  double f1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double f2[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          f1[iblock][k][j][i] = 0.0;
          f2[iblock][k][j][i] = 0.0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double chkslp = - sqrt(pow(e[iblock][0][j][k][i], 2) 
              + pow(e[iblock][1][j][k][i], 2)) * slmxr; 
          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          const double slopemod = sqrt(pow(e[iblock][0][j][k][i], 2)
              + pow(e[iblock][1][j][k][i], 2))
                  / fabs(e[iblock][2][j][k][i] + eps);
          f1[iblock][k][j][i] = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));
          const double nondimr = - zkt[k] / (rrd2[iblock][j][i]
              * (slopemod + eps));
          if (nondimr >= 1.0) {
            f2[iblock][k][j][i] = 1.0;
          } else {
            f2[iblock][k][j][i] = 0.5 * (1.0 + sin(PI * (nondimr - 0.5));
          }
          k2[iblock][0][j][k+1][i] = (- e[iblock][1][j][k][i] 
              * e[iblock][2][j][k][i] * f1[iblock][k][j][i] * f2[iblock][k][j][i])
                  / (pow(e[iblock][2][j][k][i], 2) + eps);
        }
      }
    }
  }

#else  // LDD97
  
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double chkslp = - sqrt(pow(e[iblock][0][j][k][i], 2) 
              + pow(e[iblock][1][j][k][i], 2)) * slmxr;
          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          // Note: allocate(k2(imt, 0:km, jmt, 3:3, max_blocks_clinic))
          k2[iblock][0][j][k+1][i] = (- e[iblock][1][j][k][i]
              * e[iblock][2][j][k][i] * fzisop[k])
                  / (pow(e[iblock][2][j][k][i], 2) + eps);
        }
      }
    }
  }

#endif // LDD97

  return ;
}

void k1_3 (const double (&rhoi)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&e)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&k1)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT]) {

	using CppDomain::nblocks_clinic;
  using CppGrid::dyu;
  using CppGrid::hun;
  using CppGrid::kmt;
  // using CppIsopycMod::e;
  // using CppIsopycMod::k1;
  using CppIsopycMod::dzr;
  using CppIsopycMod::dzw;
  using CppIsopycMod::dzwr;
  // using CppIsopycMod::rhoi;
  using CppIsopycMod::tmask;
  using CppIsopycMod::slmxr;
  using CppIsopycMod::fzisop;
  using CppIsopycMod::kisrpl;

  const double c1    = 1.0;
  const double c0    = 0.0;
  const double p5    = 0.5;
  const double p25   = 0.25;
  const double eps   = 1.0e-25;
  const double c1e10 = 1.0e10;

  int m;
  double fxa, fxb, fxc, fxd, fxe;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 1; k < KM-1; ++k) {
        m = static_cast<int>(kisrpl[k]) - 1; 
        fxd = c1e10 * p25 * dzr[k]; 
        for (int i = 0; i < IMT-1; ++i) {
          e[iblock][2][j][k][i] = fxd 
              * (rhoi[iblock][m][j][k][i  ] - rhoi[iblock][m][j][k+2][i  ]
               + rhoi[iblock][m][j][k][i+1] - rhoi[iblock][m][j][k+2][i+1]);
        }
      }
    }
  }

  fxd = c1e10 * dzr[0];
  fxe = dzw[0] + dzw[1];
  m = static_cast<int>(kisrpl[0]) - 1;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 0; i < IMT-1; ++i) {
        fxa = p5 * (rhoi[iblock][m][j][2][i] + rhoi[iblock][m][j][2][i+1]);
        fxb = p5 * (rhoi[iblock][m][j][1][i] + rhoi[iblock][m][j][1][i+1]);
        fxc = dzwr[1] * (fxb * fxe - fxa * dzw[0]);
        e[iblock][2][j][0][i] = fxd * (fxc - p5 * (fxa + fxb));
      }
    }
  }

  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        e[iblock][2][j][KM-1][i] = c0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 0; i < IMT-1; ++i) {
        const int k = std::min(kmt[iblock][j][i], kmt[iblock][j][i+1]) - 1;
        if (k != -1) {
          fxe = dzw[k] + dzw[k+1];
          m = static_cast<int>(kisrpl[k]) - 1;
          fxa = p5 * (rhoi[iblock][m][j][k  ][i] + rhoi[iblock][m][j][k  ][i+1]);
          fxb = p5 * (rhoi[iblock][m][j][k+1][i] + rhoi[iblock][m][j][k+1][i+1]);
          fxc = dzwr[k] * (fxb * fxe - fxa * dzw[k+1]);
          e[iblock][2][j][k][i] = dzr[k] * c1e10 * (p5 * (fxa + fxb) - fxc);
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        m = static_cast<int>(kisrpl[k]) - 1;
        for (int i = 1; i < IMT-1; ++i) {
          e[iblock][0][j][k][i] = tmask[iblock][j][k][i] 
              * tmask[iblock][j][k][i+1] / hun[iblock][j][i+1] * c1e10
                  * (rhoi[iblock][m][j][k+1][i+1] - rhoi[iblock][m][j][k+1][i]);

          e[iblock][1][j][k][i] = p5 * c1e10 
              / (dyu[iblock][j-1][i+1] + dyu[iblock][j][i+1])
              * (rhoi[iblock][m][j+1][k+1][i  ] - rhoi[iblock][m][j-1][k+1][i  ]
               + rhoi[iblock][m][j+1][k+1][i+1] - rhoi[iblock][m][j-1][k+1][i+1]);
      
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < IMT-1; ++i) {
          const double olmask = tmask[iblock][j-1][k][i] 
              * tmask[iblock][j+1][k][i] * tmask[iblock][j-1][k][i+1]
                  * tmask[iblock][j+1][k][i+1];
          if (olmask < c1) {
            e[iblock][1][j][k][i] = c0;
          }
        }
      }
    }
  }
#ifdef LDD97

  double f1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double f2[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          f1[iblock][k][j][i] = 1.0;
          f2[iblock][k][j][i] = 1.0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < IMT-1; ++i) {
          const double chkslp = - sqrt(pow(e[iblock][0][j][k][i], 2)
              + pow(e[iblock][1][j][k][i], 2)) * slmxr;
          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          const double slopemod = sqrt(pow(e[iblock][0][j][k][i], 2)
              + pow(e[iblock][1][j][k][i], 2))
                  / fabs(e[iblock][2][j][k][i] + eps);
          f1[iblock][k][j][i] = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));
          const double nondimr = - zkt[k] / (rrd1[iblock] * (slopemod + eps));
          if (nondimr >= 1.0) {
            f2[iblock][k][j][i] = 1.0;
          } else {
            f2[iblock][k][j][i] = 0.5 * (1.0 + sin(PI * (nondimr - 0.5)));
          }
          k1[iblock][0][j][k+1][i] = ( - e[iblock][0][j][k][i]
              * e[iblock][2][j][k][i] * f1[iblock][k][j][i] * f2[iblock][k][j][i]) 
                  / (pow(e[iblock][2][j][k][i], 2) + eps);
        }
      }
    }
  }
#else  // LDD97

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < IMT-1; ++i) {
          const double chkslp = - sqrt(pow(e[iblock][0][j][k][i], 2)
              + pow(e[iblock][1][j][k][i], 2)) * slmxr;
          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          // Note: allocate(k1(imt, 0:km, jmt, 3:3, max_blocks_clinic))
          k1[iblock][0][j][k+1][i] = ( - e[iblock][0][j][k][i]
              * e[iblock][2][j][k][i] * fzisop[k])
                  / (pow(e[iblock][2][j][k][i], 2) + eps);
        }
      }
    }
  }

#endif // LDD97

  return ;
}

void k3_123 (const double (&rhoi)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&e)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&k3)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT]) {

	using CppDomain::nblocks_clinic;
  using CppGrid::hue;
  using CppGrid::hun;
  // using CppIsopycMod::e;
  // using CppIsopycMod::k3;
  using CppIsopycMod::dzwr;
  // using CppIsopycMod::rhoi;
  using CppIsopycMod::tmask;
  using CppIsopycMod::slmxr;
  using CppIsopycMod::kisrpl;

  const double c0    = 0.0;
  const double p25   = 0.25;
  const double eps   = 1.0e-25;
  const double c1e10 = 1.0e10;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 1; k < KM; ++k) {
        const int m = static_cast<int>(kisrpl[k]) - 1;
        for (int i = 1; i < IMT-1; ++i) {
          e[iblock][0][j][k-1][i] = p25 * c1e10
              * (tmask[iblock][j][k-1][i-1] * tmask[iblock][j][k-1][i  ]
                  * (rhoi[iblock][m][j][k  ][i  ] - rhoi[iblock][m][j][k  ][i-1]) 
                      / hun[iblock][j][i  ]
               + tmask[iblock][j][k-1][i  ] * tmask[iblock][j][k-1][i+1]
                  * (rhoi[iblock][m][j][k  ][i+1] - rhoi[iblock][m][j][k  ][i  ]) 
                      / hun[iblock][j][i+1]
               + tmask[iblock][j][k  ][i-1] * tmask[iblock][j][k  ][i  ]
                  * (rhoi[iblock][m][j][k+1][i  ] - rhoi[iblock][m][j][k+1][i-1]) 
                      / hun[iblock][j][i  ]
               + tmask[iblock][j][k  ][i  ] * tmask[iblock][j][k  ][i+1]
                  * (rhoi[iblock][m][j][k+1][i+1] - rhoi[iblock][m][j][k+1][i  ]) 
                      / hun[iblock][j][i+1]);

          e[iblock][1][j][k-1][i] = p25 * c1e10
              * (tmask[iblock][j-1][k-1][i] * tmask[iblock][j  ][k-1][i]
                  * (rhoi[iblock][m][j  ][k  ][i] - rhoi[iblock][m][j-1][k  ][i]) 
                      / hue[iblock][j-1][i]
               + tmask[iblock][j  ][k-1][i] * tmask[iblock][j+1][k-1][i]
                  * (rhoi[iblock][m][j+1][k  ][i] - rhoi[iblock][m][j  ][k  ][i]) 
                      / hue[iblock][j  ][i]
               + tmask[iblock][j-1][k  ][i] * tmask[iblock][j][k  ][i  ]
                  * (rhoi[iblock][m][j  ][k+1][i] - rhoi[iblock][m][j-1][k+1][i]) 
                      / hue[iblock][j-1][i]
               + tmask[iblock][j  ][k  ][i] * tmask[iblock][j+1][k  ][i]
                  * (rhoi[iblock][m][j+1][k+1][i] - rhoi[iblock][m][j  ][k+1][i]) 
                      / hue[iblock][j  ][i]);
          
          e[iblock][2][j][k-1][i] = dzwr[k] 
              * tmask[iblock][j][k-1][i] * tmask[iblock][j][k][i] * c1e10
                  * (rhoi[iblock][m][j][k][i] - rhoi[iblock][m][j][k+1][i]);
        }
      }
      for (int i = 0; i < IMT; ++i) {
        e[iblock][0][j][KM-1][i] = c0;
        e[iblock][1][j][KM-1][i] = c0;
        e[iblock][2][j][KM-1][i] = c0;
      }
    }
  }
#ifdef LDD97

  double f1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double f2[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          f1[iblock][k][j][i] = 0.0;
          f2[iblock][k][j][i] = 0.0;
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 0; i < IMT-1; ++i) {
          const double chkslp = - sqrt (
              pow(e[iblock][0][j][k][i], 2) + pow(e[iblock][1][j][k][i], 2))
                  * slmxr;

          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          const double slopemod = sqrt(
              pow(e[iblock][0][j][k][i], 2) + pow(e[iblock][1][j][k][i], 2))
                  / fabs(e[iblock][2][j][k][i] + eps);

          f1[iblock][k][j][i] = 0.5 * (1.0 + tanh((0.004 - slopemod) / 0.001));

          const double nondimr = - zkp[k] / (rrd1[iblock][j][i] 
              * (slopemod + eps));

          if (nondimr >= 1.0) {
            f2[iblock][k][j][i] = 1.0;
          } else {
            f2[iblock][k][j][i] = 0.5 * (1.0 + sin(PI * (nondimr - 0.5)));
          }
          const double ahfctr = 1.0 / (pow(e[iblock][2][j][k][i], 2) + eps)
              * f1[iblock][k][j][i] * f2[iblock][k][j][i];

          k3[iblock][0][j][k+1][i] = - e[iblock][2][j][k][i] 
                                     * e[iblock][0][j][k][i] * ahfctr;
               
          k3[iblock][1][j][k+1][i] = - e[iblock][2][j][k][i] 
                                     * e[iblock][1][j][k][i] * ahfctr;

          k3[iblock][2][j][k+1][i] = (pow(e[iblock][0][j][k][i], 2) 
                                    + pow(e[iblock][1][j][k][i], 2));
        }
      }
    }
  }

#else //  LDD97

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 1; j < JMT-1; ++j) {
      for (int k = 0; k < KM; ++k) {
        for (int i = 1; i < IMT-1; ++i) {
          const double chkslp = - sqrt (
              pow(e[iblock][0][j][k][i], 2) + pow(e[iblock][1][j][k][i], 2))
                  * slmxr;

          if (e[iblock][2][j][k][i] > chkslp) {
            e[iblock][2][j][k][i] = chkslp;
          }
          const double ahfctr = 1.0 / (pow(e[iblock][2][j][k][i], 2) + eps);

          // Note: allocate(k3(imt, 0:km, jmt, 1:3, max_blocks_clinic))
          k3[iblock][0][j][k+1][i] = - e[iblock][2][j][k][i]
                                     * e[iblock][0][j][k][i] * ahfctr;
                                    
          k3[iblock][1][j][k+1][i] = - e[iblock][2][j][k][i]
                                     * e[iblock][1][j][k][i] * ahfctr;

          k3[iblock][2][j][k+1][i] = (pow(e[iblock][0][j][k][i], 2)
                                    + pow(e[iblock][1][j][k][i], 2)) * ahfctr;
        }
      }
    }
  }

#endif // LDD97
  return ;
}
#endif // ISO

#endif // LICOM_ENABLE_FORTRAN