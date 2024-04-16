#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_blocks.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pmix_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_work_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

using CppDomain::nblocks_clinic;
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMM1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;
using CppParamMod::JST;
using CppParamMod::JET;

void cpp_readyt() {

  using CppDynMod::dlu;
  using CppDynMod::dlv;
  using CppDynMod::gg;
  using CppDynMod::h0;
  using CppDynMod::h0l;
  using CppDynMod::h0f;
  using CppDynMod::u;
  using CppDynMod::v;
  using CppDynMod::utl;
  using CppDynMod::vtl;
  using CppDynMod::utf;
  using CppDynMod::vtf;

  using CppForcMod::psa;
  using CppForcMod::swv;
  using CppForcMod::nswv;
  using CppForcMod::buoytur;
  using CppForcMod::buoysol;

  using CppPconstMod::to;
  using CppPconstMod::so;
  using CppPconstMod::po;
  using CppPconstMod::akt;
  using CppPconstMod::dzp;
  using CppPconstMod::ist;
  using CppPconstMod::viv;
  using CppPconstMod::vit;
  using CppPconstMod::od0;
  using CppPconstMod::od0cp;
  using CppPconstMod::odzt;
  using CppPconstMod::ohbt;
  using CppPconstMod::zkt;
  using CppPconstMod::hbx;
  using CppPconstMod::hby;

  using CppPmixMod::ric;
  using CppPmixMod::rit;
  using CppPmixMod::rict;
  using CppPmixMod::rict_replace;
  using CppPmixMod::ricdt;
  using CppPmixMod::ricdttms;

  using CppWorkMod::pax;
  using CppWorkMod::pay;
  using CppWorkMod::pxb;
  using CppWorkMod::pyb;
  using CppWorkMod::whx;
  using CppWorkMod::why;
  using CppWorkMod::wgp;
  using CppWorkMod::work;

  using CppTracerMod::ax;
  using CppTracerMod::ay;
  using CppTracerMod::az;
  using CppTracerMod::at;
  using CppTracerMod::atb;
  using CppTracerMod::pdensity;

  using CppConstantMod::C0;
  using CppConstantMod::G;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int n = 0; n < NTRA; ++n) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            akt[iblock][n][k][j][i] = 0.0;
          }
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          gg[iblock][k][j][i]  = 0.0;
          dlu[iblock][k][j][i] = 0.0;
          dlv[iblock][k][j][i] = 0.0;
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          rit[iblock][k][j][i]      = 0.0;
          ric[iblock][k][j][i]      = 0.0;
          rict[iblock][k][j][i]     = 0.0;
          ricdt[iblock][k][j][i]    = 0.0;
          ricdttms[iblock][k][j][i] = 0.0;
        }
      }
    }
  }
  
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = JST - 1; j < JET; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h0l[iblock][j][i] = h0f[iblock][j][i];
        h0f[iblock][j][i] = h0 [iblock][j][i];
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = JST - 1; j < JET; ++j) {
        for (int i = 0; i < IMT; ++i) {
          utl[iblock][k][j][i] = utf[iblock][k][j][i];
          vtl[iblock][k][j][i] = vtf[iblock][k][j][i];
          utf[iblock][k][j][i] = u[iblock][k][j][i];
          vtf[iblock][k][j][i] = v[iblock][k][j][i];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      // NOTE: ATB[0:KM]
      tgrid_to_ugrid (dlu[iblock][k], atb[iblock][0][k+1], iblock);
      tgrid_to_ugrid (dlv[iblock][k], atb[iblock][1][k+1], iblock);
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          double tup = dlu[iblock][k][j][i]   - to[k+1];
          double sup = dlv[iblock][k][j][i]   - so[k+1];
          double tlo = dlu[iblock][k+1][j][i] - to[k+1];
          double slo = dlv[iblock][k+1][j][i] - so[k+1];
          double rhoup = dens(tup, sup, k+1);
          double rholo = dens(tlo, slo, k+1);
          ric[iblock][k][j][i] = viv[iblock][k+1][j][i] * 
              od0 * G * (rholo - rhoup) * odzt[k+1];
        }
      }
    }
  }
  density ();
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          // NOTE: ATB[0:KM]
          double tup = atb[iblock][0][k+1][j][i] - to[k+1];
          double sup = atb[iblock][1][k+1][j][i] - so[k+1];
          double tlo = atb[iblock][0][k+2][j][i] - to[k+1];
          double slo = atb[iblock][1][k+2][j][i] - so[k+1];

          double rhoup = dens(tup, sup, k+1);
          double rholo = dens(tlo, slo, k+1);

          rict[iblock][k][j][i] = vit[iblock][k+1][j][i] * 
              od0 * G * (rholo - rhoup) * odzt[k+1];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          rict_replace[iblock][k][j][i] = rict[iblock][k][j][i];
        }
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          gg[iblock][k][j][i] = - od0 * G * 
              pdensity[iblock][k][j][i] *
                  vit[iblock][k][j][i];
        }
      }
    }
  }
  double pp[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  double ppa[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double ppb[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double ppc[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        pp[iblock][0][j][i] = gg[iblock][0][j][i] * 0.5 *
            dzp[0] * vit[iblock][0][j][i];
        ppa[iblock][0][j][i] = psa[iblock][j][i] * 
            vit[iblock][0][j][i];
        ppb[iblock][0][j][i] = at[iblock][0][0][j][i] * 
            vit[iblock][0][j][i];
        ppc[iblock][0][j][i] = at[iblock][1][0][j][i] * 
            vit[iblock][0][j][i];
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 1; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          pp[iblock][k][j][i] = vit[iblock][k][j][i] * 
              (pp[iblock][k-1][j][i] + 0.5 * (gg[iblock][k][j][i] * dzp[k] + 
               gg[iblock][k-1][j][i] * dzp[k-1]));
          ppa[iblock][k][j][i] = vit[iblock][k][j][i] * 
              (ppa[iblock][k-1][j][i] + gg[iblock][k-1][j][i] *
              dzp[k-1]);

          ppb[iblock][k][j][i] = vit[iblock][k][j][i] * 
              (at[iblock][0][k-1][j][i] - 
              (at[iblock][0][k-1][j][i] - at[iblock][0][k][j][i]) * 
              dzp[k-1] / (dzp[k-1] + dzp[k]));

          ppc[iblock][k][j][i] = vit[iblock][k][j][i] * 
              (at[iblock][1][k-1][j][i] - 
              (at[iblock][1][k-1][j][i] - at[iblock][1][k][j][i]) * 
              dzp[k-1] / (dzp[k-1] + dzp[k]));
        }
      }
    }
  }

  double alpha[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  double beta[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  thermal (ppb, ppc, ppa, alpha, beta, vit);

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        buoytur[iblock][j][i] = vit[iblock][0][j][i] *
            nswv[iblock][j][i] * G * alpha[iblock][0][j][i] * od0cp;
        buoysol[iblock][j][i] = vit[iblock][0][j][i] *
            swv[iblock][j][i] * G * alpha[iblock][0][j][i] * od0cp;
      }
    }
  }

  //    epsln = 1.0D-25 
  const double epsln = 1.0e-25;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KMM1; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ricdttms[iblock][k][j][i] = vit[iblock][k+1][j][i] * G *
              ((at[iblock][0][k][j][i] - at[iblock][0][k+1][j][i]) * alpha[iblock][k+1][j][i] + 1000.0 * 
               (at[iblock][1][k][j][i] - at[iblock][1][k+1][j][i]) *
                beta[iblock][k+1][j][i]) * odzt[k+1];
          ricdt[iblock][k][j][i] = vit[iblock][k+1][j][i] / 
              ((at[iblock][0][k][j][i] - at[iblock][0][k+1][j][i] + epsln) * 
                alpha[iblock][k+1][j][i]) * 1000.0 *
              ((at[iblock][1][k][j][i] - at[iblock][1][k+1][j][i]) * 
                beta[iblock][k+1][j][i]);
#ifdef CANUTOMIXOUT
          CppPconstMod::alpha_canuto[iblock][k][j][i] = alpha[iblock][k][j][i];
          CppPconstMod::beta_canuto[iblock][k][j][i]  = beta[iblock][k][j][i];
#endif // CANUTOMIXOUT
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          gg[iblock][k][j][i] = - od0 * G *
              (pdensity[iblock][k][j][i] - po[k] - 1000.0) *
               vit[iblock][k][j][i];
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      grad (k, dlu[iblock][k], dlv[iblock][k], pp[iblock][k]);
    }
  }

  vinteg(dlu, pxb);
  vinteg(dlv, pyb);

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        dlu[iblock][0][j][i] = 0.0;
        dlu[iblock][1][j][i] = 0.0;
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          double abcd = gg[iblock][k][j][i] * 
              ohbt[iblock][j][i] * dzp[k];
          dlu[iblock][0][j][i] = dlu[iblock][0][j][i] + 
              abcd; 
          dlu[iblock][1][j][i] = dlu[iblock][1][j][i] + 
              abcd * zkt[k]; 
        }
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        dlv[iblock][0][j][i] = (dlu[iblock][0][j][i] +
            dlu[iblock][1][j][i] * ohbt[iblock][j][i]) / G;
        dlv[iblock][1][j][i] = dlu[iblock][1][j][i] *
            ohbt[iblock][j][i] * ohbt[iblock][j][i];
      }
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    tgrid_to_ugrid (wgp[iblock], dlv[iblock][0], iblock);
    tgrid_to_ugrid (work[iblock], dlv[iblock][1], iblock);
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        whx[iblock][j][i] = hbx[iblock][j][i] * 
            work[iblock][j][i] * viv[iblock][0][j][i];
        why[iblock][j][i] = hby[iblock][j][i] * 
            work[iblock][j][i] * viv[iblock][0][j][i];
      }
    }
  }
  double work1[MAX_BLOCKS_CLINIC][JMT][IMT];
  double work2[MAX_BLOCKS_CLINIC][JMT][IMT];
  for (int iblock = 0; iblock < MAX_BLOCKS_CLINIC; ++iblock) {
    grad (0, work1[iblock], work2[iblock], psa[iblock]);
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        pay[iblock][j][i] = 
            - od0 * work2[iblock][j][i];
        pax[iblock][j][i] = 
            - od0 * work1[iblock][j][i];
      }
    }
  }

  fortran_mpi_barrier_();

  if (ist == 0) {
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int n = 0; n < NTRA; ++n) {
        for (int k = 0; k < KM; ++k) {
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              ax[iblock][n][k][j][i] = C0;
              ay[iblock][n][k][j][i] = C0;
              az[iblock][n][k][j][i] = C0;
            }
          }
        }
      }
    }
  }
  return ;
}

#endif // LICOM_ENABLE_FORTRAN
