#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_canuto_mod.h"
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_param_mod.h"

#include "../head/fortran_canuto_mod.h"

#include <cstdio>
#include <cmath>
#include <vector>

using CppParamMod::KM;

static void formld(const double *, const double *,
    double &, const int &);

static void formld_te(const double *, const double *, 
    const double &, double &, const int &);

static void formld_rh(const double *, const double *, 
    const double &, double &, const int &);

static void interp1d_expabs(double &, const double *, 
    const double *,  const double *, const double *,
        const double *, double &, double &, double &, double &,
            const int &, const int &, const int &,
                const double &, const double &);

static void interp2d_expabs(double &, double &,
    double &, double &, double &, double &);

static double eplatidepend_(const double &, const double &);

void turb_2(
    const double (&z)[KM],        // wp8
    const double (&t)[KM],        // wp1
    const double (&s)[KM],        // wp2
    const double (&rh)[KM],       // wp3
          double (&ri)[KM],       // wp4
    const double (&rid)[KM],      // wp5
          double (&s2)[KM],       // wp6
    const double &fricmx,         // DFRICMX * 1.0e+4
    const double &wndmix,         // DWNDMIX * 1.0e+4
          double (&v_back)[KM-1], // akm_back
          double (&t_back)[KM-1], // akt_back
          double (&s_back)[KM-1], // aks_back
    const double (&an2)[KM],      // wp7
    const double &ustar_,         // wp9
    const double &buoytur,        // wp10
    const double &buoysol,        // wp11
    const double &coriol,         // fcort[iblock][j][i]
          double &amld,           // mldtmp
          double (&akm)[KM-1],    // wk1
          double (&akh)[KM-1],    // wk2
          double (&aks)[KM-1],    // wk3
    const int    &n,              // iwk
    const int    &na,             // iwk1
    const int    &nmax,           // km - 1
    const int    &isurfuse,       // 1
    const int    &ifextermld,     // 0
    const int    &ifoutput,       // 0
    const int    &ii,             // i
    const int    &jj              /* j */) {

  using CppCanutoMod::b1;
  using CppCanutoMod::dri;
  using CppCanutoMod::back_l_0;
  using CppCanutoMod::back_ra_r;
  using CppCanutoMod::slq2a;
  using CppCanutoMod::slq2_r1;
  using CppCanutoMod::sha;
  using CppCanutoMod::sma;
  using CppCanutoMod::sma1;
  using CppCanutoMod::sha1;
  using CppCanutoMod::ssa1;
  using CppCanutoMod::sh_r1;
  using CppCanutoMod::ss_r1;
  using CppCanutoMod::sm_r1;
  using CppCanutoMod::ria;
  using CppCanutoMod::rib;
  using CppCanutoMod::rnd2on2;
  using CppCanutoMod::sisamax;
  using CppCanutoMod::dand2on2;
  using CppCanutoMod::and2on2a1;
  using CppCanutoMod::amtaun2a1;
  using CppCanutoMod::deltheta_r;
  using CppCanutoMod::theta_rcrp;
  using CppCanutoMod::theta_rcrn;
  using CppCanutoMod::eps_bot__under;


  using CppParamMod::mytid;

  using std::vector;

  double an = 0.0;
  double al_back = 0.0;
  double sh_back = 0.0;
  double sm_back = 0.0;
  double ss_back = 0.0;
  double slq2_back = 0.0;

  double ri1   = 0.0;
  double rid1  = 0.0;
  double rimax = 0.0;


  int ifbelow;
  int ifpureshear = 0;

  int ifnofsmall    = 0;
  //int ilmldpoint    = 0;
  //int ilmldpointneg = 0;

  if (nmax > NBIG) {
    if (mytid == 0) {
      printf("\n\n");
      printf("****************************\n");
      printf("PROBLEM IN TURB_2 ROUTINE.\n");
      printf("Number of model levels exceeds nbig.\n");
      printf("nbig = %d,\tnmax = %d\n", NBIG, nmax);
      printf("If want to use this many levels increase nbig to be bigger than nmax.\n");
      printf("Program is stopping now.\n");
      exit(0);
    }
  }
  int nb = std::max(0, n);
  int kbot = n + 1;
  double visc_cbu_limit1 = fricmx;
  double diff_cbt_limit1 = fricmx;
  double buoytot = buoytur + buoysol;

  if (ifextermld == 0) {
    if (n > 0) {
      if (IDEFMLD == 0) {
        formld(z, t, amld, n);
      } else if (IDEFMLD == 1) {
        formld_te(z, t, DELTEMLD, amld, n);
      } else if (IDEFMLD == 2) {
        formld_rh(z, rh, DELRHMLD, amld, n);
      }
    } else if (n == 0) {
      amld = z[0];
    } else {
      amld = 0.0;
    }
  }

  double al0 = 0.17 * amld;

  if (IFEPSON2 == 2) {
    ifbelow = 0;
  }


  int icall(0), ipoint(0);
  int iproblem(0), inegproblem(0);

  double eps_bot;

  double sm, sh, ss;
  double slq2, epson2;

  double ss_last, sh_last;

  vector<double> epsy(na);
  vector<bool> lifepsy(na);

//------------------
  for (int k = 0; k < n; ++k) {

    double amtaun2(0.0);

    if ((IFEPSON2 == 2) && (IFDEEPLAT > 0)) {
      ifnofsmall = 0;
      if (an2[k] >= 0.0) {
        an = sqrt(an2[k]);
      }
      if ((an / fabs(coriol)) < 1.0) {
        ifnofsmall = 1;
      }
    }

    ri1 = ri[k];

    if (IFSALI == 0) {
      if (ri1 > rimax) {
        ri[k] = rimax;
        if (IFBACK == 0) {
          akm[k] = v_back[k];
          akh[k] = t_back[k];
          aks[k] = akh[k];
          continue;
        } else {
          sm = sma[NTBL - 1];
          sh = sha[NTBL - 1];
          ss = sh;
          slq2 = slq2a[NTBL - 1];
        }
      } else if (ri1 <= RI0) {
        sm = sma[0];
        sh = sha[0];
        ss = sh;
        slq2 = slq2a[0] * ria[0] / ri1;
      } else {
        const int m = static_cast<int>((ri1 - RI0) / dri) + 1;
        const double tmp = (ri1 - ria[m-1]) / (ria[m] - ria[m-1]);
        sm = sma[m-1] + (sma[m] - sma[m-1]) * tmp;
        sh = sha[m-1] + (sha[m] - sha[m-1]) * tmp;
        slq2 = slq2a[m-1] + (slq2a[m] - slq2a[m-1]) * tmp;
	      ss = sh;
      } 
    } else if (IFSALI == 1) {
      rid1 = rid[k];
      const double and2 = rid[k] / (ri[k] + 1.0e-25) * an2[k];
      double and2on2 = and2 / (an2[k] + 1.0e-25);

      if ((IFZEROSHEAR) && (an2[k] < rib[0] * s2[k])) {
        ifpureshear = 1;
      } else {
        ifpureshear = 0;
      }
      if (ifpureshear == 1) {
        const int imax = MT;
        if (IFASTEXPABS == 0) {
          //interp1d();
        } else if (IFASTEXPABS == 1) {
          interp1d_expabs(and2on2, and2on2a1, amtaun2a1,
              sma1, sha1, ssa1, amtaun2, sm, sh, ss, imax, 
                  MT, MT0, dand2on2, rnd2on2);
        }
        slq2 = (-amtaun2) / ((b1 * b1) * ri[k] + 1.0e-25);
      }
      if (ifpureshear != 1) {
        if (IFASTEXPABS == 0) {
          //interp2d();
        } else if (IFASTEXPABS == 1) {
          interp2d_expabs(ri1, rid1, slq2, sm, sh, ss);
        }
      }
    }

    if (ifpureshear != 1) {
      if (slq2 < 0.0) {
        if (mytid == 0) {
          printf("************************************************\n");
          printf("Error detected in turbulence module.\n");
          printf("'slq2' negative in turb_2 subroutine after interpolation.\n");
          printf("k = %d, slq2 = %f\n", k, slq2);
          printf("sm = %f, sh = %f, ss = %f\n", sm, sh, ss);
          printf("ri1 = %f, rid1 = %f\n", ri1, rid1);
          printf("dri = %f\n", dri);
          printf("Program will stop.");
          printf("************************************************");
        }
        exit(0);
      }
      if (slq2 == 0) {
        ifbelow = 1;
      }
    }

    bool lifupper = (((IFEPSON2  < 2) || (ifbelow == 0) || 
        ((IFDEEPLAT > 0) && (ifnofsmall == 1)) || 
        ((ri1 < 0.0) && (k <= 2))) && (slq2 > 0.0));

    if (lifupper) {
      //ilmldpoint = 1;
      if (ri1 < 0.0) {
        //ilmldpointneg = 1;
      }
    }

    ++ipoint;
    if (k == 0) {
      ++icall;
    }
    epsy[k] = 0.0;
    lifepsy[k] = false;

    if ((isurfuse == 1) && n > 0) {
      if (slq2 == 0.0) {
        epsy[k] = 0.0;
      } else {
        epsy[k] = - buoytot
            / ((1.0 / ((b1 * b1) * slq2)) - 0.5 * sm);
      }
      lifepsy[k] = ((epsy[k] >= 0.0) && lifupper);
      if ((epsy[k] < 0.0) && lifupper) {
        ++iproblem;
        if (ri1 < 0.0) {
          ++inegproblem;
        }
      }
    }

    double akz = 0.4 * z[k];
    double al = akz * al0 / (al0 + akz);

    double rlomega(0.0);

    if (ILOMEGA == 1) {
      if ((buoytot < 0.0) && (amld >= AMLDMINLOM)) {
        rlomega = sqrt(pow(coriol, 3) / (-buoytot));
      } else {
        rlomega = 0.0;
      }
      const double rlblackadar = 1.0 / al;
      const double rl = rlblackadar + rlomega;
      al = 1.0 / rl;
    }

    double al2 = al * al;

    if (!(((IFEPSON2 == 2) && (ifbelow == 1)) || lifepsy[k])) {
      if (ri1 > 0.0) {
        const double anlq2 = slq2 * ri1;
        if ((anlq2 > 0.281) && (ICONDEAR >= 0)) {
          al2 = 0.281 / anlq2 * al2;
          if (ICONDEAR == 0) {
            slq2 = 0.281 / (ri1 + 1.0e-20);
          }
        }
      }
    }

    double epson2_;
    if (an2[k] < 0.0) {
      epson2_ = EPSON2__;
    } else {
      if (IFDEEPLAT == 0) {
        epson2_ = EPSON2__;
      } else if ((IFDEEPLAT == 1) || (IFDEEPLAT == 2)) {
        double eplatidepend;
        if (ifnofsmall == 1) {
          eplatidepend = 0.0;
        } else {
	        eplatidepend = eplatidepend_(fabs(coriol), an);
        }

        eplatidepend = fmax(eplatidepend, EPLATIDEPENDMIN);

        epson2_ = EPSON2__ * eplatidepend;
      }
    }

    if (IFEPSON2 >= 1) {
      if (IFBOTENHANCE == 0) {
        epson2 = epson2_;
      } else if (IFBOTENHANCE == 1) {
        eps_bot = EPS_BOT0 
            * exp((z[k] - z[kbot - 1]) / SCALE_BOT);
        const double epson2_bot = eps_bot / (ri[k] * s2[k]);
        epson2 = fmax(epson2_, epson2_bot);
      }
    }

    if ((IFBACK >= 4) && (IFSALI == 0)) {
      double back_ri1;
      if (IFBACK == 4) {
        back_ri1 = RI_INTERNAL;
      } else {
        back_ri1 = BACKFRAC * rimax;
      }
      const int m = static_cast<int>((back_ri1 - RI0) / dri) + 1;
      const double tmp = (back_ri1 - ria[m-1]) / (ria[m] - ria[m-1]);
      sm_back = sma[m-1] + (sma[m] - sma[m-1]) * tmp;
      sh_back = sha[m-1] + (sha[m] - sha[m-1]) * tmp;
      slq2_back = slq2a[m-1] + (slq2a[m] - slq2a[m-1]) * tmp;
      ss_back = sh_back;

      if (mytid == 0) {
        printf("i = %d, j = %d, m = %d\n", ii, jj, m);
        printf("\n");
        printf("tmp = %f\n", tmp);
        printf("\n");
        printf("sm_back = %f, sma(m) = %f, sma(m+1) = %f\n",
            sm_back, sma[m-1], sma[m]);
        printf("\n");
        printf("sh_back = %f, sha(m) = %f, sha(m+1) = %f\n",
            sh_back, sha[m-1], sha[m]);
        printf("\n");
        printf("slq2_back = %f, slq2a(m) = %f, slq2a(m+1) = %f\n",
            slq2_back, slq2a[m-1], slq2a[m]);
        printf("\n");
      }

      if (slq2_back < 0) {
        if (mytid == 0) {
          printf("************************************************\n");
          printf("Error detected in turbulence module.\n");
          printf("'slq2_back' negative in turb_2 subroutine after 1D interpolation of background.\n");
          printf("k = %d, slq2_back = %f\n", k, slq2_back);
          printf("sm_back = %f, sh_back = %f\n", sm_back, sh_back);
          printf("back_ri1 = %f\n", back_ri1);
          printf("dri = %f\n", dri);
          printf("Program will stop.\n");
          printf("************************************************\n");
        }
        exit(0);
      }

      double s2_back = (ri1 / back_ri1) * s2[k];

      if (ri1 <= 0.0) {
        s2_back = 0.0;
      }
      if (ri1 < 0.0) {
        sm_back = 0.0;
        sh_back = 0.0;
        ss_back = 0.0;
      }
      double tmp_back;
      if (IFEPSON2 == 0) {
        const double al0_back = back_l_0;
        akz = 0.4 * z[k];
        al_back  = akz * al0_back / (al0_back + akz);
        double al2_back = al_back * al_back;
        if (back_ri1 > 0.0) {
          const double anlq2_back = slq2_back * back_ri1;
          if ((anlq2_back > 0.281) && (ICONDEAR >= 0)) {
            al2_back = 0.281 / anlq2_back * al2_back;
            if (ICONDEAR == 0) {
              slq2_back = 0.281 / (back_ri1 + 1.0e-20); 
            }
          }
        }
        tmp_back = 0.5 * b1 * al2_back *
            sqrt(s2_back / (slq2_back + 1.0e-40));
      } else if (IFEPSON2 >= 1) {
        al_back = 0.0;
        tmp_back = 0.5 * pow(b1, 2) * back_ri1 * slq2_back * epson2;
      }
      v_back[k] = tmp_back * sm_back;
      t_back[k] = tmp_back * sh_back;
      s_back[k] = tmp_back * ss_back;
    }

        int ifrafglt = 0;

        double rit(0.0), ric(0.0);
        double s2_back(0.0);
        double back_ra_r1;
        double back_ri1,  back_ric1;
        double back_rid1, back_rit1;
        double theta_r(0);
        double deltheta_r1;
        int itheta_r0(0), itheta_r1(0);
        int jtheta_r0(0), jtheta_r1(0);
        double theta_r_deg(0);
        double ra_r;

    if (IFSALI > 0) {
      if (IFSALBACK == 1) {
        if (slq2 != 0.0) {
          sh_last = sh;
          ss_last = ss;
        } else if (k == 0) {
          s_back[k] = t_back[k];
        }
        if (k != 0){
          s_back[k] = (ss_last / sh_last) * t_back[k];
        }

      } else if (IFSALBACK == 2) {
        double sisa;
        if (slq2 != 0.0) {
          sisa = ss / sh;
        } else {
          rit = (ri[k] + rid[k]) / 2.0;
          ric = (ri[k] - rid[k]) / 2.0;
          
          if (rit == 0.0) {
            if (ric == 0.0) {
              theta_r = atan(1.0);
            } else {
              theta_r = PI / 2.0;
            }
          } else {
            theta_r = atan(ric / rit);
          }
          if (fabs(theta_r) > (PI / 2.0)) {
            exit(0);
          }
          if (theta_r < (-PI) / 4.0) {
            theta_r = theta_r + PI;
          }
          jtheta_r0 = static_cast<int>(
              (theta_r + (PI / 4.0)) / deltheta_r);
          itheta_r0 = jtheta_r0 - N_THETA_R_OCT;
          itheta_r1 = itheta_r0 + 1;
          const double theta_r0  = itheta_r0 * deltheta_r;
          theta_r_deg = theta_r * 180.0 / PI;

          if ((itheta_r1 > 3 * N_THETA_R_OCT) ||
              (itheta_r0 < - N_THETA_R_OCT)) {
            if (mytid == 0) {
              printf("************************************************\n");
              printf("Problem in turbulence module!\n");
              printf("Unrealizability outside Ri>0 region.\n");
              printf("slq2 = %f, sm = %f, sh = %f, ss = %f\n",
                  slq2, sm, sh, ss);
              printf("i = %d, j = %d, k = %d, ri(k) = %f, rid(k) = %f\n",
                  ii, jj, k, ri[k], rid[k]);
              printf("rit = %f, ric = %f, theta_r = %f\n",
                  rit, ric, theta_r);
              printf("theta_r_deg = %f\n", theta_r_deg);
              printf("itheta_r0 = %d, itheta_r1 = %d\n",
                  itheta_r0, itheta_r1);
              printf("n_theta_r_oct = %d\n", N_THETA_R_OCT);
              printf("\n");
              printf("Program will stop.\n");
              printf("\n");
              printf("************************************************\n");
            }
            exit(0);
          }
          deltheta_r1 = theta_r - theta_r0;
          const double delsisa_r = 
              sisamax[N_THETA_R_OCT + itheta_r1] 
                  - sisamax[N_THETA_R_OCT + itheta_r0];
          const double dsisa_o_dtheta = delsisa_r / deltheta_r;
          sisa = sisamax[N_THETA_R_OCT + itheta_r0] 
              + deltheta_r1 * dsisa_o_dtheta;
        }
        s_back[k] = sisa * t_back[k];
      } else if (IFSALBACK == 3) {
        back_ri1  = ri1  * s2[k] * BACK_SM2;
        back_rid1 = rid1 * s2[k] * BACK_SM2;
        if (IFASTEXPABS == 0) {
          //interp2d();
        } else if (IFASTEXPABS == 1) {
          interp2d_expabs(back_ri1, back_rid1,
              slq2_back, sm_back, sh_back, ss_back);
        }
        if (slq2_back < 0) {
          if (mytid == 0) {
            printf("************************************************\n");
            printf("Error detected in turbulence module.\n");
            printf("'slq2_back' negative in turb_2 subroutine after interpolation of background.\n");
            printf("k = %d, slq2_back = %f\n", k, slq2_back);
            printf("sm_back = %f, sh_back = %f, ss_back = %f\n",
                sm_back, sh_back, ss_back);
            printf("back_ri1 = %f, back_rid1 = %f\n",
                back_ri1, back_rid1);
            printf("dri = %f\n", dri);
            printf("Program will stop.\n");
            printf("************************************************\n");
          }
          exit(0);
        }
        double tmp_back;
        if (IFEPSON2 == 0) {
          const double al0_back = back_l_0;
          akz = 0.4 * z[k];
          al_back = akz * al0_back / (al0_back + akz);
          double al2_back = al_back * al_back;
          if (back_ri1 > 0.0) {
            const double anlq2_back = slq2_back * back_ri1;
            if ((anlq2_back > 0.281) && (ICONDEAR >= 0)) {
              al2_back = 0.281 / anlq2_back * al2_back;
              if (ICONDEAR == 0) {
                slq2_back = 0.281 / (back_ri1 + 1.0e-20);
              }
            }
          }
          tmp_back = 0.5 * b1 * al2_back *
              sqrt(BACK_S2 / (slq2_back + 1.0e-40));
        } else if (IFEPSON2 > 0) {
          tmp_back = 0.5 * pow(b1, 2) 
              * back_ri1 * slq2_back * epson2;
        }
        v_back[k] = tmp_back * sm_back + V_BACK0;
        t_back[k] = tmp_back * sh_back + T_BACK0;
        s_back[k] = tmp_back * ss_back + S_BACK0;
      } else if (IFSALBACK >= 4) {

        if (IFSALBACK == 4) {
          back_ri1 = RI_INTERNAL;
          back_rid1 = (rid1 / ri1) * RI_INTERNAL;
        } else {
          if (ri[k] <= 0.0) {
            back_ra_r1 = 0.0;
            back_rit1  = 0.0;
            back_ric1  = 0.0;
            back_ri1   = 0.0;
            back_rid1  = 0.0;
          } else {
            rit = (ri[k] + rid[k]) / 2.0;
            ric = (ri[k] - rid[k]) / 2.0;
            ra_r = sqrt((rit * rit) + (ric * ric));

            if (rit == 0.0) {
              if (ric == 0.0) {
                theta_r = atan(1.0);
              } else {
                theta_r = PI / 2.0;
              }
            } else {
              theta_r = atan(ric / rit);
            }
            if (fabs(theta_r) > (PI / 2.0)) {
              if (mytid == 0) {
                printf("fabs(theta_r) > (PI / 2.0)\n");
                printf("err in cpp version of turb_2. line 689\n");
              }
              exit(0);
            }
            if (theta_r < (-PI) / 4.0) {
              theta_r += PI;
            }

            jtheta_r0 = static_cast<int>(
                (theta_r + (PI / 4.0)) / deltheta_r);
            jtheta_r1 = jtheta_r0 + 1;

            itheta_r0 = jtheta_r0 - N_THETA_R_OCT;
            itheta_r1 = itheta_r0 + 1;

            double theta_r0 = itheta_r0 * deltheta_r;
            double theta_r1 = itheta_r1 * deltheta_r;

            if ((theta_r0 <= theta_rcrp) &&
                 (theta_r > theta_rcrp)) {
              theta_r    = theta_r1;
              theta_r0   = theta_r1;
              itheta_r0  = itheta_r1;
              itheta_r1 += 1;
              theta_r1  += deltheta_r;
            } else if ((theta_r1 >= theta_rcrn) &&
                       (theta_r  < theta_rcrn)) {
              theta_r    = theta_r0;
              theta_r1   = theta_r0;
              itheta_r1  = itheta_r0;
              itheta_r0 -= 1;
              theta_r0  -= deltheta_r;
            }

            theta_r_deg = theta_r * 180 / PI;
            if ((itheta_r1 > 3 * N_THETA_R_OCT) ||
                (itheta_r0 < - N_THETA_R_OCT)) {
              if (mytid == 0) {
                printf("************************************************\n");
                printf("Problem in turbulence module!\n");
                printf("Unrealizability outside Ri>0 region. \n");
                printf("slq2 = %f, sm = %f, sh = %f, ss = %f\n",
                    slq2, sm, sh, ss);
                printf("mytid = %d, ii = %d, jj = %d, k = %d, ri(k) = %f, rid(k) = %f\n",
                    mytid, ii, jj, k, ri[k], rid[k]);
                printf("rit = %f, ric = %f, theta_r = %f\n",
                    rit, ric, theta_r);
                printf("theta_r_deg = %f\n", theta_r_deg);
                printf("itheta_r0 = %d, itheta_r1 = %d\n",
                    itheta_r0, itheta_r1);
                printf("n_theta_r_oct = %d\n", N_THETA_R_OCT);
                printf("pi = %f\n", PI);
                printf("\n");
                printf("Program will stop.\n");
                printf("************************************************\n");
              }
              exit(0);
            }

            deltheta_r1 = theta_r - theta_r0;

            const double delback_ra_r = 
                back_ra_r[N_THETA_R_OCT + itheta_r1] 
                    - back_ra_r[N_THETA_R_OCT + itheta_r0];

            const double dback_ra_r_o_dtheta = 
                delback_ra_r / deltheta_r;

            back_ra_r1 = back_ra_r[N_THETA_R_OCT + itheta_r0] 
                + deltheta_r1 * dback_ra_r_o_dtheta;

            ifrafglt = 0;
            if (IFRAFGMAX == 1) {
              if ((theta_r <= theta_rcrp) ||
                  (theta_r >= theta_rcrn)) {
                if (back_ra_r1 > ra_r) {
                  ifrafglt = 1;
                  back_ra_r1 = ra_r;
                }
              }
            }

            if (back_ra_r1 < 0.0) {
              if (mytid == 0) {
                printf("************************************************\n");
                printf("Problem in turbulence module!\n");
                printf("Negative bg ra_r \\equiv (Ri_T^2+Ri_C^2)^(1/2)\n");
                printf("back_ra_r1 = %f\n", back_ra_r1);
                printf("theta_r = %f\n", theta_r);
                printf("\n");
                printf("slq2 = %f, sm = %f, sh = %f, ss = %f\n",
                    slq2, sm, sh, ss);
                printf("k = %d, ri(k) = %f, rid(k) = %f\n",
                    k, ri[k], rid[k]);
                printf("rit = %f, ric = %f\n", rit, ric);
                printf("itheta_r0 = %d, itheta_r1 = %d\n",
                    itheta_r0, itheta_r1);
                printf("jtheta_r0 = %d, jtheta_r1 = %d\n",
                    jtheta_r0, jtheta_r1);
                printf("theta_r_deg = %f\n", theta_r_deg);
                printf("n_theta_r_oct = %d\n", N_THETA_R_OCT);
                printf("\n");
                printf("\n");
                printf("Program will stop.\n");
                printf("************************************************\n");
              }
              exit(0);
            }
            back_rit1 = cos(theta_r) * back_ra_r1;
            back_ric1 = sin(theta_r) * back_ra_r1;
            back_ri1  = back_rit1 + back_ric1;
            back_rid1 = back_rit1 - back_ric1;
          }
        }

        if ((IFSALBACK != 4) && (ri[k] > 0.0)) {

          if ((IFBG_THETA_INTERP == 0) ||
              (ifrafglt == 1)) {
            if (IFASTEXPABS == 0) {
              // interp2d();
            } else if (IFASTEXPABS == 1) {
              interp2d_expabs(back_ri1, back_rid1,
                  slq2_back, sm_back, sh_back, ss_back);
            }
          } else if(IFBG_THETA_INTERP == 1) {
            deltheta_r1 = theta_r - itheta_r0 * deltheta_r;

            const double delsm_back = sm_r1[N_THETA_R_OCT + itheta_r1]
                - sm_r1[N_THETA_R_OCT + itheta_r0];

            const double dsm_back_o_dtheta = delsm_back / deltheta_r;

            sm_back = sm_r1[N_THETA_R_OCT + itheta_r0]
                + deltheta_r1 * dsm_back_o_dtheta;

            const double delsh_back = sh_r1[N_THETA_R_OCT + itheta_r1]
                - sh_r1[N_THETA_R_OCT + itheta_r0];

            const double dsh_back_o_dtheta = delsh_back / deltheta_r;

            sh_back = sh_r1[N_THETA_R_OCT + itheta_r0]
                + deltheta_r1 * dsh_back_o_dtheta;

            const double delss_back = ss_r1[N_THETA_R_OCT + itheta_r1]
                - ss_r1[N_THETA_R_OCT + itheta_r0];

            const double dss_back_o_dtheta = delss_back / deltheta_r;

            ss_back = ss_r1[N_THETA_R_OCT + itheta_r0]
                + deltheta_r1 * dss_back_o_dtheta;

            const double delslq2_back = slq2_r1[N_THETA_R_OCT + itheta_r1]
                - slq2_r1[N_THETA_R_OCT + itheta_r0];

            const double dslq2_back_o_dtheta = delslq2_back / deltheta_r;

            slq2_back = slq2_r1[N_THETA_R_OCT + itheta_r0]
                + deltheta_r1 * dslq2_back_o_dtheta;
          } else {
            if (mytid == 0) {
              printf("Problem with choice of background interpolation.\n");
              printf("ifbg_theta_interp = %d\n", IFBG_THETA_INTERP);
              printf("ifrafglt = %d\n", ifrafglt);
              printf("Program is stopping.\n");
            }
            exit(0);
          }
          if (slq2_back < 0.0) {
            if (mytid == 0) {
              printf("************************************************\n");
              printf("Error detected in turbulence module.\n");
              printf("'slq2_back' negative in turb_2 subroutine after interpolation of background.\n");
              printf("k = %d, slq2_back = %f\n", k, slq2_back);
              printf("sm_back = %f, sh_back = %f, ss_back = %f\n",
                sm_back, sh_back, ss_back);
              printf("back_ri1 = %f, back_rid1 = %f\n",
                  back_ri1, back_rid1);
              printf("dri = %f\n", dri);
              printf("Program will stop.\n");
              printf("************************************************\n");
            }
            exit(0);
          }
          s2_back = (ri1 / back_ri1) * s2[k];
        } // not go to 19
        // 19 continue
        if (ri1 <= 0.0) {
          s2_back = 0.0;
        }
        if (ri1 < 0.0) {
          sm_back = 0.0;
          sh_back = 0.0;
          ss_back = 0.0;
        }

        if ((sm_back < 0.0) || (sh_back < 0.0) ||
            (ss_back < 0.0)) {
          if (mytid == 0) {
            printf("************************************************\n");
            printf("Problem in turbulence module!\n");
            printf("Negative Structure Function in Background.\n");
            printf("i = %d, j = %d\n", ii, jj);
            printf("temperature = %f\n", t[k]);
            printf("\n");
            printf("salinity = %f\n", s[k]);
            printf("\n");
            printf("density = %f\n", rh[k]);
            printf("\n");
            printf("richardson = %f\n", ri[k]);
            printf("\n");
            printf("rid = %f\n", rid[k]);
            printf("\n");
            printf("shear = %f\n", s2[k]);
            printf("\n");
            printf("BV frequency = %f\n", an2[k]);
            printf("\n");
            printf("ustart = %f\n", ustar_);
            printf("\n");
            printf("buoytur = %f\n", buoytur);
            printf("\n");
            printf("buoysol = %f\n", buoysol);
            printf("slq2_back = %f\n", slq2_back);
            printf("sm_back = %f, sh_back = %f, ss_back = %f\n", 
                sm_back, sh_back, ss_back);
            printf("\n");
            printf("back_ra_r1 = %f\n", back_ra_r1);
            printf("theta_r = %f\n", theta_r);
            printf("back_rit1 = %f, back_ric1 = %f\n", 
                back_rit1, back_ric1);
            printf("back_ri1 = %f, back_rid1 = %f\n", 
                back_ri1, back_rid1);
            printf("\n");
            printf("itheta_r0 = %d, itheta_r1 = %d\n",
                itheta_r0, itheta_r1);
            printf("jtheta_r0 = %d, jtheta_r1 = %d\n",
                jtheta_r0, jtheta_r1);
            printf("theta_r_deg = %f\n", theta_r_deg);
            printf("n_theta_r_oct = %d\n", N_THETA_R_OCT);
            printf("\n");
            printf("k = %d, ri(k) = %f, rid(k) = %f\n",
                k, ri[k], rid[k]);
            printf("rit = %f, ric = %f\n", rit, ric);
            printf("\n");
            printf("\n");
            printf("Program will stop.\n");
            printf("************************************************\n");
          }
          exit(0);
        }
        // TODO
        double tmp_back;
        if (IFEPSON2 == 0) {
          const double al0_back = back_l_0;
          akz = 0.4 * z[k];
          al_back = akz * al0_back
              / (al0_back + akz);

          double al2_back = al_back * al_back;

          if (back_ri1 > 0.0) {
            const double anlq2_back = slq2_back * back_ri1;
            if ((anlq2_back > 0.281) && (ICONDEAR >= 0)) {
              al2_back = 0.281 / anlq2_back * al2_back;
              if (ICONDEAR == 0) {
                slq2_back = 0.281 / (back_ri1 + 1.0e-20);
              }
            }
          }
          tmp_back = 0.5 * b1 * al2_back
              * sqrt(s2_back / (slq2_back + 1.0e-40));

        } else if (IFEPSON2 > 0) {
          tmp_back = 0.5 * b1 * b1
              * back_ri1 * slq2_back * epson2;
        }
        v_back[k] = tmp_back * sm_back;
        t_back[k] = tmp_back * sh_back;
        s_back[k] = tmp_back * ss_back;

        if ((v_back[k] < 0.0) || (t_back[k] < 0.0)
            || (s_back[k] < 0.0)) {
          if (mytid == 0) {
            printf("************************************************\n");
            printf("Problem in turbulence module!\n");
            printf("Negative Background Diffusivity.\n");
            printf("i = %d, j = %d\n", ii, jj);
            printf("temperature = %f\n", t[k]);
            printf("\n");
            printf("salinity = %f\n", s[k]);
            printf("\n");
            printf("density = %f\n", rh[k]);
            printf("\n");
            printf("richardson = %f\n", ri[k]);
            printf("\n");
            printf("rid = %f\n", rid[k]);
            printf("\n");
            printf("shear = %f\n", s2[k]);
            printf("\n");
            printf("BV frequency = %f\n", an2[k]);
            printf("\n");
            printf("ustart = %f\n", ustar_);
            printf("\n");
            printf("buoytur = %f\n", buoytur);
            printf("\n");
            printf("buoysol = %f\n", buoysol);
            printf("slq2_back = %f\n", slq2_back);
            printf("v_back = %f, t_back = %f, s_back = %f\n", 
                v_back[k], t_back[k], s_back[k]);
            printf("\n");
            printf("slq2_back = %f\n", slq2_back);
            printf("sm_back = %f, sh_back = %f, ss_back = %f\n", 
                sm_back, sh_back, ss_back);
            printf("\n");
            printf("back_ra_r1 = %f\n", back_ra_r1);
            printf("theta_r = %f, theta_r_deg = %f\n", 
                theta_r, theta_r_deg);
            printf("back_rit1 = %f, back_ric1 = %f\n", 
                back_rit1, back_ric1);
            printf("back_ri1 = %f, back_rid1 = %f\n", 
                back_ri1, back_rid1);
            printf("\n");
            printf("k = %d, ri(k) = %f, rid(k) = %f\n",
                k, ri[k], rid[k]);
            printf("rit = %f, ric = %f\n", rit, ric);
            printf("\n");
            printf("\n");
            printf("Program will stop.\n");
            printf("************************************************\n");
          }
          exit(0);
        }
        if((ri[k] > 0.0) && ((v_back[k] == 0.0)
            || (t_back[k] == 0.0) || (s_back[k] == 0.0))) {
          if (mytid == 0) {
            printf("************************************************\n");
            printf("Problem in turbulence module!\n");
            printf("Zero Background Diffusivity in stable case.\n");
            printf("i = %d, j = %d\n", ii, jj);
            printf("temperature = %f\n", t[k]);
            printf("\n");
            printf("salinity = %f\n", s[k]);
            printf("\n");
            printf("density = %f\n", rh[k]);
            printf("\n");
            printf("richardson = %f\n", ri[k]);
            printf("\n");
            printf("rid = %f\n", rid[k]);
            printf("\n");
            printf("shear = %f\n", s2[k]);
            printf("\n");
            printf("BV frequency = %f\n", an2[k]);
            printf("\n");
            printf("ustart = %f\n", ustar_);
            printf("\n");
            printf("buoytur = %f\n", buoytur);
            printf("\n");
            printf("buoysol = %f\n", buoysol);
            printf("v_back = %f, t_back = %f, s_back = %f\n", 
                v_back[k], t_back[k], s_back[k]);
            printf("\n");
            printf("slq2_back = %f\n", slq2_back);
            printf("sm_back = %f, sh_back = %f, ss_back = %f\n", 
                sm_back, sh_back, ss_back);
            printf("\n");
            printf("slq2_r1(itheta_r0) = %f, slq2_r1(itheta_r1) = %f\n",
                slq2_r1[N_THETA_R_OCT + itheta_r0], 
                    slq2_r1[N_THETA_R_OCT + itheta_r1]);
            printf("sm_r1(itheta_r0) = %f, sm_r1(itheta_r1) = %f\n",
                sm_r1[N_THETA_R_OCT + itheta_r0], 
                    sm_r1[N_THETA_R_OCT + itheta_r1]);
            printf("sh_r1(itheta_r0) = %f, sh_r1(itheta_r1) = %f\n",
                sh_r1[N_THETA_R_OCT + itheta_r0], 
                    sh_r1[N_THETA_R_OCT + itheta_r1]);
            printf("ss_r1(itheta_r0) = %f, ss_r1(itheta_r1) = %f\n",
                ss_r1[N_THETA_R_OCT + itheta_r0], 
                    ss_r1[N_THETA_R_OCT + itheta_r1]);
            printf("\n");
            printf("back_ra_r1 = %f\n", back_ra_r1);
            printf("theta_r = %f\n", theta_r);
            printf("theta_r_deg = %f\n", theta_r_deg);
            printf("back_rit1 = %f, back_ric1 = %f\n", 
                back_rit1, back_ric1);
            printf("back_ri1 = %f, back_rid1 = %f\n", 
                back_ri1, back_rid1);
            printf("itheta_r0 = %d, itheta_r1 = %d\n",
                itheta_r0, itheta_r1);
            printf("jtheta_r0 = %d, jtheta_r1 = %d\n",
                jtheta_r0, jtheta_r1);
            printf("n_theta_r_oct = %d\n", N_THETA_R_OCT);
            printf("deltheta_r = %f\n", deltheta_r);
            printf("\n");
            printf("\n");
            printf("k = %d, ri(k) = %f, rid(k) = %f\n",
                k, ri[k], rid[k]);
            printf("rit = %f, ric = %f\n", rit, ric);
            printf("\n");
            printf("\n");
            printf("Program will stop.\n");
            printf("************************************************\n");
          }
          exit(0);
        }
      } // ifsalback >= 4
      // 20 continue
    } // end ifsali.GT.0

    if (ifoutput == 1) {
      if (mytid == 0) {
        if (isurfuse == 0) {
          printf("z(k) = %f\n", z[k]);
          printf("al = %f\n", al);
          printf("slq2 = %f\n", slq2);
          printf("ri1 = %f\n", ri1);
          printf("rid1 = %f\n", rid1);
          printf("sm = %f\n", sm);
          printf("sh = %f\n", sh);
          printf("ss = %f\n", ss);
          printf("v_back(k) = %f\n", v_back[k]);
          printf("t_back(k) = %f\n", t_back[k]);
          printf("s_back(k) = %f\n", s_back[k]);
        } else if (isurfuse == 1) {
          printf("z(k) = %f\n", z[k]);
          printf("al = %f\n", al);
          printf("slq2 = %f\n", slq2);
          printf("ri1 = %f\n", ri1);
          printf("rid1 = %f\n", rid1);
          printf("sm = %f\n", sm);
          printf("sh = %f\n", sh);
          printf("ss = %f\n", ss);
          printf("epsy(k) = %f\n", epsy[k]);
          printf("v_back(k) = %f\n", v_back[k]);
          printf("t_back(k) = %f\n", t_back[k]);
          printf("s_back(k) = %f\n", s_back[k]);
        }
      }
      if ((IFBACK == 4) || (IFSALBACK == 4)) {
        if (mytid == 0) {
          if (k == 0) {
            printf("\n");
            printf("z[cm]      l_back[cm] Ri-table   Ri_d-table Ri_d_back  s2_back    slq2_back  sm_back    sh_back    ss_back   \n");
          }
          printf("z(k) = %f\n", z[k]);
          printf("al_back = %f\n", al_back);
          printf("ri1 = %f\n", ri1);
          printf("rid1 = %f\n", rid1);
          printf("back_rid1 = %f\n", back_rid1);
          printf("s2_back = %f\n", s2_back);
          printf("slq2_back = %f\n", slq2_back);
          printf("sm_back = %f\n", sm_back);
          printf("sh_back = %f\n", sh_back);
          printf("ss_back = %f\n", ss_back);
        }
      } else if ((IFBACK > 4) || (IFSALBACK > 4)) {
        if (mytid == 0) {
          if (k == 0) {
            printf("\n");
            if (IFEPSON2 == 0) {
            printf("z[cm]      l_back[cm] Ri-table   Ri_d-table \n");
            printf("Ri_back    Ri_d_back  ra_r_back \n");
            printf("s2_back    slq2_back  sm_back    sh_back    ss_back   \n");
            printf("N^2        S^2        \n");
            } else {
              if (IFBOTENHANCE == 0) {
                printf("z[cm]      l_back[cm] Ri-table   Ri_d-table \n");
                printf("Ri_back    Ri_d_back  ra_r_back \n");
                printf("s2_back    slq2_back  sm_back    sh_back    ss_back   \n");
                printf("N^2        S^2        \n");
                printf("epson2     \n");
              } else {
                printf("z[cm]      l_back[cm] Ri-table   Ri_d-table \n");
                printf("Ri_back    Ri_d_back  ra_r_back \n");
                printf("s2_back    slq2_back  sm_back    sh_back    ss_back   \n");
                printf("N^2        S^2        \n");
                printf("epson2     eps_bot\n");
              }
            }
          }
        }
        if (IFEPSON2 == 0) {
          if (mytid == 0) {
            printf("z(k) = %f\n", z[k]);
            printf("al_back = %f\n", al_back);
            printf("ri1 = %f\n", ri1);
            printf("rid1 = %f\n", rid1);
            printf("back_ri1 = %f\n", back_ri1);
            printf("back_rid1 = %f\n", back_rid1);
            printf("back_ra_r1 = %f\n", back_ra_r1);
            printf("s2_back = %f\n", s2_back);
            printf("slq2_back = %f\n", slq2_back);
            printf("sm_back = %f\n", sm_back);
            printf("sh_back = %f\n", sh_back);
            printf("ss_back = %f\n", ss_back);
            printf("ri(k) * s2(k) = %f\n", ri[k] * s2[k]);
            printf("s2(k) = %f\n", s2[k]);
          }
        } else {
          if (IFBOTENHANCE == 0) {
            if (mytid == 0) {
              printf("z(k) = %f\n", z[k]);
              printf("al_back = %f\n", al_back);
              printf("ri1 = %f\n", ri1);
              printf("rid1 = %f\n", rid1);
              printf("back_ri1 = %f\n", back_ri1);
              printf("back_rid1 = %f\n", back_rid1);
              printf("back_ra_r1 = %f\n", back_ra_r1);
              printf("s2_back = %f\n", s2_back);
              printf("slq2_back = %f\n", slq2_back);
              printf("sm_back = %f\n", sm_back);
              printf("sh_back = %f\n", sh_back);
              printf("ss_back = %f\n", ss_back);
              printf("ri(k) * s2(k) = %f\n", ri[k] * s2[k]);
              printf("s2(k) = %f\n", s2[k]);
              printf("epson2 = %f\n", epson2);
            }
          } else {
            if (mytid == 0) {
              printf("z(k) = %f\n", z[k]);
              printf("al_back = %f\n", al_back);
              printf("ri1 = %f\n", ri1);
              printf("rid1 = %f\n", rid1);
              printf("back_ri1 = %f\n", back_ri1);
              printf("back_rid1 = %f\n", back_rid1);
              printf("back_ra_r1 = %f\n", back_ra_r1);
              printf("s2_back = %f\n", s2_back);
              printf("slq2_back = %f\n", slq2_back);
              printf("sm_back = %f\n", sm_back);
              printf("sh_back = %f\n", sh_back);
              printf("ss_back = %f\n", ss_back);
              printf("ri(k) * s2(k) = %f\n", ri[k] * s2[k]);
              printf("s2(k) = %f\n", s2[k]);
              printf("epson2 = %f\n", epson2);
              printf("eps_bot = %f\n", eps_bot);
            }
            if ((IFCHECKBOTTOMEPS == 1) && (k == n)) {
              eps_bot__under = EPS_BOT0
                  * exp((z[k+1] - z[kbot - 1]) / SCALE_BOT);
              if (mytid == 0) {
                printf("z(k+1) = %f\n", z[k+1]);
                printf("eps_bot__under = %f\n", eps_bot__under);
              }
            }
          }
        }
      }
    }

    double aldeep[NBIG];

    aldeep[k] = 0.0;

    double tmp(0.0);

    if ((IFEPSON2 == 2) && (ifbelow == 1)) {
      if ((ri1 >= 0.0) || ((IFDEEPLAT == 2) && (ifnofsmall == 1))) {
        tmp = 0.5 * b1 * b1 * ri1 * slq2 * epson2;
      } else if (k > 1) {
        double delz, delrh, del2rh;
        if (k == n - 1) {
          delz = z[k] - z[k-1];
          delrh = rh[k] - rh[k-1];
          del2rh = rh[k] - 2.0 * rh[k-1] + rh[k-2];
        } else {
          delz = z[k+1] - z[k-1];
          delrh = rh[k+1] - rh[k-1];
          del2rh = rh[k+1] - 2.0 * rh[k] + rh[k-1];
        }

        const double dzrh = delrh / delz;
        const double d2zrh = 4.0 * del2rh / (delz * delz);
        const double rdzlndzrh = dzrh / d2zrh;
        const double al0deep = 0.17 * fabs(rdzlndzrh);
        akz = 0.4 * z[k];
        aldeep[k] = akz * al0deep / (al0deep + akz);
        al2 = aldeep[k] * aldeep[k];

        if (ifpureshear == 1) {
          // GO TO 21
          tmp = 0.5 * (b1 * b1) * al2
              * sqrt(-an2[k] / amtaun2);
     
          akm[k] = fmin(tmp * sm + v_back[k],
              visc_cbu_limit1);
          akh[k] = fmin(tmp * sh + t_back[k],
              diff_cbt_limit1);
          aks[k] = fmin(tmp * ss + s_back[k],
              diff_cbt_limit1);
          continue ;
        } else if (IFSHEARMIN) {
          s2[k] = fmax(s2[k], S2MIN);
        }
        tmp = 0.5 * b1 * al2
            * sqrt(s2[k] / (slq2 + 1.0e-40));
      } else {
        if (ifpureshear == 1) {
          // GO TO 21
          tmp = 0.5 * (b1 * b1) * al2
              * sqrt(-an2[k] / amtaun2);
     
          akm[k] = fmin(tmp * sm + v_back[k],
              visc_cbu_limit1);
          akh[k] = fmin(tmp * sh + t_back[k],
              diff_cbt_limit1);
          aks[k] = fmin(tmp * ss + s_back[k],
              diff_cbt_limit1);
          continue ;
        } else if (IFSHEARMIN) {
          s2[k] = fmax(s2[k], S2MIN);
        }
        if (lifepsy[k]) {
          tmp = 0.5 * epsy[k] / (s2[k] + 1.0e-40);
        } else {
          tmp = 0.5 * b1 * al2
              * sqrt(s2[k] / (slq2 + 1.0e-40));
        }
      }
    } else {
      if (ifpureshear == 1) {
          if (ifpureshear == 1) {
            tmp = 0.5 * (b1 * b1) * al2
                * sqrt(-an2[k] / amtaun2);
          }
          akm[k] = fmin(tmp * sm + v_back[k],
              visc_cbu_limit1);
          akh[k] = fmin(tmp * sh + t_back[k],
              diff_cbt_limit1);
          aks[k] = fmin(tmp * ss + s_back[k],
              diff_cbt_limit1);
          continue ;
      } else if (IFSHEARMIN) {
        s2[k] = fmax(s2[k], S2MIN);
      }
      if (lifepsy[k]) {
        tmp = 0.5 * epsy[k] / (s2[k] + 1.0e-40);
      } else {
        tmp = 0.5 * b1 * al2
            * sqrt(s2[k] / (slq2 + 1.0e-40));
      }
    }

    // 21
    if (ifpureshear == 1) {
      tmp = 0.5 * (b1 * b1) * al2
          * sqrt(-an2[k] / amtaun2);
    }

    akm[k] = fmin(tmp * sm + v_back[k],
        visc_cbu_limit1);
    akh[k] = fmin(tmp * sh + t_back[k],
        diff_cbt_limit1);
    aks[k] = fmin(tmp * ss + s_back[k],
        diff_cbt_limit1);

  }

  for (int k = nb+1; k < nmax; ++k) {
    akm[k] = 0.0;
    akh[k] = 0.0;
    aks[k] = 0.0;
  }

  if (n > 0) {
    if (akm[0] < wndmix) {
      akm[0] = wndmix;
    }
    if (akh[0] < wndmix) {
      akh[0] = wndmix;
    }
    if (aks[0] < wndmix) {
      aks[0] = wndmix;
    }
  }

  // Fortran: COMPLEX*16 zlomega
  //zlomega = sqrt(CMPLX(-buoytot / (pow(coriol, 3))));

  if (ifoutput == 1) {
    if (mytid == 0) {
      // print something
    }
    for (int k = 0; k < n; ++k) {
      if (mytid == 0) {
        // print something
      }
      if ((akm[k] < 0.0) || (akh[k] < 0.0) || (akm[k] < 0.0)) {
        if (mytid == 0) {
          // print something
        }
        exit(0);
      }
    }
  }
  return ;
} 

inline int sign_int(const double &x, const double &y) {
  return y >= 0 ? std::abs(static_cast<int>(x)) 
      : -std::abs(static_cast<int>(x));
}

inline float sign_float(const double &x, const double &y) {
  return y >= 0 ? std::abs(static_cast<float>(x)) 
      : -std::abs(static_cast<float>(x));
}

inline double sign_double(const double &x, const double &y) {
  return y >= 0.e0 ? fabs(x) : -fabs(x);
}

static void formld(const double *z, const double *t,
    double &amld, const int &n) {
  for (int k = 0; k < n; ++k) {
    if (fabs(t[k] - t[0]) > 0.1) {
#ifdef D_PRECISION
      const double tm = t[0] - sign_double(0.1e0, t[0] - t[k]);
#else
      const float tm = t[0] - sign_float(0.1, t[0] - t[k]);
#endif // D_PRECISION
      amld = z[k] + (z[k-1] - z[k]) *
          (tm - t[k]) / (t[k-1] - t[k] + 1.e-20);
      return ;
    }
  }
  amld = z[n-1];
  return ;
}

static void formld_te(const double *z, const double *t, 
    const double &delte, double &amld, const int &n) {
  for (int k = 0; k < n; ++k) {
    if (fabs(t[k] - t[0]) > delte) {
      const double tm = t[0] - sign_double(delte, t[0] - t[k]);
      amld = z[k] + (z[k-1] - z[k]) *
          (tm - t[k]) / (t[k-1] - t[k] + 1.e-20);
      return ;
    }
  }
  amld = z[n-1];
  return ;
}

static void formld_rh(const double *z, const double *t, 
    const double &delrh, double &amld, const int &n) {
  for (int k = 0; k < n; ++k) {
    if (fabs(t[k] - t[0]) > delrh) {
      const double tm = t[0] - sign_double(delrh, t[0] - t[k]);
      amld = z[k] + (z[k-1] - z[k]) * 
          (tm- t[k]) / (t[k-1] - t[k] + 1.e-20);
      return ;
    }
  }
  amld = z[n-1];
  return ;
}

static void interp1d_expabs(
    double &x,            // and2on2
    const double *x_1,    // and2on2a1
    const double *slq2_1, // amtaun2a1
    const double *sm_1,   // sma1
    const double *sh_1,   // sha1
    const double *ss_1,   // ssa1
    double &slq2,         // amtaun2
    double &sm,           // sm
    double &sh,           // sh
    double &ss,           // ss
    const int &ixmax,     // imax
    const int &m,         // MT
    const int &m0,        // MT0
    const double &delta,  // dand2on2
    const double &rat     /* rnd2on2 */ ) {

  int lx0(0), lx1(0);
  
  if (x > x_1[MT + ixmax]) {
    x = x_1[MT + m];
  } else if (x < x_1[0]) {
    x = x_1[0];
  }

  if (fabs(x) < x_1[MT + m0]) {
#ifdef D_PRECISION
    lx1 = static_cast<int>(x / delta) + round(sign_double(static_cast<double>(1.0), x)); 
#else  // D_PRECISION
    lx1 = static_cast<int>(x / delta) + round(sign_float(static_cast<float>(1.0), x)); 
#endif // D_PRECISION
  } else if (fabs(x) >= x_1[MT + m]) {
#ifdef D_PRECISION
    lx0 = round(sign_double(static_cast<double>(m), x));
#else  // D_PRECISION
    lx0 = round(sign_float(static_cast<float>(m), x));
#endif // D_PRECISION
    lx1 = lx0;
  } else {
#ifdef D_PRECISION
    const double tabindx = sign_double(static_cast<double>(m0)
        + ((log(fabs(x)) - log(x_1[MT + m0])) / log(rat)), x);
#else  // D_PRECISION
    const float tabindx = sign_double(static_cast<float>(m0)
        + ((log(fabs(x)) - log(x_1[MT + m0])) / log(rat)), x);
#endif // D_PRECISION

#ifdef D_PRECISION
    lx1 = static_cast<int>(tabindx) + round(sign_double(
        static_cast<double>(1.0), x));
#else  // D_PRECISION
    lx1 = static_cast<int>(tabindx) + round(sign_float(
        static_cast<float>(1.0), x));
#endif // D_PRECISION
  }

  if (!(fabs(x) >= x_1[MT + m])) {
    if (fabs(x_1[MT + lx1]) < fabs(x)) {
      lx1 += sign_int(1, lx1);
    } else if (fabs(x_1[MT + lx1 - sign_int(1, lx1)]) > fabs(x)) {
      lx1 -= sign_int(1, lx1);
    }

#ifdef D_PRECISION
    lx0 = lx1 - round(sign_double(
        static_cast<double>(1.0), x));
#else  // D_PRECISION
    lx0 = lx1 - round(sign_float(
        static_cast<float>(1.0), x));
#endif // D_PRECISION
    if (x == 0.0) {
      lx1 = 1;
    }
  }
  if ((x > 0.0 && (x < x_1[MT + lx0] || x > x_1[MT + lx1])) ||
      (x < 0.0 && (x > x_1[MT + lx0] || x < x_1[MT + lx1]))) {
    if (CppParamMod::mytid == 0) {
      printf("x is outside interpolation range in interp1d_expabs.\n");
      printf("delta = %f\n", delta);
      printf("m0 = %d, m = %d, rat = %f\n", m0, m, rat);
      printf("x = %f, lx0 = %d, lx1 = %d\n", x, lx0, lx1);
      printf("x_1(lx0) = %f, x_1(lx1) = %f\n", x_1[MT + lx0], x_1[MT + lx1]);
      printf("Program is stopping.\n");
    }
    exit(0);
  }

  const double deltaxta = x_1[MT + lx1] - x_1[MT + lx0];

  const double deltax = x - x_1[MT + lx0];

  double dslq2_x;
  if (lx1 == lx0) {
    dslq2_x = 0.0;
  } else {
    dslq2_x = (slq2_1[MT + lx1] - slq2_1[MT + lx0])
        / deltaxta;
  }

  slq2 = slq2_1[MT + lx0] + dslq2_x * deltax;

  double dsm_x;
  if (lx1 == lx0) {
    dsm_x = 0.0;
  } else {
    dsm_x = (sm_1[MT + lx1] - sm_1[MT + lx0])
        / deltaxta;
  }
  sm = sm_1[MT + lx0] + dsm_x * deltax;

  double dsh_x;
  if (lx1 == lx0) {
    dsh_x = 0.0;
  } else {
    dsh_x = (sh_1[MT + lx1] - sh_1[MT + lx0])
        / deltaxta;
  }
  sh = sh_1[MT + lx0] + dsh_x * deltax;

  double dss_x;
  if (lx1 == lx0) {
    dss_x = 0.0;
  } else {
    dss_x = (ss_1[MT + lx1] - ss_1[MT + lx0])
        / deltaxta;
  }
  ss = ss_1[MT + lx0] + dss_x * deltax;
  
  return ;  
}

static void interp2d_expabs(double &ri, double &rid,
    double &slq2, double &sm, double &sh, double &ss) {

  using CppCanutoMod::dri;
  using CppCanutoMod::shb;
  using CppCanutoMod::smb;
  using CppCanutoMod::ssb;
  using CppCanutoMod::rri;
  using CppCanutoMod::rib;
  using CppCanutoMod::ridb;
  using CppCanutoMod::slq2b;
  using CppCanutoMod::irimax;

  if (ri > rib[MT + MT]) {
    if (fabs(rid) <= ri) {
      rid = rib[MT + MT] * (rid / ri);
      ri  = rib[MT + MT];
    } else if (rid > ri) {
      ri  = ridb[MT + MT] * (ri / rid);
      rid = ridb[MT + MT];
    } else if (rid > -ri) {
      ri  = ridb[0] * (ri / rid);
      rid = ridb[0];
    }
  } else if (ri < rib[0]) {
    if (fabs(rid) < -ri) {
      rid = rib[0] * (rid / ri);
      ri  = rib[0];
    } else if (rid > -ri) {
      ri  = ridb[MT + MT] * (ri / rid);
      rid = ridb[MT + MT];
    } else if (rid < ri) {
      ri  = ridb[0] * (ri / rid);
      rid = ridb[0];
    }
  } else if (rid > ridb[MT + MT]) {
    ri  = ridb[MT + MT] * (ri / rid);
    rid = ridb[MT + MT];
  } else if (rid < ridb[0]) {
    ri  = ridb[0] * (ri / rid);
    rid = ridb[0];
  }

  int lrid0(0), lrid1(0);

  if (fabs(rid) < ridb[MT + MT0]) {
#ifdef D_PRECISION
    lrid1 = static_cast<int>(rid / dri)
        + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
    lrid1 = static_cast<int>(rid / dri)
        + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION
  } else if (fabs(rid) >= ridb[MT + MT]) {
#ifdef D_PRECISION
    lrid0 = round(sign_double(static_cast<double>(MT), rid));
#else // D_PRECISION
    lrid0 = round(sign_float(static_cast<float>(MT), rid));
#endif // D_PRECISION
    lrid1 = lrid0;
  } else {
#ifdef D_PRECISION
    const double tabindrid = sign_double(static_cast<double>(MT0)
        + ((log(fabs(rid)) - log(ridb[MT + MT0]))
            / log(rri)), ri);
#else // D_PRECISION
    const double tabindrid = sign_float(static_cast<float>(MT0)
        + ((log(fabs(rid)) - log(ridb[MT + MT0]))
            / log(rri)), ri);
#endif // D_PRECISION

#ifdef D_PRECISION
    lrid1 = static_cast<int>(tabindrid)
        + round(sign_double(static_cast<double>(1.0), rid));
#else // D_PRECISION
    lrid1 = static_cast<int>(tabindrid)
        + round(sign_float(static_cast<float>(1.0), rid));
#endif // D_PRECISION
  }

  if (!(fabs(rid) >= ridb[MT + MT])) {
    if (fabs(ridb[MT + lrid1]) < fabs(rid)) {
      lrid1 += sign_int(1, lrid1);
    } else if (fabs(ridb[MT + lrid1 - sign_int(1, lrid1)]) > 
        fabs(rid)) {
      lrid1 -= sign_int(1, lrid1);
    }

#ifdef D_PRECISION
    lrid0 = lrid1 - round(sign_double(
        static_cast<double>(1.0), rid));
#else // D_PRECISION
    lrid0 = lrid1 - round(sign_float(
        static_cast<float>(1.0), rid));
#endif // D_PRECISION
    if (rid == 0.0) {
      lrid1 = 1;
    }
  }

  if ((rid > 0.0 && (rid < ridb[MT + lrid0] || rid > ridb[MT + lrid1])) ||
      (rid < 0.0 && (rid > ridb[MT + lrid0] || rid < ridb[MT + lrid1]))) {
    if (CppParamMod::mytid == 0) {
      printf("Ri_d is outside interpolation range in interp2d.\n");
      printf("rid = %f, lrid0 = %d, lrid1 = %d\n", rid, lrid0, lrid1);
      printf("rid_1(lrid0) = %f, rid_1(lrid1) = %f\n", 
          ridb[MT + lrid0], ridb[MT + lrid1]);
      printf("Program is stopping.\n");
    }
    exit(0);
  }

  if (ri > fmin(rib[MT + irimax[MT + lrid0]],
                rib[MT + irimax[MT + lrid1]])) {
    slq2 = 0.0;
    sm   = 0.0;
    sh   = 0.0;
    ss   = 0.0;
    return ;
  }

  int lri0(0), lri1(0);

  if (fabs(ri) < rib[MT + MT0]) {
#ifdef D_PRECISION
    lri1 = static_cast<int>(ri / dri)
        + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
    lri1 = static_cast<int>(ri / dri)
        + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
  } else if (fabs(ri) >= rib[MT + MT]) {
#ifdef D_PRECISION
    lri0 = round(sign_double(static_cast<double>(MT), ri));
#else // D_PRECISION
    lri0 = round(sign_float(static_cast<float>(MT), ri));
#endif // D_PRECISION
    lri1 = lri0;
  } else {
#ifdef D_PRECISION
    const double tabindri = sign_double(static_cast<double>(MT0)
        + ((log(fabs(ri)) - log(rib[MT + MT0]))
            / log(rri)), ri);
#else // D_PRECISION
    const double tabindri = sign_float(static_cast<float>(MT0)
        + ((log(fabs(ri)) - log(rib[MT + MT0]))
            / log(rri)), ri);
#endif // D_PRECISION

#ifdef D_PRECISION
    lri1 = static_cast<int>(tabindri)
        + round(sign_double(static_cast<double>(1.0), ri));
#else // D_PRECISION
    lri1 = static_cast<int>(tabindri)
        + round(sign_float(static_cast<float>(1.0), ri));
#endif // D_PRECISION
  }

  if (!(fabs(ri) >= rib[MT + MT])) {
    if (fabs(rib[MT + lri1]) < fabs(ri)) {
      lri1 += sign_int(1, lri1);
    } else if (fabs(rib[MT + lri1 - sign_int(1, lri1)]) >
        fabs(ri)) {
      lri1 -= sign_int(1, lri1);
    }

#ifdef D_PRECISION
    lri0 = lri1 - round(sign_double(
        static_cast<double>(1.0), ri));
#else // D_PRECISION
    lri0 = lri1 - round(sign_float(
        static_cast<float>(1.0), ri));
#endif // D_PRECISION

    if (ri == 0.0) {
      lri1 = 1;
    }
  }

  if ((ri > 0.0 && (ri < rib[MT + lri0] || ri > rib[MT + lri1])) ||
      (ri < 0.0 && (ri > rib[MT + lri0] || ri < rib[MT + lri1]))) {
    if (CppParamMod::mytid == 0) {
      printf("Ri is outside interpolation range in interp2d.\n");
      printf("ri = %f, lri0 = %d, lri1 = %d\n", ri, lri0, lri1);
      printf("ri_1(lri0) = %f, ri_1(lri1) = %f\n", 
          rib[MT + lri0], rib[MT + lri1]);
      printf("Program is stopping.\n");
    }
    exit(0);
  }

  const double deltaridta = ridb[MT + lrid1] - ridb[MT + lrid0];
  const double deltarita = rib[MT + lri1] - rib[MT + lri0];
  const double deltarid = rid - ridb[MT + lrid0];
  const double deltari = ri - rib[MT + lri0];

  double dslq2_rid;

  if (lrid1 == lrid0) {
    dslq2_rid = 0.0;
  } else {
    dslq2_rid = (slq2b[MT+lrid1][MT+lri0] - slq2b[MT+lrid0][MT+lri0])
        / deltaridta;
  }

  double dslq2_ri;
  if (lri1 == lri0) {
    dslq2_ri = 0.0;
  } else {
    dslq2_ri = (slq2b[MT+lrid0][MT+lri1] - slq2b[MT+lrid0][MT+lri0])
        / deltarita;
  }

  slq2 = slq2b[MT+lrid0][MT+lri0]
      + dslq2_ri * deltari + dslq2_rid * deltarid;

  double dsm_rid;
  if (lrid1 == lrid0) {
    dsm_rid = 0.0;
  } else {
    dsm_rid = (smb[MT+lrid1][MT+lri0] - smb[MT+lrid0][MT+lri0])
        / deltaridta;
  }

  double dsm_ri;

  if (lri1 == lri0) {
    dsm_ri = 0.0;
  } else {
    dsm_ri = (smb[MT+lrid0][MT+lri1] - smb[MT+lrid0][MT+lri0])
        / deltarita;
  }

  sm = smb[MT+lrid0][MT+lri0] 
      + dsm_ri * deltari + dsm_rid * deltarid;

  double dsh_rid;
  if (lrid1 == lrid0) {
    dsh_rid = 0.0;
  } else {
    dsh_rid = (shb[MT+lrid1][MT+lri0] - shb[MT+lrid0][MT+lri0])
        / deltaridta;
  }
  double dsh_ri;
  if (lri1 == lri0) {
    dsh_ri = 0.0;
  } else {
    dsh_ri = (shb[MT+lrid0][MT+lri1] - shb[MT+lrid0][MT+lri0])
        / deltarita;
  }

  sh = shb[MT+lrid0][MT+lri0] 
      + dsh_ri * deltari + dsh_rid * deltarid;
  
  double dss_rid;
  if (lrid1 == lrid0) {
    dss_rid = 0.0;
  } else {
    dss_rid = (ssb[MT+lrid1][MT+lri0] - ssb[MT+lrid0][MT+lri0])
        / deltaridta;
  }
  double dss_ri;
  if (lri1 == lri0) {
    dss_ri = 0.0;
  } else {
    dss_ri = (ssb[MT+lrid0][MT+lri1] - ssb[MT+lrid0][MT+lri0])
        / deltarita;
  }

  ss = ssb[MT+lrid0][MT+lri0] 
      + dss_ri * deltari + dss_rid * deltarid;

  return ;
}

inline double acosh1(const double &x) {
  return log(x + sqrt((x * x) - 1.0));
}
inline double wavelat(const double &xf, const double &yn) {
  return xf * acosh1(yn / xf);
}

static double eplatidepend_(const double &f, const double &an) {

  double f_30;
  double anum, den;
  double pi1, omega;
  const double an0 = 5.24e-3;

  pi1    = 4.0 * atan(1.0);
  omega  = pi1 / 43082.0e0;

  f_30   = omega;

  den    = wavelat(f_30, an0);

  anum   = wavelat(f, an);
  return anum / den;
}

#endif // LICOM_ENABLE_FORTRAN