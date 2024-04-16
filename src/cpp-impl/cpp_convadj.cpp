#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_output_mod.h"

#include <cstdio>
#include <cstdlib>
#include <string>

inline double dens(const double &, const double &, const int &);

void cpp_convadj(){
  using CppDomain::nblocks_clinic;
  using CppParamMod::NTRA;
  using CppParamMod::KM;  
  using CppParamMod::JMT;
  using CppParamMod::IMT;
  using CppParamMod::mytid;
  using CppPconstMod::to;
  using CppPconstMod::so;
  using CppPconstMod::dzp;
  using CppPconstMod::adv_tracer;
  using CppPconstMod::ist;
  using CppPconstMod::dts;
  using CppPconstMod::vit;
  using CppGrid::kmt;
  using CppTracerMod::at;
  using CppTracerMod::atb; 
  using CppTracerMod::dt_conv;
  using CppTracerMod::tend;
#ifdef LOWRES
  using CppOutputMod::icmon;
#endif // LOWRES
  double rhoup[KM], rholo[KM], trasum[2];
  double tup, sup, tlo, slo, dztsum, tramix;// ek0, ek1;
  int    kcon, lctot, lcven, l1, l, lcon, lcona, lconb, lmix;// n2;
  double c2dtts; 

  for (int k = 0; k < KM; ++k) {
    rhoup[k] = 0.0;
    rholo[k] = 0.0;
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0 ; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        kcon = kmt[iblock][j][i]; 
        lctot = 0;
        lcven = 0;
        if(kcon == 0){
          continue ;
        }//end if
        lcven =  1;
        lcon  = -1;
      
        for (l = 0 ; l < KM - 1 ; ++l) {
          l1 = l + 1;
          tup = at[iblock][0][l1][j][i] - to[l1];
          sup = at[iblock][1][l1][j][i] - so[l1];
          tlo = at[iblock][0][l][j][i]  - to[l1];
          slo = at[iblock][1][l][j][i]  - so[l1];
          rhoup[l1] = dens(tup, sup, l1);
          rholo[l]  = dens(tlo, slo, l1);
        }// end for
        for (int k = kcon - 1; k >= 1; --k) { //could be opt
          if (rholo[k - 1] > rhoup[k]) {
            lcon = k - 1;
          } 
        }//end for
        if (lcon < 0) {
          continue ;
        }
        for (;;) { // start CONV_1
          lcona = lcon;
          lconb = lcon + 1;
          dztsum = dzp[lcona] + dzp[lconb];
          for (int n = 0; n < 2; ++n) {
            trasum[n] = at[iblock][n][lcona][j][i] * dzp[lcona]
                      + at[iblock][n][lconb][j][i] * dzp[lconb];

            tramix = trasum[n] / dztsum;

            at[iblock][n][lcona][j][i] = tramix; 
            at[iblock][n][lconb][j][i] = tramix;
          } // end for
          for (;;) {//start CONV_2 
            if (lconb != (kcon - 1)) {
              l1 = lconb + 1;
              rholo[lconb] = dens(at[iblock][0][lconb][j][i] - to[l1],
                                  at[iblock][1][lconb][j][i] - so[l1],
                                      l1);

              if (rholo[lconb] > rhoup[l1]) {
                ++lconb;
                dztsum += dzp[lconb];
                for (int n = 0; n < 2; ++n) {
                  trasum[n] += at[iblock][n][lconb][j][i] * dzp[lconb];
                  tramix = trasum[n] / dztsum ; 
                  for (lmix = lcona; lmix <= lconb; ++lmix) {
                    at[iblock][n][lmix][j][i] = tramix;
                  }
                }
                continue;
              }// end if
            }//end if     
            if (lcona > 0) {
              l1 = lcona - 1;
              rholo[l1] = dens(at[iblock][0][l1][j][i] - to[lcona],
                               at[iblock][1][l1][j][i] - so[lcona],
                                   lcona);
              rhoup[lcona] = dens(at[iblock][0][lcona][j][i] - to[lcona],
                                  at[iblock][1][lcona][j][i] - so[lcona],
                                      lcona);

              if (rholo[lcona - 1] > rhoup[lcona]) {
                --lcona;
                dztsum += dzp[lcona];
                for (int n = 0; n < 2; ++n) {
                  trasum[n] += at[iblock][n][lcona][j][i] * dzp[lcona];
                  tramix = trasum[n] / dztsum ;
                  for (lmix = lcona; lmix <= lconb; ++lmix) {
                    at[iblock][n][lmix][j][i] = tramix ;
                  }
                }
                continue;
              }// end if
            }//end if
            break;
          }//end CONV_2
          lctot += lconb - lcona + 1;
          if (lcona == 0) {
            lcven = lconb - lcona + 1;
          }
          if (lconb == kcon - 1) {
#ifdef LOWRES 
            icmon[iblock][0][j][i] += lctot;
            icmon[iblock][1][j][i] += lcven;  
#endif // LOWRES 
            goto III;// cycle III
          }// end if
          lcon = lconb ;
          while(true) {//start conv_3 
            lcon = lcon + 1;
            if (lcon == kcon - 1) {
#ifdef LOWRES
              icmon[iblock][0][j][i] += lctot;
              icmon[iblock][0][j][i] += lctot;
              icmon[iblock][1][j][i] += lcven;
#endif // LOWRES
              goto III ;//cycle III 
            }
            if (lcon == (KM-1)) {
              break ;
            }
            if (rholo[lcon] <= rhoup[lcon + 1]) {
              continue ;//cycle CONV_3
            }
            break ; //exit CONV_3
          }// end CONV_3 
        }// end CONV_1
        III:;
      }// end III
    }// end JJJ
  }// end nblocks

  std::string str_adv_tracer(adv_tracer);
  if (str_adv_tracer.find("centered") != 
      str_adv_tracer.npos) {
    if (ist >= 1) {
      c2dtts = dts * 2.0;
    } else {
      c2dtts = dts ;
    }
  } else if (str_adv_tracer.find("tspas") !=
      str_adv_tracer.npos) {
    c2dtts = dts ;
  } else {
    if (mytid == 0) {
      printf("error in convadj\n");
      printf("The false advection option for tracer in convadj\n");
    }
    exit(0);
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int n = 0; n < NTRA; ++n) {
      for (int k = 0; k < KM; ++k) {  
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            dt_conv[iblock][n][k][j][i] = (at[iblock][n][k][j][i] -
                atb[iblock][n][k + 1][j][i]) / 
                    c2dtts * vit[iblock][k][j][i];

            tend[iblock][n][k][j][i] += dt_conv[iblock][n][k][j][i];
          }
        }
      }
    }
  }

  if (str_adv_tracer.find("tspas") !=
      str_adv_tracer.npos) {
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int n = 0; n < NTRA; ++n) {
        for (int k = 0; k < KM; ++k) {  
          for (int j = 0; j < JMT; ++j) {
            for (int i = 0; i < IMT; ++i) {
              atb[iblock][n][k + 1][j][i] = 
                  at[iblock][n][k][j][i];
            }
          }
        }
      }
    }//end for
  } 
  return ;
}// end convadj

inline double dens(const double &tq, const double &sq, const int &kk) {
  using CppPconstMod::c;
  double dens;
  dens = (c[0][kk] + (c[3][kk] + c[6][kk] * sq) * sq +
         (c[2][kk] +  c[7][kk] * sq + c[5][kk] * tq) * tq) * tq +
         (c[1][kk] + (c[4][kk] + c[8][kk] * sq) * sq) * sq;
  return dens;
}
#endif // LICOM_ENABLE_FORTRAN
