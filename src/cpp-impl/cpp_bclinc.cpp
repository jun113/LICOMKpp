#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_constant_mod.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_work_mod.h"

#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

using CppConstantMod::G;
using CppConstantMod::C0;
using CppConstantMod::P5;
using CppConstantMod::P25;
using CppDynMod::sbcx;//M J I
using CppDynMod::sbcy;//M J I
using CppDynMod::bbcx;//M J I
using CppDynMod::bbcy;//M J I
using CppDynMod::up;//M K J I
using CppDynMod::vp;//M K J I
using CppDynMod::dlu;//K J I
using CppDynMod::dlv;//K J I
using CppDynMod::h0bf;//M J I
using CppDynMod::h0bl;//M J I
using CppDynMod::gg;//K J I
using CppDynMod::u;//M K J I
using CppDynMod::v;//M K J I
using CppDynMod::ub;//M J I
using CppDynMod::vb;//M J I
using CppDynMod::utf;//M K J I
using CppDynMod::vtf;//M K J I
using CppDomain::nblocks_clinic;
using CppForcMod::su;//M J I
using CppForcMod::sv;//M J I
using CppForcMod::psa;//M J I
using CppGrid::kmu;//M NY NX
using CppGrid::fcor;//M NY NX
using CppGrid::dyur; //M NY NX
using CppGrid::dxur; //M NY NX
using CppGrid::horiz_grid_opt;
using CppParamMod::KM;  
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JET;
using CppParamMod::JST;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppPconstMod::od0;
using CppPconstMod::c0f;
using CppPconstMod::cag;
using CppPconstMod::sag;
using CppPconstMod::isc;
using CppPconstMod::onbb;
using CppPconstMod::snlat;//M J I
using CppPconstMod::vit;//M k J I
using CppPconstMod::viv;//M K J I
using CppPconstMod::dzp;//K
using CppPconstMod::afc1;
using CppPconstMod::afc2;
using CppPconstMod::dtc;
using CppPconstMod::dtc2;
using CppPconstMod::epea;//M J I
using CppPconstMod::epeb;//M J I
using CppPconstMod::epla;//M J I
using CppPconstMod::eplb;//M J I
using CppPconstMod::odzp;//K
using CppPconstMod::odzt;//K
using CppPconstMod::ohbu;//M J I
using CppPconstMod::ohbt;//M J I
using CppPconstMod::akmu;//M K j I
using CppPconstMod::zkt;//K
using CppWorkMod::wka;//M K J I
using CppWorkMod::wkk;//K+1
using CppWorkMod::work;//M J I

// BCLINC
void cpp_bclinc(){

#ifdef  LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BCLINC
#define LICOM_ENABLE_TEST_BCLINC
#endif  // LICOM_ENABLE_TEST_TIME

#ifdef LICOM_ENABLE_TEST_BCLINC
    using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_BCLINC

#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC

  int kmb;
  // int errorcode;
  double aidif = 0.0;

#ifdef CANUTO
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock){
    for(int j = 1; j < JMT-1; ++j){
      for(int i = 1; i < IMT-1; ++i){
        kmb = kmu[iblock][j][i];
        if ( kmb >= 1) {
          sbcx[iblock][j][i] = su[iblock][j][i] * od0;
          sbcy[iblock][j][i] = sv[iblock][j][i] * od0; 

          bbcx[iblock][j][i] = c0f * 
              sqrt(up[iblock][kmb-1][j][i] * up[iblock][kmb-1][j][i] 
                 + vp[iblock][kmb-1][j][i] * vp[iblock][kmb-1][j][i]) 
                * (up[iblock][kmb-1][j][i] * cag + snlat[iblock][j][i] 
                 * vp[iblock][kmb-1][j][i] * sag); 
                  
          bbcy[iblock][j][i] = c0f * 
              sqrt(up[iblock][kmb-1][j][i] * up[iblock][kmb-1][j][i] 
                 + vp[iblock][kmb-1][j][i] * vp[iblock][kmb-1][j][i]) 
                * (vp[iblock][kmb-1][j][i] * cag - snlat[iblock][j][i] 
                 * up[iblock][kmb-1][j][i] * sag);
        }//end if
      }
    }
  }//end for
  aidif = 0.5; 
#endif // CANUTO


  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {  
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          dlu[iblock][k][j][i] -= fcor[iblock][j][i] 
              * vp[iblock][k][j][i];
          dlv[iblock][k][j][i] += fcor[iblock][j][i] 
              * up[iblock][k][j][i];
        }
      }   
    }
  }//end for 
  double aa = 0.0;
  if (isc != 0) {
    aa = 0.5;
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {  
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h0bf[iblock][j][i] *= onbb;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {  
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        work[iblock][j][i] = aa * h0bf[iblock][j][i] 
            + (1.0 - aa) * h0bl[iblock][j][i];
      }
    }
  }
 
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {  
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wkk[0] = (psa[iblock][j][i] * od0 + work[iblock][j][i] * G) 
            * vit[iblock][0][j][i];
        for (int k = 0; k < KM; ++k) { 
          wkk[k+1] = wkk[k] - gg[iblock][k][j][i] * 
              dzp[k] * vit[iblock][k][j][i];
        }
        for (int k = 0; k < KM; ++k) { 
          wka[iblock][k][j][i] = P5 * (wkk[k+1] + wkk[k]);
        }
      }
    }
  }//end for bclinc_4

  double  gradx[JMT][IMT];
  double  grady[JMT][IMT];

  for (int j = 0; j < JMT; ++j) {
    for (int i = 0; i < IMT; ++i) {
      gradx[j][i] = 0.0;
      grady[j][i] = 0.0;
    }
  }

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      grad(k, gradx, grady, wka[iblock][k]); 
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          dlv[iblock][k][j][i] -= grady[j][i];
          dlu[iblock][k][j][i] -= gradx[j][i]; 
        }
      }   
    }
  }//end for 5

   
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wka[iblock][k][j][i] = (1.0 + ohbt[iblock][j][i] * zkt[k]) 
              * work[iblock][j][i] * vit[iblock][k][j][i];
        }
      }   
    }
  }//end for 6

  double ggu;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) { 
      grad(k, gradx, grady, wka[iblock][k]); 
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          ggu = P25 * (gg[iblock][k][j  ][i] + gg[iblock][k][j  ][i-1]
                     + gg[iblock][k][j+1][i] + gg[iblock][k][j+1][i-1]);

          dlv[iblock][k][j][i] += ggu * grady[j][i];
          dlu[iblock][k][j][i] += ggu * gradx[j][i];
        }
      }   
    }
  }//end for 7
  double  wk1, wk2;
  
  if (isc == 0) {
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 1; j < JMT-1; ++j) {
          for (int i = 1; i < IMT-1; ++i) {
            wk1 = epea[iblock][j][i] * dlv[iblock][k][j][i] 
                + epeb[iblock][j][i] * dlu[iblock][k][j][i];
            wk2 = epea[iblock][j][i] * dlu[iblock][k][j][i] 
                - epeb[iblock][j][i] * dlv[iblock][k][j][i];

            dlv[iblock][k][j][i] = wk1 * viv[iblock][k][j][i];
            dlu[iblock][k][j][i] = wk2 * viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 8
  } else {
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 1; j < JMT-1; ++j) {
          for (int i = 1; i < IMT-1; ++i) {
            wk1 = epla[iblock][j][i] * dlv[iblock][k][j][i] 
                + eplb[iblock][j][i] * dlu[iblock][k][j][i];
            wk2 = epla[iblock][j][i] * dlu[iblock][k][j][i] 
                - eplb[iblock][j][i] * dlv[iblock][k][j][i];

            dlv[iblock][k][j][i] = wk1 * viv[iblock][k][j][i];
            dlu[iblock][k][j][i] = wk2 * viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 9
  }

  /*
  double fil_lat1 = 63.0;
  double fil_lat2 = 63.0;
  if (horiz_grid_opt == "lat_lon") {
    call POP_HaloUpdate(DLU, POP_haloClinic, POP_gridHorzLocSWcorner , &
        POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

    call POP_HaloUpdate(DLV, POP_haloClinic, POP_gridHorzLocSWcorner , &
        POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

    smuv_3d(dlv, viv, fil_lat1);
    smuv_3d(dlv, viv, fil_lat1);
  }
  */
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
 
  if (isc < 1) {
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 2; j < JMT - 2; ++j) {
          for (int i = 2; i < IMT - 2; ++i) {
            v[iblock][k][j][i] = vp[iblock][k][j][i] 
                + dlv[iblock][k][j][i] * dtc;
            u[iblock][k][j][i] = up[iblock][k][j][i] 
                + dlu[iblock][k][j][i] * dtc;
          }
        }   
      }
    } //end for 10

    invtriu(u, sbcx, bbcx, akmu, aidif, dtc);
    invtriu(v, sbcy, bbcy, akmu, aidif, dtc);   
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
    my_time.testTime_start("bclinc halo u v");
#endif // LICOM_ENABLE_TEST_BCLINC
    // pop_haloupdate_bclinc2_(&errorcode);
    CppPOPHaloMod::pop_halo_update_3dr8(u[0],
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
    CppPOPHaloMod::pop_halo_update_3dr8(v[0],
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc halo u v");
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC

    vinteg(u, work);

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            u[iblock][k][j][i] = (u[iblock][k][j][i] - 
                work[iblock][j][i] + ub[iblock][j][i]) * 
                    viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 11

    vinteg(v, work); 

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            v[iblock][k][j][i] = (v[iblock][k][j][i] - 
                work[iblock][j][i] + vb[iblock][j][i]) * 
                    viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 12
    ++isc;
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
  } else {
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 2; j < JMT - 2; ++j) {
          for (int i = 2; i < IMT - 2; ++i) {
            wka[iblock][k][j][i] = vp[iblock][k][j][i] 
                + dlv[iblock][k][j][i] * dtc2;
          }
        }   
      }
    }//end for 13

    invtriu(wka, sbcy, bbcy, akmu, aidif, dtc2);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
    my_time.testTime_start("bclinc halo wka");
#endif // LICOM_ENABLE_TEST_BCLINC
    // pop_haloupdate_bclinc3_(&errorcode);
    CppPOPHaloMod::pop_halo_update_3dr8(wka[0],
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc halo wka");
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
    vinteg(wka, work);

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {  
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            wka[iblock][k][j][i] = (wka[iblock][k][j][i] 
                - work[iblock][j][i] + vb[iblock][j][i]) 
                    * viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 14

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {  
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            vp[iblock][k][j][i] = afc2 * v[iblock][k][j][i] 
                + afc1 * (vp[iblock][k][j][i] + wka[iblock][k][j][i]);

            v[iblock][k][j][i] = wka[iblock][k][j][i];
          }
        }   
      }
    }//end for 15

    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {  
        for (int j = 2; j < JMT - 2; ++j) {
          for (int i = 2; i < IMT - 2; ++i) {
            wka[iblock][k][j][i] = up[iblock][k][j][i] 
                + dlu[iblock][k][j][i] * dtc2;
          }
        }   
      }
    }//end for 16
   
    invtriu(wka, sbcx, bbcx, akmu, aidif, dtc2);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
    my_time.testTime_start("bclinc halo wka");
#endif // LICOM_ENABLE_TEST_BCLINC
    // pop_haloupdate_bclinc3_(&errorcode);
    CppPOPHaloMod::pop_halo_update_3dr8(wka[0],
        CppDomain::POP_haloClinic_C, 
        CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
        CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc halo wka");
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
    vinteg(wka, work);
 
    for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
      for (int k = 0; k < KM; ++k) {
        for (int j = 0; j < JMT; ++j) {
          for (int i = 0; i < IMT; ++i) {
            wka[iblock][k][j][i] = (wka[iblock][k][j][i] 
                - work[iblock][j][i] + ub[iblock][j][i]) 
                    * viv[iblock][k][j][i];
          }
        }   
      }
    }//end for 17

    for(int iblock = 0; iblock < nblocks_clinic; ++iblock){
      for(int k = 0; k < KM; ++k){  
        for(int j = 0; j < JMT; ++j){
          for(int i = 0; i < IMT; ++i){
            up[iblock][k][j][i] = afc2 * u[iblock][k][j][i] 
                + afc1 * (up[iblock][k][j][i] + wka[iblock][k][j][i]);

            u[iblock][k][j][i] = wka[iblock][k][j][i];
          }
        }   
      }
    }//end for 18

    /*
    if (trim(horiz_grid_opt) == 'lat_lon') then
      IF (MOD(ISC,160)==1) THEN
        CALL SMUV_3D (U,VIV,fil_lat2)
        CALL SMUV_3D (V,VIV,fil_lat2)
        CALL SMUV_3D (UP,VIV,fil_lat2)
        CALL SMUV_3D (VP,VIV,fil_lat2)
      END IF
    end if
    */
    
    ++isc;
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
  } //end if else
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_start("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {  
      for (int j = JST - 1; j < JET; ++j) {
        for (int i = 0; i < IMT; ++i) {
          utf[iblock][k][j][i] += u[iblock][k][j][i];
          vtf[iblock][k][j][i] += v[iblock][k][j][i];
        }
      }   
    }
  }//end for 19

  //check part
  fortran_mpi_barrier_();
  //mpi_barrier
#ifdef LICOM_ENABLE_TEST_BCLINC
    my_time.testTime_stop("bclinc calc");
#endif // LICOM_ENABLE_TEST_BCLINC
  return ;
}// end cpp_bclinc


#endif // LICOM_ENABLE_FORTRAN
