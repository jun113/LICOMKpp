#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

#include "cmath"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

void thermal (const double (&tt)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
							const double (&ss)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
							const double (&pp)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
							      double (&aa)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
							      double (&bb)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
							const double (&mask)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]) {

	using CppDomain::nblocks_clinic;
  using CppPconstMod::od0;
  using CppPconstMod::vit;

  double tmp[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          tmp[iblock][k][j][i] = - pp[iblock][k][j][i] / od0 /
          10000.0 * mask[iblock][k][j][i];
        }
      }
    }

    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {

          bb[iblock][k][j][i] = (0.785567e-3 - 0.301985e-5 * tt[iblock][k][j][i] + 
              0.555579e-7 * pow(tt[iblock][k][j][i], 2) - 0.415613e-9 * pow(tt[iblock][k][j][i], 3) +
             (ss[iblock][k][j][i] * 1000.0) * (-0.356603e-6 + 0.788212e-8 *
              tt[iblock][k][j][i] + 0.408195e-10 * tmp[iblock][k][j][i] -  
              0.602281e-15 * pow(tmp[iblock][k][j][i], 2)) + pow(ss[iblock][k][j][i] * 1000.0, 2) * 
              (0.515032e-8) + tmp[iblock][k][j][i] * (-0.121555e-7 + 0.192867e-9 * tt[iblock][k][j][i] - 0.2131127e-11 * 
              pow(tt[iblock][k][j][i], 2)) + pow(tmp[iblock][k][j][i], 2) * (0.176621e-12 - 0.175379e-14 * 
              tt[iblock][k][j][i]) + pow(tmp[iblock][k][j][i], 3) * (0.12155e-17)) * 
              vit[iblock][k][j][i];

          aa[iblock][k][j][i] = (0.665157e-1 + 0.170907e-1 * tt[iblock][k][j][i] - 
              0.203814e-3 * pow(tt[iblock][k][j][i], 2) + 0.298357e-5 * 
              pow(tt[iblock][k][j][i], 3) - 0.255019e-7 * pow(tt[iblock][k][j][i], 4) +
              (ss[iblock][k][j][i] * 1000.0) * (0.378110e-2 - 0.846960e-4 * tt[iblock][k][j][i] - 
              0.164759e-6 * tmp[iblock][k][j][i] - 0.251520e-11 * pow(tmp[iblock][k][j][i], 2)) +
              pow(ss[iblock][k][j][i] * 1000.0, 2) * (-0.678662e-5) + 
              tmp[iblock][k][j][i] * (0.380374e-4 - 0.933746e-6 * tt[iblock][k][j][i] + 0.791325e-8 * 
              pow(tt[iblock][k][j][i], 2)) + 0.512857e-12 * pow(tmp[iblock][k][j][i], 2) * 
              pow(tt[iblock][k][j][i], 2) - 0.302285e-13 * pow(tmp[iblock][k][j][i], 3)) * 
              bb[iblock][k][j][i] * vit[iblock][k][j][i];
        }
      }
    }

	}
	return ;
}

#endif // LICOM_ENABLE_FORTRAN