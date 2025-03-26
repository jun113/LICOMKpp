#include "../head/def-undef.h"
#include "../head/cpp_domain.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

void cpp_time_interplate_chlorophyll () {
	using CppForcMod::sss;
	using CppForcMod::sss3;
#if (defined SOLARCHLORO)
	using CppForcMod::chloro;
	using CppForcMod::chloro3;
#endif // SOLARCHLORO
  using CppPconstMod::mon0;
  using CppPconstMod::nmonth;
  using CppPconstMod::dd_licom;
  using CppDomain::nblocks_clinic;

  int ipt1, ipt2;
	double factor;
  if (dd_licom <= 15) {
    ipt1 = mon0 - 1;
  	if (ipt1 == 0) {
      ipt1 = 12;
    }
    ipt2 = mon0;
    factor = static_cast<double> (dd_licom - 15) / static_cast<double> (nmonth[ipt1-1]) + 1;
  } else {
  	ipt1 = mon0;
    ipt2 = mon0%12 + 1;
    factor = static_cast<double> (dd_licom - 15) / static_cast<double> (nmonth[ipt1-1]);
  }

	ipt1 -= 1;
	ipt2 -= 1;

	for (int iblock = 0; iblock < CppDomain::nblocks_clinic; ++iblock) {
  	for (int j = 0; j < CppParamMod::JMT; ++j) {
    	for (int i = 0; i < CppParamMod::IMT; ++i) {
				sss[iblock][j][i] = (sss3[iblock][ipt2][j][i] - sss3[iblock][ipt1][j][i])
						* factor + sss3[iblock][ipt1][j][i];
    	}
  	}
	}
#if (defined SOLARCHLORO)
	for (int iblock = 0; iblock < CppDomain::nblocks_clinic; ++iblock) {
  	for (int j = 0; j < CppParamMod::JMT; ++j) {
    	for (int i = 0; i < CppParamMod::IMT; ++i) {
        chloro[iblock][j][i] = (chloro3[iblock][ipt2][j][i] - chloro3[iblock][ipt1][j][i])
          * factor + chloro3[iblock][ipt1][j][i]
    	}
  	}
	}
  CALL  SW_ABSOR
#endif // SOLARCHLORO
	/*
	// Fortran code
  IF ( dd_licom <= 15) THEN
    IPT1 = MON0-1
    IF (IPT1 == 0 ) IPT1 = 12
    IPT2 = MON0
    FACTOR = FLOAT (dd_licom -15)/ FLOAT (NMONTH (IPT1)) + 1
  ELSE
    IPT1 = MON0
    IPT2 = MOD (MON0,12) +1
    FACTOR = FLOAT (dd_licom -15)/ FLOAT (NMONTH (IPT1))
  END IF
  !!$OMP PARALLEL DO PRIVATE (J,I)
  DO iblock = 1, nblocks_clinic
    DO J = 1,JMT
      DO I = 1,IMT
        SSS(I,J,iblock)=(SSS3(I,J,IPT2,iblock)-SSS3(I,J,IPT1,iblock))*FACTOR+SSS3(I,J,IPT1,iblock)
#if (defined SOLARCHLORO)
        chloro(I,J,iblock)= ( (chloro3 (I,J,IPT2,iblock) - chloro3 (I,J,IPT1,iblock))     &
          * FACTOR + chloro3 (I,J,IPT1,iblock))
#endif
      END DO
    END DO
  END DO

#if (defined SOLARCHLORO)
       CALL  SW_ABSOR
#endif
	*/

	return ;
}