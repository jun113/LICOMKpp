#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_msg_mod.h"
#include "../head/cpp_output_mod.h"
#include "../head/cpp_tracer_mod.h"

#include "../head/fortran_extern_functions.h"

#include "netcdf.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

template<typename T>
static inline T sign (const T &x, const T &y) {
  return y >= static_cast<T>(0.0) ? std::abs(x) : -std::abs(x);
}

static void ncar_ocean_fluxes_jra(const double *u_del, const double *t, 
		const double *ts, const double *q, const double *qs, const double *z, 
				const double *avail, double *sh, double *lh, double *tau, double *ustar) {
	// !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	// ! Over-ocean fluxes following Large and Yeager (used in NCAR models) 
	// ! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
	// !
	// ! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
	// ! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
	// !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	using CppParamMod::IMT;
	using CppParamMod::JMT;
	double cd[JMT][IMT];
	double ch[JMT][IMT];
	double ce[JMT][IMT];
	// memset(&(cd[0][0]), 0, sizeof(double) * IMT * JMT);
	// memset(&(ch[0][0]), 0, sizeof(double) * IMT * JMT);
	// memset(&(ce[0][0]), 0, sizeof(double) * IMT * JMT);
	// double bstar[JMT][IMT];
  const double grav(9.80), vonkarm(0.40);
	const double reciprocal_vonkarm = 2.5;

	double cd_n10    = 0.0;
	double cd_n10_rt = 0.0;
	double ce_n10    = 0.0;
	double stab      = 0.0;
	double ch_n10    = 0.0;

	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			if (avail[j * IMT + i] > 0.5) {
				const double tv = t[j * IMT + i] * (1.0 + 0.608 * q[j * IMT + i]);
				const double u  = std::max(u_del[j * IMT + i], 0.5);
				double u10      = u;

				cd_n10    = (2.7/u10 + 0.142 + 0.0764*u10) * 0.001;
				cd_n10_rt = sqrt(cd_n10);
				ce_n10    = 34.6 * cd_n10_rt * 0.001;
				stab      = 0.5 + sign(0.5, t[j * IMT + i] - ts[j * IMT + i]);
				ch_n10    = (18.0 * stab + 32.7 * (1 - stab)) * cd_n10_rt * 0.001;

				cd[j][i] = cd_n10;
				ch[j][i] = ch_n10;
				ce[j][i] = ce_n10;
				double bstar;
				for (int jj = 0; jj < 2; ++jj) {
					// for: 1, n_itts. n_itts = 2
					// Monin-Obukhov iteration
					const double cd_rt = sqrt(cd[j][i]);
					ustar[j * IMT + i] = cd_rt * u;
					const double tstar = (ch[j][i] / cd_rt) * (t[j*IMT + i] - ts[j*IMT + i]);
					const double qstar = (ce[j][i] / cd_rt) * (q[j*IMT + i] - qs[j*IMT + i]);
					// bstar[j][i] = grav * (tstar/tv + qstar/(q[j*IMT + i] + 1/0.608));
					bstar = grav * (tstar/tv + qstar/(q[j*IMT + i] + 1.0/0.608));

					// double zeta = vonkarm * bstar[j][i] * z[j*IMT + i] / 
					// 		(ustar[j*IMT + i] * ustar[j*IMT + i]);
					double zeta = vonkarm * bstar * z[j*IMT + i] / 
							(ustar[j*IMT + i] * ustar[j*IMT + i]);
					zeta = sign(std::min(std::abs(zeta), 10.0), zeta);

					double x2 = sqrt(std::abs(1.0 - 16.0 * zeta));
					x2 = std::max(x2, 1.0);
					const double x = sqrt(x2);
					double psi_m;
					double psi_h;
					if (zeta > 0) {
						psi_m = - 5.0 * zeta;
						psi_h = psi_m;
					} else {
						psi_m = log((1.0 + 2.0*x + x2) * (1.0+x2)*0.125) - 2.0 * (atan(x) - atan(1.0));
						psi_h = 2.0 * log((1.0 + x2) * 0.5);
					}
					u10 = u / (1.0 + cd_n10_rt * (log(z[j*IMT + i] * 0.1) - psi_m) * reciprocal_vonkarm);

					cd_n10    = (2.7/u10 + 0.142 + 0.0764*u10) * 0.001;
					cd_n10_rt = sqrt(cd_n10);
					ce_n10    = 34.6 * cd_n10_rt * 0.001;
					stab      = 0.5 + sign(0.5, zeta);
					ch_n10    = (18.0*stab + 32.7*(1-stab)) * cd_n10_rt * 0.001;
					// diagnostic
					// z0        = 10 * exp(-vonkarm / cd_n10_rt);
					double xx = (log(z[j*IMT + i] * 0.1) - psi_m) * reciprocal_vonkarm;
					const double tmp = (1.0 + cd_n10_rt * xx);
					cd[j][i]  = cd_n10 / (tmp * tmp);
					xx = (log(z[j * IMT + i] * 0.1) - psi_h) * reciprocal_vonkarm;

					ch[j][i] = ch_n10 / (1.0 + ch_n10 * xx / cd_n10_rt) * sqrt(cd[j][i] / cd_n10);
					ce[j][i] = ce_n10 / (1.0 + ce_n10 * xx / cd_n10_rt) * sqrt(cd[j][i] / cd_n10);
				}
			} else {
				cd[j][i] = 0.0;
				ch[j][i] = 0.0;
				ce[j][i] = 0.0;
			}
		}
	}

	const double L  = 2.5e6;
	const double R0 = 1.22;
	const double CP = 1000.5;
	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			sh[j*IMT + i]  = R0 * CP * ch[j][i] * 
					(t[j*IMT + i] - ts[j*IMT + i]) * u_del[j*IMT + i];
			lh[j*IMT + i]  = R0 * ce[j][i] * L * 
					(q[j*IMT + i] - qs[j*IMT + i]) * u_del[j*IMT + i];
			tau[j*IMT + i] = R0 * cd[j][i] * u_del[j*IMT+ i];
		}
	}

	return ;
}

void read_jra (const int &start_day, const int &num_day, 
		const char* fname, float* const val) {
	using CppParamMod::S_JMT;
	using CppParamMod::S_IMT;

	const size_t start[3] = {(size_t)(start_day-1), (size_t)(0), 		 (size_t)(0)};
	const size_t count[3] = {(size_t)(num_day),     (size_t)(S_JMT), (size_t)(S_IMT)};

	int ncid;
	// int nc_status;
	nc_open(fname, NC_NOWRITE, &ncid);
	// check error
	nc_get_vara_float(ncid, 4, &(start[0]), &(count[0]), val);
	// check error
	nc_close(ncid);
	// check error

	return ;
}

void read_jra (const char* fname, double* const buffer, 
    double* const s_lon, double* const s_lat) {
  using CppParamMod::S_JMT;
  using CppParamMod::S_IMT;

	int ncid;
	// int nc_status;
	nc_open(fname, NC_NOWRITE, &ncid);
	// check error
	nc_get_var_double(ncid, 2, buffer);
	// check error
	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_lon[j*S_IMT + i] = buffer[i];
		}
	}
	nc_get_var_double(ncid, 3, buffer);
	// check error
	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_lat[j*S_IMT + i] = buffer[j];
		}
	}
	// check error
	nc_close(ncid);
	// check error
	return ;
}

void read_jra (const int &start_day, const int &num_day, const char* fname,
		double* const buffer, double* const s_lon, double* const s_lat, float* const val) {
	using CppParamMod::S_JMT;
	using CppParamMod::S_IMT;

	const size_t start[3] = {
			(size_t)(start_day-1), (size_t)(0), 		(size_t)(0)};
	const size_t count[3] = {
			(size_t)(num_day), 	   (size_t)(S_JMT), (size_t)(S_IMT)};

	int ncid;
	// int nc_status;
	nc_open(fname, NC_NOWRITE, &ncid);
	// check error
	nc_get_var_double(ncid, 2, buffer);
	// check error
	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_lon[j*S_IMT + i] = buffer[i];
		}
	}
	nc_get_var_double(ncid, 3, buffer);
	// check error
	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_lat[j*S_IMT + i] = buffer[j];
		}
	}
	nc_get_vara_float(ncid, 4, &(start[0]), &(count[0]), val);
	// check error
	nc_close(ncid);
	// check error
	return ;
}

static void near(const double *a, const int &mx, const int &my, 
		const double *alon, const double *alat, double* const b) {
	// ! input : a
	// ! output: b
	// ! A bi-linear interpolation will be used for the initial guess, then
	// ! a refilling procedure be called to redefine the missing data.
	// ! mx,my            x and y grids number of a
	// ! alon,alat        lontitude and latitude of a
	// ! spval            missing flag
	// ! nf               1 for T grid; 0 for U grid
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	using CppParamMod::IMT_GLOBAL;
	using CppParamMod::JMT_GLOBAL;
	using CppPconstMod::lon_o;
	using CppPconstMod::lat_o;

	for (int j = 0; j < JMT_GLOBAL; ++j) {
		for (int i = 0; i < IMT_GLOBAL; ++i) {
			b[j*IMT_GLOBAL + i] = CppOutputMod::SPVAL;
		}
	}

	double lon_tmp;
	for (int j = 0; j < JMT_GLOBAL; ++j) {
		for (int i = 0; i < IMT_GLOBAL; ++i) {
			// ic = static_cast<int>(isp);
			// jc = static_cast<int>(isp);
			int ic = INT32_MAX;
			int jc = INT32_MAX;

			// ! ---find adjacent two grids on x direction (ic)
			if (lon_o[j][i] < 0.0) {
				lon_tmp = static_cast<double>(lon_o[j][i]) + 360.0;
			} else{
				lon_tmp = static_cast<double>(lon_o[j][i]);
			}
			for (int ip = 1; ip < mx; ++ip) {
				if ((alon[ip-1] <= lon_tmp) && (alon[ip] >= lon_tmp)) {
					ic = ip;
					break;
				}
			}
			// ---find adjacent two grids on y direction (jc)
			for (int jp = 1; jp < my; ++jp) {
				if ((alat[jp * mx] <= lat_o[j][i]) && 
						(alat[(jp-1) * mx] >= lat_o[j][i])) {
					jc = jp;
					break;
				}
			}

			if ((ic == INT32_MAX) || (jc == INT32_MAX)) {
				if (CppParamMod::mytid == 0) {
					printf("break adjacent grids has no found\n");
				}
				exit(EXIT_FAILURE);
			}
			// Bilinear interpolater
			b[j * IMT_GLOBAL + i] = a[jc * mx + ic];
		}
	}
	for (int j = 0; j < JMT_GLOBAL; ++j) {
		b[j*IMT_GLOBAL + IMT_GLOBAL-2] = b[j*IMT_GLOBAL    ];
		b[j*IMT_GLOBAL + IMT_GLOBAL-1] = b[j*IMT_GLOBAL + 1];
	}
	return;
}

static void interplation_nearest(const float *source, double* const object,
		const double *s_lon, const double *s_lat) {
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;

	const int iwk = S_IMT + 2;
	const int jwk = S_JMT + 2;

	double s_wx[jwk][iwk];
	double s_wy[jwk][iwk];
	double s_work[jwk][iwk];

	// both dx and dy are 0.5625
	const double dx_dy = 0.5625;

	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_wx[j+1][i+1] = s_lon[j*S_IMT + i];
		}
	}

	for (int j = 0; j < S_JMT; ++j) {
		s_wx[j+1][0    ] = s_wx[j+1][1    ] - dx_dy;
		s_wx[j+1][iwk-1] = s_wx[j+1][iwk-2] + dx_dy;
	}

	for (int i = 0; i < S_IMT; ++i) {
		s_wx[0    ][i+1] = s_wx[1][i+1];
		s_wx[jwk-1][i+1] = s_wx[0][i+1];
	}

	s_wx[0    ][0    ] = s_wx[1    ][0    ];
	s_wx[0    ][iwk-1] = s_wx[1    ][iwk-1];
	s_wx[jwk-1][0    ] = s_wx[jwk-2][0    ];
	s_wx[jwk-1][iwk-1] = s_wx[jwk-2][iwk-1];

	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_wy[j+1][i+1] = s_lat[j*S_IMT + i];
		}
	}
	for (int j = 0; j < S_JMT; ++j) {
		s_wy[j+1][0    ] = s_wy[j+1][1    ];
		s_wy[j+1][iwk-1] = s_wy[j+1][iwk-2];
	}
	s_wy[0][0    ] = s_wy[1][0    ] + dx_dy;
	s_wy[0][iwk-1] = s_wy[1][iwk-1] + dx_dy;
	for (int i = 0; i < S_IMT; ++i) {
		s_wy[0    ][i+1] = s_wy[1    ][i+1] + dx_dy;
		s_wy[jwk-1][i+1] = s_wy[jwk-2][i+1] - dx_dy;
	}
	s_wy[jwk-1][0    ] = s_wy[jwk-2][0    ] - dx_dy;
	s_wy[jwk-1][iwk-1] = s_wy[jwk-2][iwk-1] - dx_dy;

	for (int j = 0; j < S_JMT; ++j) {
		for (int i = 0; i < S_IMT; ++i) {
			s_work[j+1][i+1] = source[j*S_IMT + i];
		}
	}
	for (int j = 0; j < S_JMT; ++j) {
		s_work[j+1][0    ] = s_work[j+1][1    ];
		s_work[j+1][iwk-1] = s_work[j+1][iwk-2];
	}
	for (int i = 0; i < iwk; ++i) {
		s_work[0    ][i] = s_work[1    ][i];
		s_work[jwk-1][i] = s_work[jwk-2][i];
	}

	near(&(s_work[0][0]), iwk, jwk, &(s_wx[0][0]), &(s_wy[0][0]), object);

	return ;
}

static inline void interplationNearestThenScatterGlobal(const int &tid, 
		const double *s_lon, const double *s_lat, const float *input, 
				double *output, double *buffer) {
	if (CppParamMod::mytid == tid) {
		interplation_nearest(input, buffer, s_lon, s_lat);
	}
	int tmp_tid = tid;
	scatter_global_jra_r8_(output, buffer, &tmp_tid);
	pop_halo_update(output, CppParamMod::IMT, CppDomain::POP_haloClinic_C,
			CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
			CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
	return ;
}
static inline void interplationNearestThenScatterGlobal(const int &tid, const int &irec,
		const double *s_lon, const double *s_lat, const float *input, 
				double *output, double *buffer) {
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	if (CppParamMod::mytid == tid) {
		const int stride = (irec - 1) * S_JMT * S_IMT;
		interplation_nearest(&input[stride], buffer, s_lon, s_lat);
	}
	int tmp_tid = tid;
	scatter_global_jra_r8_(output, buffer, &tmp_tid);
	pop_halo_update(output, CppParamMod::IMT, CppDomain::POP_haloClinic_C,
			CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
			CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
	return ;
}

void cpp_init_jra_daily_low () {
	using CppParamMod::IMT;
	using CppParamMod::JMT;
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	using CppParamMod::IMT_GLOBAL;
	using CppParamMod::JMT_GLOBAL;
	using CppParamMod::MAX_BLOCKS_CLINIC;

	using CppPconstMod::s_lon;
	using CppPconstMod::s_lat;

	using CppForcMod::t10;
	using CppForcMod::u10;
	using CppForcMod::v10;
	using CppForcMod::slp;
	using CppForcMod::q10;
	using CppForcMod::swhf;
	using CppForcMod::lwhf;
	using CppForcMod::precr;
	using CppForcMod::precs;
	using CppForcMod::rf;
	using CppForcMod::si;

	using CppForcMod::buffer;

	using CppForcMod::tsa3;
	using CppForcMod::wspdu3;
	using CppForcMod::wspdv3;
	using CppForcMod::psa3;
	using CppForcMod::qar3;
	using CppForcMod::swv3;
	using CppForcMod::lwv3;
	using CppForcMod::rain3;
	using CppForcMod::snow3;
	using CppForcMod::runoff3;
	using CppForcMod::seaice3;
	using CppForcMod::wspd3;

	t10   = new float[S_IMT * S_JMT];
	u10   = new float[S_IMT * S_JMT];
	v10   = new float[S_IMT * S_JMT];
	slp   = new float[S_IMT * S_JMT];
	q10   = new float[S_IMT * S_JMT];
	swhf  = new float[S_IMT * S_JMT];
	lwhf  = new float[S_IMT * S_JMT];
	precr = new float[S_IMT * S_JMT];
	precs = new float[S_IMT * S_JMT];
	rf    = new float[S_IMT * S_JMT];
	si    = new float[S_IMT * S_JMT];

	s_lon = new double[S_IMT * S_JMT];
	s_lat = new double[S_IMT * S_JMT];

	buffer = new double[std::max(S_IMT, std::max(S_JMT, 
			IMT_GLOBAL * JMT_GLOBAL))];

	tsa3    = new double [JMT * IMT];
	wspdu3  = new double [JMT * IMT];
	wspdv3  = new double [JMT * IMT];
	psa3    = new double [JMT * IMT];
	qar3    = new double [JMT * IMT];
	swv3    = new double [JMT * IMT];
	lwv3    = new double [JMT * IMT];
	rain3   = new double [JMT * IMT];
	snow3   = new double [JMT * IMT];
	runoff3 = new double [JMT * IMT];
	seaice3 = new double [JMT * IMT];
	wspd3   = new double [JMT * IMT];

	return ;
}
void cpp_init_jra_daily_high () {
	using CppParamMod::IMT;
	using CppParamMod::JMT;
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	using CppParamMod::IMT_GLOBAL;
	using CppParamMod::JMT_GLOBAL;
	using CppParamMod::MAX_BLOCKS_CLINIC;

	using CppPconstMod::s_lon;
	using CppPconstMod::s_lat;

	if (CppMsgMod::nproc < 11) {
		if (CppParamMod::mytid == 0) {
			printf("The number of processes must be greater than 11 at high-resolution, currently %d !\n", 
					CppMsgMod::nproc);
		}
		exit(EXIT_FAILURE);
	}

	s_lon  = new double[S_IMT * S_JMT];
	s_lat  = new double[S_IMT * S_JMT];
	CppForcMod::buffer = new double[std::max(S_IMT, std::max(S_JMT, 
			IMT_GLOBAL * JMT_GLOBAL))];

	const int irec = CppPconstMod::number_day;
	if (CppParamMod::mytid == 0) {
		printf("jra daily start irec = %d\n", irec);
	}
	
	// ! read in jra data
	switch(CppParamMod::mytid) {
		case 0:
			CppForcMod::rf = new float[S_IMT * S_JMT * 365];
			// update s_lon s_lat
			read_jra(irec, 365, "runoff_all.2016_daily.nc", 
					CppForcMod::buffer, s_lon, s_lat, CppForcMod::rf);
			break;
		case 1:
			CppForcMod::t10 = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "t_10.2016_daily.nc", CppForcMod::t10);
			break;
		case 2:
			CppForcMod::u10 = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "u_10.2016_daily.nc", CppForcMod::u10);
			break;
		case 3:
			CppForcMod::v10 = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "v_10.2016_daily.nc", CppForcMod::v10);
			break;
		case 4:
			CppForcMod::slp = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "slp.2016_daily.nc",  CppForcMod::slp);
			break;
		case 5:
			CppForcMod::q10 = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "q_10.2016_daily.nc", CppForcMod::q10);
			break;
		case 6:
			CppForcMod::swhf = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "rsds.2016_daily.nc", CppForcMod::swhf);
			break;
		case 7:
			CppForcMod::lwhf = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "rlds.2016_daily.nc", CppForcMod::lwhf);
			break;
		case 8:
			CppForcMod::precr = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "rain.2016_daily.nc", CppForcMod::precr);
			break;
		case 9:
			CppForcMod::precs = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "snow.2016_daily.nc", CppForcMod::precs);
			break;
		case 10:
			CppForcMod::si = new float[S_IMT * S_JMT * 365];
			read_jra(irec, 365, "ice.2016_daily.nc", CppForcMod::si);
			break;
		default: break;
	}

	// Broadcasting the s_lon and s_lat
	MPI_Request request_lon, request_lat;
	MPI_Ibcast(s_lon, S_IMT * S_JMT, MPI_DOUBLE, 0,
			CppDomain::POP_haloClinic_C.communicator, &request_lon);
	MPI_Ibcast(s_lat, S_IMT * S_JMT, MPI_DOUBLE, 0,
			CppDomain::POP_haloClinic_C.communicator, &request_lat);

	CppForcMod::tsa3    = new double [JMT * IMT];
	CppForcMod::wspdu3  = new double [JMT * IMT];
	CppForcMod::wspdv3  = new double [JMT * IMT];
	CppForcMod::psa3    = new double [JMT * IMT];
	CppForcMod::qar3    = new double [JMT * IMT];
	CppForcMod::swv3    = new double [JMT * IMT];
	CppForcMod::lwv3    = new double [JMT * IMT];
	CppForcMod::rain3   = new double [JMT * IMT];
	CppForcMod::snow3   = new double [JMT * IMT];
	CppForcMod::runoff3 = new double [JMT * IMT];
	CppForcMod::seaice3 = new double [JMT * IMT];
	CppForcMod::wspd3   = new double [JMT * IMT];

	MPI_Wait(&request_lon, MPI_STATUS_IGNORE);
	MPI_Wait(&request_lat, MPI_STATUS_IGNORE);
	return ;
}

void cpp_jra_daily_low(const int &iday) {
	using CppDynMod::u;
	using CppDynMod::v;

	using CppForcMod::t10;
	using CppForcMod::u10;
	using CppForcMod::v10;
	using CppForcMod::slp;
	using CppForcMod::q10;
	using CppForcMod::swhf;
	using CppForcMod::lwhf;
	using CppForcMod::precr;
	using CppForcMod::precs;
	using CppForcMod::rf;
	using CppForcMod::si;

	using CppForcMod::buffer;

	using CppForcMod::tsa3;
	using CppForcMod::wspdu3;
	using CppForcMod::wspdv3;
	using CppForcMod::psa3;
	using CppForcMod::qar3;
	using CppForcMod::swv3;
	using CppForcMod::lwv3;
	using CppForcMod::rain3;
	using CppForcMod::snow3;
	using CppForcMod::runoff3;
	using CppForcMod::seaice3;
	using CppForcMod::wspd3;

	using CppForcMod::su;
	using CppForcMod::sv;
	using CppForcMod::lwv;
	using CppForcMod::swv;
	using CppForcMod::nswv;
	using CppForcMod::fresh;
	using CppForcMod::sshf;
	using CppForcMod::lthf;
	using CppForcMod::ustar;
	using CppForcMod::runoff;
	using CppForcMod::seaice;
	using CppParamMod::IMT;
	using CppParamMod::JMT;
	using CppParamMod::IMM;
	using CppParamMod::JST;
	using CppParamMod::JSM;
	using CppParamMod::JEM;
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	using CppParamMod::IMT_GLOBAL;
	using CppParamMod::JMT_GLOBAL;
	using CppParamMod::MAX_BLOCKS_CLINIC;
	using CppPconstMod::vit;
	using CppPconstMod::viv;
	using CppPconstMod::s_lon;
	using CppPconstMod::s_lat;
	using CppTracerMod::at;
	using CppPOPHaloMod::pop_halo_update_2dr8;

	int irec = iday;

	if (CppParamMod::mytid == CppParamMod::master_task) {
		printf("jra daily, irec = %d\n", irec);

		// ! read in jra data
		read_jra(irec, 1, "t_10.2016_daily.nc", t10);
		read_jra(irec, 1, "u_10.2016_daily.nc", u10);
		read_jra(irec, 1, "v_10.2016_daily.nc", v10);
		read_jra(irec, 1, "slp.2016_daily.nc",  slp);
		read_jra(irec, 1, "q_10.2016_daily.nc", q10);
		read_jra(irec, 1, "rsds.2016_daily.nc", swhf);
		read_jra(irec, 1, "rlds.2016_daily.nc", lwhf);
		read_jra(irec, 1, "rain.2016_daily.nc", precr);
		read_jra(irec, 1, "snow.2016_daily.nc", precs);
		read_jra(irec, 1, "ice.2016_daily.nc",  si);
		// update s_lon s_lat
		read_jra(irec, 1, "runoff_all.2016_daily.nc", 
				buffer, s_lon, s_lat, rf);
	}

	// ! interplate to T grid
	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, t10, tsa3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, u10, wspdu3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, v10, wspdv3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, slp, psa3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, q10, qar3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, swhf, swv3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, lwhf, lwv3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, precr, rain3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, precs, snow3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, rf, runoff3, buffer);

	interplationNearestThenScatterGlobal(CppParamMod::master_task, 
			s_lon, s_lat, si, seaice3, buffer);
	
	double uu[JMT][IMT];
	double vv[JMT][IMT];
	const double epsln = 1e-25;
	for (int j = JSM - 1; j < JEM; ++j) {
		for (int i = 1; i < IMM; ++i) {
			const double tmp = vit[0][0][j][i] / (
					viv[0][0][j  ][i] + viv[0][0][j  ][i+1]
				+ viv[0][0][j-1][i] + viv[0][0][j-1][i+1] + epsln);
			uu[j][i] = tmp * (
					u[0][0][j  ][i] + u[0][0][j  ][i+1] 
				+ u[0][0][j-1][i] + u[0][0][j-1][i+1]);
			vv[j][i] = tmp * (
					v[0][0][j  ][i] + v[0][0][j  ][i+1] 
				+ v[0][0][j-1][i] + v[0][0][j-1][i+1]);
		}
	}
	double windx[JMT][IMT];
	double windy[JMT][IMT];
	double model_sst[JMT][IMT];
	double zz[JMT][IMT];
	double qs[JMT][IMT];
	double theta[JMT][IMT];
	const double tok = 273.15;
	// ! transfer core data to what the subroutine need
	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			// ! relative speed to surface currents
			windx[j][i] = (wspdu3[j*IMT + i] - uu[j][i]) * vit[0][0][j][i];
			windy[j][i] = (wspdv3[j*IMT + i] + vv[j][i]) * vit[0][0][j][i];

			// !  1.0 is from mom4
			wspd3[j*IMT + i] = sqrt(windx[j][i] * windx[j][i] 
													+ windy[j][i] * windy[j][i] + 1.0) * vit[0][0][j][i];

			// ! using a transient temperature, not daily mean
			model_sst[j][i] = (at[0][0][0][j][i] + tok) * vit[0][0][j][i];
			zz[j][i] = 10.0;
			qs[j][i] = (0.98 * 640380 * exp(-5107.4 / model_sst[j][i]) 
					/ 1.22) * vit[0][0][j][i];

			// ! temperature to potential temperature
			theta[j][i] = tsa3[j*IMT + i] * pow(100000.0 / psa3[j*IMT + i], 0.286)
					* vit[0][0][j][i];
			runoff[0][j][i] = runoff3[j*IMT + i] * vit[0][0][j][i];
			seaice[0][j][i] = seaice3[j*IMT + i] * vit[0][0][j][i];
		}
	}

	double core_sensible[JMT][IMT];
	double core_latent[JMT][IMT];
	double core_tau[JMT][IMT];

	ncar_ocean_fluxes_jra(wspd3, &(theta[0][0]), &(model_sst[0][0]), 
			qar3, &(qs[0][0]), &(zz[0][0]), &(vit[0][0][0][0]), 
					&(core_sensible[0][0]), &(core_latent[0][0]), &(core_tau[0][0]), 
							&(ustar[0][0][0]));

	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			const double vit_times_one_minus_seaice = vit[0][0][j][i] * (1.0 - seaice[0][j][i]);

			sshf[0][j][i] = core_sensible[j][i] * vit_times_one_minus_seaice;
			lthf[0][j][i] = (core_latent[j][i] - snow3[j*IMT + i] * 3.335e5) 
					* vit_times_one_minus_seaice;

			double tmp_sst = model_sst[j][i] * model_sst[j][i];
			tmp_sst *= tmp_sst;
			lwv[0][j][i] = (0.95*lwv3[j*IMT + i] - 0.95*5.67e-8*tmp_sst)
					* vit_times_one_minus_seaice;

			qs[j][i] =   core_tau[j][i] * windx[j][i] * vit[0][0][j][i];
			zz[j][i] = - core_tau[j][i] * windy[j][i] * vit[0][0][j][i];

			nswv[0][j][i] = lwv[0][j][i] + sshf[0][j][i] + lthf[0][j][i];
			swv[0][j][i] = 0.934 * swv3[j*IMT + i] * vit_times_one_minus_seaice;

			fresh[0][j][i] = - (core_latent[j][i] / (2.5e+6) 
					+ rain3[j*IMT + i] + snow3[j*IMT + i] + runoff3[j*IMT + i])
							* vit_times_one_minus_seaice;

			ustar[0][j][i] = ustar[0][j][i] * vit[0][0][j][i];
		}
	}
	// ! tau to U/V grid
	for (int j = JST-1; j < JEM; ++j) {
		for (int i = 1; i < IMM; ++i) {
			su[0][j][i] = 0.25 * (
					qs[j  ][i] + qs[j  ][i-1] 
				+ qs[j+1][i] + qs[j+1][i-1]) * viv[0][0][j][i];
			sv[0][j][i] = 0.25 * (
					zz[j  ][i] + zz[j  ][i-1] 
				+ zz[j+1][i] + zz[j+1][i-1]) * viv[0][0][j][i];
		}
	}
}

void cpp_jra_daily_high(const int &iday) {
	using CppDynMod::u;
	using CppDynMod::v;


	using CppForcMod::buffer;

	using CppForcMod::tsa3;
	using CppForcMod::wspdu3;
	using CppForcMod::wspdv3;
	using CppForcMod::psa3;
	using CppForcMod::qar3;
	using CppForcMod::swv3;
	using CppForcMod::lwv3;
	using CppForcMod::rain3;
	using CppForcMod::snow3;
	using CppForcMod::runoff3;
	using CppForcMod::seaice3;
	using CppForcMod::wspd3;

	using CppForcMod::su;
	using CppForcMod::sv;
	using CppForcMod::lwv;
	using CppForcMod::swv;
	using CppForcMod::nswv;
	using CppForcMod::fresh;
	using CppForcMod::sshf;
	using CppForcMod::lthf;
	using CppForcMod::ustar;
	using CppForcMod::runoff;
	using CppForcMod::seaice;
	using CppParamMod::IMT;
	using CppParamMod::JMT;
	using CppParamMod::IMM;
	using CppParamMod::JST;
	using CppParamMod::JSM;
	using CppParamMod::JEM;
	using CppParamMod::S_IMT;
	using CppParamMod::S_JMT;
	using CppParamMod::IMT_GLOBAL;
	using CppParamMod::JMT_GLOBAL;
	using CppParamMod::MAX_BLOCKS_CLINIC;
	using CppPconstMod::vit;
	using CppPconstMod::viv;
	using CppPconstMod::s_lon;
	using CppPconstMod::s_lat;
	using CppMsgMod::nproc;
	using CppTracerMod::at;
	using CppPOPHaloMod::pop_halo_update_2dr8;

	int irec = iday;

	if (CppParamMod::mytid == 0) {
		printf("jra daily, irec = %d\n", irec);
	}

	// ! interplate to T grid
	interplationNearestThenScatterGlobal(0, irec, s_lon, s_lat, 
			CppForcMod::rf, runoff3, buffer);
	interplationNearestThenScatterGlobal(1, irec, s_lon, s_lat, 
			CppForcMod::t10, tsa3, buffer);
	interplationNearestThenScatterGlobal(2,  irec, s_lon, s_lat, 
			CppForcMod::u10, wspdu3, buffer);
	interplationNearestThenScatterGlobal(3,  irec, s_lon, s_lat, 
			CppForcMod::v10, wspdv3, buffer);
	interplationNearestThenScatterGlobal(4,  irec, s_lon, s_lat, 
			CppForcMod::slp, psa3, buffer);
	interplationNearestThenScatterGlobal(5,  irec, s_lon, s_lat, 
			CppForcMod::q10, qar3, buffer);
	interplationNearestThenScatterGlobal(6,  irec, s_lon, s_lat, 
			CppForcMod::swhf, swv3, buffer);
	interplationNearestThenScatterGlobal(7,  irec, s_lon, s_lat, 
			CppForcMod::lwhf, lwv3, buffer);
	interplationNearestThenScatterGlobal(8,  irec, s_lon, s_lat, 
			CppForcMod::precr, rain3, buffer);
	interplationNearestThenScatterGlobal(9,  irec, s_lon, s_lat, 
			CppForcMod::precs, snow3, buffer);
	interplationNearestThenScatterGlobal(10, irec, s_lon, s_lat, 
			CppForcMod::si, seaice3, buffer);
	
	double uu[JMT][IMT];
	double vv[JMT][IMT];
	const double epsln = 1e-25;
	for (int j = JSM - 1; j < JEM; ++j) {
		for (int i = 1; i < IMM; ++i) {
			const double tmp = vit[0][0][j][i] / (
					viv[0][0][j  ][i] + viv[0][0][j  ][i+1]
				+ viv[0][0][j-1][i] + viv[0][0][j-1][i+1] + epsln);
			uu[j][i] = tmp * (
					u[0][0][j  ][i] + u[0][0][j  ][i+1] 
				+ u[0][0][j-1][i] + u[0][0][j-1][i+1]);
			vv[j][i] = tmp * (
					v[0][0][j  ][i] + v[0][0][j  ][i+1] 
				+ v[0][0][j-1][i] + v[0][0][j-1][i+1]);
		}
	}
	double windx[JMT][IMT];
	double windy[JMT][IMT];
	double model_sst[JMT][IMT];
	double zz[JMT][IMT];
	double qs[JMT][IMT];
	double theta[JMT][IMT];
	const double tok = 273.15;
	// ! transfer core data to what the subroutine need
	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			// ! relative speed to surface currents
			windx[j][i] = (wspdu3[j*IMT + i] - uu[j][i]) * vit[0][0][j][i];
			windy[j][i] = (wspdv3[j*IMT + i] + vv[j][i]) * vit[0][0][j][i];

			// !  1.0 is from mom4
			wspd3[j*IMT + i] = sqrt(windx[j][i] * windx[j][i] 
													+ windy[j][i] * windy[j][i] + 1.0) * vit[0][0][j][i];

			// ! using a transient temperature, not daily mean
			model_sst[j][i] = (at[0][0][0][j][i] + tok) * vit[0][0][j][i];
			zz[j][i] = 10.0;
			qs[j][i] = (0.98 * 640380 * exp(-5107.4 / model_sst[j][i]) 
					/ 1.22) * vit[0][0][j][i];

			// ! temperature to potential temperature
			theta[j][i] = tsa3[j*IMT + i] * pow(100000.0 / psa3[j*IMT + i], 0.286)
					* vit[0][0][j][i];
			runoff[0][j][i] = runoff3[j*IMT + i] * vit[0][0][j][i];
			seaice[0][j][i] = seaice3[j*IMT + i] * vit[0][0][j][i];
		}
	}

	double core_sensible[JMT][IMT];
	double core_latent[JMT][IMT];
	double core_tau[JMT][IMT];

	ncar_ocean_fluxes_jra(wspd3, &(theta[0][0]), &(model_sst[0][0]), 
			qar3, &(qs[0][0]), &(zz[0][0]), &(vit[0][0][0][0]), 
					&(core_sensible[0][0]), &(core_latent[0][0]), &(core_tau[0][0]), 
							&(ustar[0][0][0]));

	for (int j = 0; j < JMT; ++j) {
		for (int i = 0; i < IMT; ++i) {
			const double vit_times_one_minus_seaice = vit[0][0][j][i] * (1.0 - seaice[0][j][i]);

			sshf[0][j][i] = core_sensible[j][i] * vit_times_one_minus_seaice;
			lthf[0][j][i] = (core_latent[j][i] - snow3[j*IMT + i] * 3.335e5) 
					* vit_times_one_minus_seaice;

			double tmp_sst = model_sst[j][i] * model_sst[j][i];
			tmp_sst *= tmp_sst;
			lwv[0][j][i] = (0.95*lwv3[j*IMT + i] - 0.95*5.67e-8*tmp_sst)
					* vit_times_one_minus_seaice;

			qs[j][i] =   core_tau[j][i] * windx[j][i] * vit[0][0][j][i];
			zz[j][i] = - core_tau[j][i] * windy[j][i] * vit[0][0][j][i];

			nswv[0][j][i] = lwv[0][j][i] + sshf[0][j][i] + lthf[0][j][i];
			swv[0][j][i] = 0.934 * swv3[j*IMT + i] * vit_times_one_minus_seaice;

			fresh[0][j][i] = - (core_latent[j][i] / (2.5e+6) 
					+ rain3[j*IMT + i] + snow3[j*IMT + i] + runoff3[j*IMT + i])
							* vit_times_one_minus_seaice;

			ustar[0][j][i] = ustar[0][j][i] * vit[0][0][j][i];
		}
	}
	// ! tau to U/V grid
	for (int j = JST-1; j < JEM; ++j) {
		for (int i = 1; i < IMM; ++i) {
			su[0][j][i] = 0.25 * (
					qs[j  ][i] + qs[j  ][i-1] 
				+ qs[j+1][i] + qs[j+1][i-1]) * viv[0][0][j][i];
			sv[0][j][i] = 0.25 * (
					zz[j  ][i] + zz[j  ][i-1] 
				+ zz[j+1][i] + zz[j+1][i-1]) * viv[0][0][j][i];
		}
	}
}

#endif // LICOM_ENABLE_FORTRAN
